import os
import sys
import json
import click
from ..utils.chemistry import CHEMISTRY
from ..utils.helper import logger, include_introns_callback, check_path
from ..utils.wrappers import cmd_execute


_steps = {
    "estep1": [],
    "estep2": ["STAR", "SortByPos", "FeatureCounts", "SortByName"],
    "estep3": [],
    "astep1": [],
    "astep2": [],
    "astep3": [],
    "astep4": [],
}

@click.group(help="pipeline can only analyze Gene Expression and ATAC data together.")
@click.option("--steps", default=None, type=click.Path(), help="json format.")
@click.pass_obj
def arc(obj, steps):
    if  steps:
        with open(steps) as fh:
            obj["steps"] = json.load(fh)
    else:
        obj["steps"] = _steps



@arc.command(help="Gene Expression extract barcode and umi.")
@click.option("--fq1", "fq1", required=True, type=click.Path(), multiple=True, help="Expression Read1 fq file, can specify multiple times.")
@click.option("--fq2", "fq2", required=True, type=click.Path(), multiple=True, help="Expression Read2 fq file, can specify multiple times.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--shift", is_flag=True, default=False, show_default=True, help="Shift, used to describe read1 structure.")
@click.option("--pattern", "shift_pattern", default="A", help="Anchor sequence, used to describe read1 structure.")
@click.option("--barcode", multiple=True, help="Barcode white list file, can specify multiple times.")
@click.option("--structure", help="Used to describe read1 structure.")
@click.option("--linker", multiple=True, help="Linker white list file, can specify multiple times.")
@click.option("--skip_misB", "do_B_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of barcode.")
@click.option("--skip_misL", "do_L_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of linker.")
@click.option("--skip_multi", "use_multi", is_flag=True, default=True, show_default=True, help="Do not rescue barcode match multi when do correction.")
@click.option("--skip_len", "use_short_read", is_flag=True, default=False, show_default=True, help="Skip filtering short reads after adapter filter, short reads will be used.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
# @click.option("--chemistry", type=click.Choice(["DDV1", "DDV2", "DDVS", "DD5V1", "MM", "MM-D", "DD-Q", "custom"]), help="DDV1, DDV2, DDVS, DD5V1, MM, MM-D, DD-Q.")
@click.pass_obj
def estep1(obj, **kwargs):
    if kwargs["barcode"]:
        print(f'atac barcode path:{kwargs["barcode"]}')
        kwargs["chemistry"] = "custom"
    else:
        kwargs["chemistry"] = "DD-Q"
    # kwargs["chemistry"] = "DD-Q"
    from ..utils.barcode import check_rna_options
    chemistry_kwargs = check_rna_options(**kwargs)
    kwargs.update(chemistry_kwargs)
    with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
        json.dump(kwargs, fh, indent=4)
    from ..utils.barcode import barcode_main
    barcode_main(**kwargs)

@arc.command(help="Gene Expression align reads to genome.")
@click.option("--fq", required=True, multiple=True, help="Expression Step1 Read2 fq file")
@click.option("--genomeDir", "genomeDir", required=True, type=click.Path(), help="Path to dir which store the genome indices.")
@click.option("--gtf", required=True, type=click.Path(), help="Path to GTF file.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--star_path", "star_path", default="STAR", help="Path to executable STAR aligner.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--include-introns", "region", is_flag=True, default=False, callback=include_introns_callback, show_default=True, help="include introns or not.")
@click.pass_obj
def estep2(obj, **kwargs):
    from .estep2 import align
    align(**kwargs)

@arc.command(help="Gene Expression quantifies.")
@click.option("--bam", required=True, help="Bam file which contain alignment info.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--gtf", required=True, type=click.Path(), help="Path to GTF file.")
@click.option("--umi_correct_method", type=click.Choice(["cluster", "adjacency", "directional"]), default="adjacency", show_default=True, help="cluster, adjacency, directional")
@click.option("--expectNum", "expectNum", default=3000, show_default=True, help="Expected number of cells that used as input in cell calling algorithm.")
@click.option("--forceCell", "forceCell", help="Force pipeline to use this number of cells, skipping cell calling algorithm.",required=False)
@click.pass_obj
def estep3(obj, **kwargs):
    from .estep3 import count, cell_calling
    count(**kwargs)
    raw_matrix = os.path.join(kwargs["outdir"], "step3", "raw_feature_bc_matrix")
    kwargs["raw_matrix"] = raw_matrix
    cell_calling(**kwargs)

@arc.command(help="cell calling")
@click.option("--raw_matrix", required=True, help="Raw Feature-barcode matrix.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--expectNum", "expectNum", default=3000, show_default=True, help="Expected number of cells that used as input in cell calling algorithm.")
@click.option("--forceCell", "forceCell", help="Force pipeline to use this number of cells, skipping cell calling algorithm.",required=False)
@click.pass_obj
def callcell(obj, **kwargs):
    from .estep3 import cell_calling
    cell_calling(**kwargs)

@arc.command(help="ATAC extract cell barcode and umi. cut reads to preserve valid sequences.")
@click.option("--fq1", "fq1", required=True, type=click.Path(), multiple=True, help="ATAC Read1 fq file, can specify multiple times.")
@click.option("--fq2", "fq2", required=True, type=click.Path(), multiple=True, help="ATAC Read2 fq file, can specify multiple times.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--shift", is_flag=True, default=False, show_default=True, help="Shift, used to describe read1 structure.")
@click.option("--pattern", "shift_pattern", default="A", help="Anchor sequence, used to describe read1 structure.")
@click.option("--barcode", multiple=True, help="Barcode white list file, can specify multiple times.")
@click.option("--structure", help="Used to describe read1 structure.")
@click.option("--linker", multiple=True, help="Linker white list file, can specify multiple times.")
@click.option("--skip_misB", "do_B_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of barcode.")
@click.option("--skip_misL", "do_L_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of linker.")
@click.option("--skip_multi", "use_multi", is_flag=True, default=True, show_default=True, help="Do not rescue barcode match multi when do correction.")
@click.option("--skip_len", "use_short_read", is_flag=True, default=False, show_default=True, help="Skip filtering short reads after adapter filter, short reads will be used.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
# @click.option("--chemistry", type=click.Choice(["DDV1", "DDV2", "DDVS", "DD5V1", "MM", "MM-D", "DD-Q", "custom"]), help="DDV1, DDV2, DDVS, DD5V1, MM, MM-D, DD-Q.")
@click.pass_obj
def astep1(obj, **kwargs):
    if kwargs["barcode"]:
        print(f'atac barcode path:{kwargs["barcode"]}')
        kwargs["chemistry"] = "custom"
    else:
        kwargs["chemistry"] = "DD_AG"
    # kwargs["chemistry"] = "DD_AG"
    from ..utils.atacbarcode import check_atac_options
    chemistry_kwargs = check_atac_options(**kwargs)
    kwargs.update(chemistry_kwargs)
    with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
        json.dump(kwargs, fh, indent=4)
    from ..utils.atacbarcode import barcode_main
    barcode_main(**kwargs)
    kwargs["atacname"] = f'{kwargs["samplename"]}_A'
    kwargs["atacjson"] = os.path.join(kwargs["outdir"], f"{kwargs['atacname']}_summary.json")
    kwargs["step1afq1"] = os.path.join(kwargs["outdir"], "step1", f"{kwargs['atacname']}_1.fq.gz")
    kwargs["step1afq2"] = os.path.join(kwargs["outdir"], "step1", f"{kwargs['atacname']}_2.fq.gz")
    from ..utils.countUtil import cutreads
    cutreads(**kwargs)

@arc.command(help="ATAC align reads to genome.")
@click.option("--afq", required=True, multiple=True, help="ATAC Step1 Read1 Read2 cut fq file")
@click.option("--genomefa", "genomefa", required=True, type=click.Path(), help="Path to genome fasta file.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.pass_obj
def astep2(obj, **kwargs):
    from .astep2 import align
    align(**kwargs)


@arc.command(help="ATAC quantifies.")
@click.option('--bam', required=True,help='atac. step2/bwa_pe/asample_mem_pe_Sort.bam')
@click.option("--outdir", default="./step3/", show_default=True, type=click.Path(), help="Output dir.")
@click.option('--samplename', required=True,help='Sample name.')
@click.option('--filtercb', required=True,help='gex. step3/filtered_feature_bc_matrix/barcodes.tsv.gz')
@click.option('--countxls', required=True,help='gex. step3/counts.xls')
@click.option('--organism', type=click.Choice(["human", "mouse"]), help="human or mouse.")
@click.option("--refpath", "refpath", required=True, type=click.Path(exists=True, resolve_path=True), help="Path to reference.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--qvalue", default=0.05, show_default=True, help="Minimum FDR (q-value) cutoff for peak detection.")
@click.option("--nolambda", is_flag=True, default=False, show_default=True, help="If True, macs3 will use the background lambda as local lambda. This means macs3 will not consider the local bias at peak candidate regions.")
@click.option("--snapshift", default=0, show_default=True, help="The shift size in MACS.")
@click.option("--extsize", default=400, show_default=True, help="The extension size in MACS.")
@click.option("--min_len", default=400, show_default=True, help="The minimum length of a called peak. If None, it is set to 'extsize'.")
@click.pass_obj
def astep3(obj, **kwargs):
    from .astep3 import runpipe
    runpipe(**kwargs)

@arc.command(help="ATAC do signac.")
@click.option("--gex_matrix", required=True, help="Filtered Feature-barcode matrix.")
@click.option("--atac_matrix", required=True, help="Filtered Peak-barcode matrix.")
@click.option("--fragpath", required=True, help="Path to fragments file.")
@click.option('--samplename', required=True,help='Sample name.')
@click.option("--rawname", default="rawname", hidden=True, help="raw name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--organism", type=click.Choice(["human", "mouse"]), help="human or mouse.")
@click.option("--anno_rds", required=True, help="Path to Anno_EnsDb.rds.")
@click.pass_obj
def astep4(obj, **kwargs):
    from .astep4 import do_signac
    do_signac(**kwargs)

@arc.command(help="report.")
@click.option('--gexjson', required=True,help="Path to gex sample json file.")
@click.option('--atacjson', required=True,help="Path to atac sample json file.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option('--samplename', required=True,help="Sample name.")
@click.option("--rawname", default="rawname", hidden=True, help="raw name.")
@click.pass_obj
def report(obj, **kwargs):
    from .report_arc import report
    report(**kwargs)


@arc.command(help="run all steps.")
@click.pass_obj
@click.option("--rnafq1", "rnafq1", required=True, type=click.Path(exists=True, resolve_path=True), multiple=True,
              help="Expression Read1 fq file, can specify multiple times.")
@click.option("--rnafq2", "rnafq2", required=True, type=click.Path(exists=True, resolve_path=True),
              multiple=True, help="Expression Read2 fq file, can specify multiple times.")
@click.option("--atacfq1", "atacfq1", required=True, type=click.Path(exists=True, resolve_path=True), multiple=True,
              help="ATAC Read1 fq file, can specify multiple times.")
@click.option("--atacfq2", "atacfq2", required=True, type=click.Path(exists=True, resolve_path=True),
              multiple=True, help="ATAC Read2 fq file, can specify multiple times.")
@click.option("--samplename", required=True,
              help="Sample name.")
@click.option("--rawname", default="rawname", hidden=True,
              help="raw name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(resolve_path=True),
              help="Output dir.")
@click.option("--shift", is_flag=True, default=False, hidden=True,
              help="Shift, used to describe read1 structure.")
@click.option("--pattern", "shift_pattern", default="A", hidden=True,
              help="Anchor sequence, used to describe read1 structure.")
@click.option("--barcode", multiple=True, hidden=True,
              help="Barcode white list file, can specify multiple times, used to describe read1 structure.")
@click.option("--structure", hidden=True,
              help="Used to describe read1 structure.")
@click.option("--linker", multiple=True, hidden=True,
              help="Linker white list file, can specify multiple times, used to describe read1 structure.")
@click.option("--skip_misB", "do_B_correction", is_flag=True, default=True, show_default=True,
              help="Not allow one base err correction in each part of barcode.")
@click.option("--skip_misL", "do_L_correction", is_flag=True, default=True, show_default=True,
              help="Not allow one base err correction in each part of linker.")
@click.option("--skip_multi", "use_multi", is_flag=True, default=True, show_default=True,
              help="Do not rescue barcode match multi when do correction.")
@click.option("--skip_len", "use_short_read", is_flag=True, default=False, show_default=True, 
              help="Skip filtering short reads after adapter filter, short reads will be used.")
@click.option("--core", default=4, show_default=True,
              help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--include-introns", "region", is_flag=True, default=False, callback=include_introns_callback, show_default=True,
              help="include introns or not.")
@click.option("--organism", type=click.Choice(["human", "mouse"]), 
              help="human or mouse") # 提示非模式
@click.option("--refpath", "refpath", required=True, type=click.Path(exists=True, resolve_path=True),
              help="Path to reference.")
@click.option("--star_path", "star_path", default="STAR",
              help="Path to executable STAR aligner.")
@click.option("--expectNum", "expectNum", default=3000, show_default=True,
              help="Expected number of cells that used as input in cell calling algorithm.")
@click.option("--forceCell", "forceCell", required=False,
              help="Force pipeline to use this number of cells, skipping cell calling algorithm.")
@click.option("--umi_correct_method", type=click.Choice(["cluster", "adjacency", "directional"]), default="adjacency", show_default=True, hidden=True,
              help="cluster, adjacency, directional")
@click.option("--qvalue", default=0.05, show_default=True, 
              help="Minimum FDR (q-value) cutoff for peak detection.")
@click.option("--nolambda", is_flag=True, default=False, show_default=True, 
              help="If True, macs3 will use the background lambda as local lambda. This means macs3 will not consider the local bias at peak candidate regions.")
@click.option("--snapshift", default=0, show_default=True, 
              help="The shift size in MACS.")
@click.option("--extsize", default=400, show_default=True, 
              help="The extension size in MACS.")
@click.option("--min_len", default=400, show_default=True, 
              help="The minimum length of a called peak. If None, it is set to 'extsize'.")

def run(obj, **kwargs):
    logger.info("Check the genomeDir path...")
    kwargs["genomeDir"] = os.path.join(kwargs["refpath"], "star")
    logger.info(check_path(kwargs["genomeDir"]))
    logger.info("Check the gtf path...")
    kwargs["gtf"] = os.path.join(kwargs["refpath"], "genes", "genes.gtf")
    logger.info(check_path(kwargs["gtf"]))
    logger.info("Check the genomefa path...")
    kwargs["genomefa"] = os.path.join(kwargs["refpath"], "fasta", "genome.fa")
    logger.info(check_path(kwargs["genomefa"]))
    logger.info("Check the anno_rds path...")
    kwargs["anno_rds"] = os.path.join(kwargs["refpath"], "anno", "Anno_EnsDb.rds")
    logger.info(check_path(kwargs["anno_rds"]))
    sampleoutdir = kwargs["outdir"]

    # GEX
    kwargs["gexname"] = f'{kwargs["samplename"]}_E'
    kwargs["outdir"] = os.path.join(sampleoutdir, 'analysis', kwargs["gexname"])
    os.makedirs(kwargs["outdir"], exist_ok=True)

    if "estep1" in obj["steps"]:
        kwargs["chemistry"] = "DD-Q"
        from ..utils.barcode import check_rna_options
        chemistry_kwargs = check_rna_options(fq1=kwargs["rnafq1"], fq2=kwargs["rnafq2"], **kwargs)
        kwargs.update(chemistry_kwargs)
        with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
            json.dump(kwargs, fh, indent=4)
        from ..utils.barcode import barcode_main
        barcode_main(fq1=kwargs["rnafq1"], fq2=kwargs["rnafq2"], **kwargs)

    fq = os.path.join(kwargs["outdir"], "step1", f"{kwargs['gexname']}_2.fq.gz")
    # paired: [r1.fq.gz, r2.fq.gz]
    kwargs["fq"] = [fq, ]

    if "estep2" in obj["steps"]:
        from .estep2 import align
        align(stpes=_steps["estep2"], **kwargs)
    gexbam = os.path.join(kwargs["outdir"], "step2", "featureCounts",  f"{kwargs['gexname']}_SortedByName.bam")
    # kwargs["bam"] = bam
   
    if "estep3" in obj["steps"]:
        from .estep3 import count
        count(bam=gexbam, **kwargs)

    raw_matrix = os.path.join(kwargs["outdir"], "step3", "raw_feature_bc_matrix")
    kwargs["raw_matrix"] = raw_matrix
    if "estep3" in obj["steps"]:
        from .estep3 import cell_calling
        cell_calling(**kwargs)
       
    matrix = os.path.join(kwargs["outdir"], "step3", "filtered_feature_bc_matrix")
    kwargs["matrix"] = matrix
    kwargs["gexjson"] = os.path.join(kwargs["outdir"], f"{kwargs['gexname']}_summary.json")

    # ATAC
    kwargs["atacname"] = f'{kwargs["samplename"]}_A'
    kwargs["outdir"] = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"])
    os.makedirs(kwargs["outdir"], exist_ok=True)

    if "astep1" in obj["steps"]:
        kwargs["chemistry"] = "DD_AG"
        from ..utils.atacbarcode import check_atac_options
        chemistry_kwargs = check_atac_options(fq1=kwargs["atacfq1"], fq2=kwargs["atacfq2"], **kwargs)
        kwargs.update(chemistry_kwargs)
        with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
            json.dump(kwargs, fh, indent=4)
        from ..utils.atacbarcode import barcode_main
        barcode_main(fq1=kwargs["atacfq1"], fq2=kwargs["atacfq2"],**kwargs)
        kwargs["atacjson"] = os.path.join(kwargs["outdir"], f"{kwargs['atacname']}_summary.json")
        kwargs["step1afq1"] = os.path.join(kwargs["outdir"], "step1", f"{kwargs['atacname']}_1.fq.gz")
        kwargs["step1afq2"] = os.path.join(kwargs["outdir"], "step1", f"{kwargs['atacname']}_2.fq.gz")
        from ..utils.countUtil import cutreads
        cutreads(**kwargs)

    cutafq1 = os.path.join(kwargs["outdir"], "step1", f"{kwargs['atacname']}_cutR1.fastq.gz")
    cutafq2 = os.path.join(kwargs["outdir"], "step1", f"{kwargs['atacname']}_cutR2.fastq.gz")
    kwargs['afq'] = [cutafq1, cutafq2]

    if "astep2" in obj["steps"]:
        from .astep2 import align
        align(**kwargs)

    atacbam = os.path.join(kwargs["outdir"], "step2", "bwa_pe",  f"{kwargs['atacname']}_mem_pe_Sort.bam")
    kwargs["filtercb"] = os.path.join(sampleoutdir, 'analysis', kwargs["gexname"], "step3", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    kwargs["countxls"] = os.path.join(sampleoutdir, 'analysis', kwargs["gexname"], "step3", "counts.xls")

    if "astep3" in obj["steps"]:
        from .astep3 import runpipe
        runpipe(bam=atacbam, **kwargs)

    kwargs["gex_matrix"] = matrix
    atacmatrix = os.path.join(kwargs["outdir"], "step3", "filtered_peaks_bc_matrix")
    kwargs["atac_matrix"] = atacmatrix
    atac_rawmatrix = os.path.join(kwargs["outdir"], "step3", "raw_peaks_bc_matrix")
    kwargs["fragpath"] = os.path.join(kwargs["outdir"], "step3", f"{kwargs['atacname']}_fragments.tsv.gz")

    if "astep4" in obj["steps"]:
        from .astep4 import do_signac
        do_signac(**kwargs)

    kwargs["outdir"] = os.path.join(sampleoutdir, 'outs')
    from .report_arc import report
    report(**kwargs)

    # summary file
    featureCountsDir = os.path.join(sampleoutdir, 'analysis', kwargs["gexname"], "step2", "featureCounts")
    bwaDir = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step2", "bwa_pe")
    astep3Dir = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3")
    peakfile = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3", f"{kwargs['atacname']}_peaks.bed")
    rds_file = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step4", f"{kwargs['samplename']}.rds")
    cmd1 = f"cd {featureCountsDir}; mv {kwargs['gexname']}_SortedByName.bam ../../../../outs/; ln -sf ../../../../outs/{kwargs['gexname']}_SortedByName.bam {kwargs['gexname']}_SortedByName.bam"
    cmd_execute(cmd1, check=True)
    cmd2 = f"cd {bwaDir}; mv {kwargs['atacname']}_mem_pe_Sort.bam ../../../../outs/; ln -sf ../../../../outs/{kwargs['atacname']}_mem_pe_Sort.bam {kwargs['atacname']}_mem_pe_Sort.bam"
    cmd_execute(cmd2, check=True)
    cmd3 = f"cd {astep3Dir}; mv {kwargs['atacname']}_fragments.tsv.gz {kwargs['atacname']}_fragments.tsv.gz.tbi ../../../outs/; ln -sf ../../../outs/{kwargs['atacname']}_fragments.tsv.gz {kwargs['atacname']}_fragments.tsv.gz; ln -sf ../../../outs/{kwargs['atacname']}_fragments.tsv.gz.tbi {kwargs['atacname']}_fragments.tsv.gz.tbi"
    cmd_execute(cmd3, check=True)
    cmd4 = f"cp -rf {kwargs['gex_matrix']} {kwargs['outdir']}; cp -rf {kwargs['raw_matrix']} {kwargs['outdir']}"
    cmd_execute(cmd4, check=True)
    cmd5 = f"cp -rf {kwargs['atac_matrix']} {kwargs['outdir']} ; cp -rf {atac_rawmatrix} {kwargs['outdir']}"
    cmd_execute(cmd5, check=True)
    cmd6 = f"cp {peakfile} {kwargs['outdir']} ; mv {rds_file} {kwargs['outdir']}"
    cmd_execute(cmd6, check=True)

@arc.command(help="retry some steps.")
@click.pass_obj
@click.option("--samplename", required=True,
              help="Sample name.")
@click.option("--rawname", default="rawname", hidden=True,
              help="raw name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(resolve_path=True),
              help="Output dir.")
@click.option("--core", default=4, show_default=True,
              help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--organism", type=click.Choice(["human", "mouse"]), 
              help="human or mouse") # 提示非模式
@click.option("--refpath", "refpath", required=True, type=click.Path(exists=True, resolve_path=True),
              help="Path to reference.")
@click.option("--forceCell", "forceCell", required=False,
              help="Force pipeline to use this number of cells, skipping cell calling algorithm.")
@click.option("--qvalue", default=0.05, show_default=True, 
              help="Minimum FDR (q-value) cutoff for peak detection.")
@click.option("--nolambda", is_flag=True, default=False, show_default=True, 
              help="If True, macs3 will use the background lambda as local lambda. This means macs3 will not consider the local bias at peak candidate regions.")
@click.option("--snapshift", default=0, show_default=True, 
              help="The shift size in MACS.")
@click.option("--extsize", default=400, show_default=True, 
              help="The extension size in MACS.")
@click.option("--min_len", default=400, show_default=True, 
              help="The minimum length of a called peak. If None, it is set to 'extsize'.")

def retry(obj, **kwargs):
    logger.info("Check the genomeDir path...")
    kwargs["genomeDir"] = os.path.join(kwargs["refpath"], "star")
    logger.info(check_path(kwargs["genomeDir"]))
    logger.info("Check the gtf path...")
    kwargs["gtf"] = os.path.join(kwargs["refpath"], "genes", "genes.gtf")
    logger.info(check_path(kwargs["gtf"]))
    logger.info("Check the genomefa path...")
    kwargs["genomefa"] = os.path.join(kwargs["refpath"], "fasta", "genome.fa")
    logger.info(check_path(kwargs["genomefa"]))
    logger.info("Check the anno_rds path...")
    kwargs["anno_rds"] = os.path.join(kwargs["refpath"], "anno", "Anno_EnsDb.rds")
    logger.info(check_path(kwargs["anno_rds"]))
    sampleoutdir = kwargs["outdir"]
    kwargs["gexname"] = f'{kwargs["samplename"]}_E'
    kwargs["atacname"] = f'{kwargs["samplename"]}_A'
    # fragpath=os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3", kwargs['atacname']) + "_fragments.tsv.gz"
    # fragindex=os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3", kwargs['atacname']) + "_fragments.tsv.gz.tbi"
    # if os.path.exists(fragpath):
    #     cmd7 = f"rm {fragpath}"
    #     cmd_execute(cmd7, check=True)
    # if os.path.exists(fragindex):
    #     cmd8 = f"rm {fragindex}"
    #     cmd_execute(cmd8, check=True)
    

    if kwargs["forceCell"] != None:
        print(f'retry forceCell : {kwargs["forceCell"]}')

        # GEX
        kwargs["outdir"] = os.path.join(sampleoutdir, 'analysis', kwargs["gexname"])
        raw_matrix = check_path(os.path.join(kwargs["outdir"], "step3", "raw_feature_bc_matrix"))
        kwargs["raw_matrix"] = raw_matrix
        from .estep3 import cell_calling
        cell_calling(**kwargs)
        
        kwargs["gex_matrix"] = os.path.join(kwargs["outdir"], "step3", "filtered_feature_bc_matrix")
        kwargs["gexjson"] = os.path.join(kwargs["outdir"], f"{kwargs['gexname']}_summary.json")

        # ATAC
        # astep3
        kwargs["outdir"] = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"])
        atacbam = check_path(os.path.join(kwargs["outdir"], "step2", "bwa_pe",  f"{kwargs['atacname']}_mem_pe_Sort.bam"))
        kwargs["filtercb"] = check_path(os.path.join(sampleoutdir, 'analysis', kwargs["gexname"], "step3", "filtered_feature_bc_matrix", "barcodes.tsv.gz"))
        kwargs["countxls"] = check_path(os.path.join(sampleoutdir, 'analysis', kwargs["gexname"], "step3", "counts.xls"))
        from .astep3 import runpipe
        runpipe(bam=atacbam, **kwargs)

        kwargs["atac_matrix"] = os.path.join(kwargs["outdir"], "step3", "filtered_peaks_bc_matrix")
        atac_rawmatrix = os.path.join(kwargs["outdir"], "step3", "raw_peaks_bc_matrix")
        kwargs["fragpath"] = os.path.join(kwargs["outdir"], "step3", f"{kwargs['atacname']}_fragments.tsv.gz")
        kwargs["atacjson"] = os.path.join(kwargs["outdir"], f"{kwargs['atacname']}_summary.json")

        # astep4
        from .astep4 import do_signac
        do_signac(**kwargs)

        kwargs["outdir"] = os.path.join(sampleoutdir, 'outs')
        from .report_arc import report
        report(**kwargs)

        # summary file
        astep3Dir = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3")
        peakfile = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3", f"{kwargs['atacname']}_peaks.bed")
        rds_file = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step4", f"{kwargs['samplename']}.rds")
        
        cmd3 = f"cd {astep3Dir}; mv {kwargs['atacname']}_fragments.tsv.gz {kwargs['atacname']}_fragments.tsv.gz.tbi ../../../outs/; ln -sf ../../../outs/{kwargs['atacname']}_fragments.tsv.gz {kwargs['atacname']}_fragments.tsv.gz; ln -sf ../../../outs/{kwargs['atacname']}_fragments.tsv.gz.tbi {kwargs['atacname']}_fragments.tsv.gz.tbi"
        cmd_execute(cmd3, check=True)
        cmd4 = f"cp -rf {kwargs['gex_matrix']} {kwargs['outdir']}; cp -rf {kwargs['raw_matrix']} {kwargs['outdir']}"
        cmd_execute(cmd4, check=True)
        cmd5 = f"cp -rf {kwargs['atac_matrix']} {kwargs['outdir']} ; cp -rf {atac_rawmatrix} {kwargs['outdir']}"
        cmd_execute(cmd5, check=True)
        cmd6 = f"cp {peakfile} {kwargs['outdir']} ; mv {rds_file} {kwargs['outdir']}"
        cmd_execute(cmd6, check=True)

    else:
        print(f'retry callpeak : Parameter : qvalue={kwargs["qvalue"]}, nolambda={kwargs["nolambda"]}, shift={kwargs["snapshift"]}, extsize={kwargs["extsize"]}, min_len={kwargs["min_len"]}, n_jobs={kwargs["core"]}')
        # astep3
        kwargs["outdir"] = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"])
        atacbam = check_path(os.path.join(kwargs["outdir"], "step2", "bwa_pe",  f"{kwargs['atacname']}_mem_pe_Sort.bam"))
        kwargs["filtercb"] = check_path(os.path.join(sampleoutdir, 'analysis', kwargs["gexname"], "step3", "filtered_feature_bc_matrix", "barcodes.tsv.gz"))
        kwargs["countxls"] = check_path(os.path.join(sampleoutdir, 'analysis', kwargs["gexname"], "step3", "counts.xls"))
        from .astep3 import runpipe
        runpipe(bam=atacbam, **kwargs)
        kwargs["gex_matrix"] = os.path.join(sampleoutdir, 'analysis', kwargs["gexname"], "step3", "filtered_feature_bc_matrix")
        kwargs["atac_matrix"] = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3", "filtered_peaks_bc_matrix")
        atac_rawmatrix = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3", "raw_peaks_bc_matrix")
        kwargs["fragpath"] = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3", f"{kwargs['atacname']}_fragments.tsv.gz")

        # astep4
        from .astep4 import do_signac
        do_signac(**kwargs)

        kwargs["outdir"] = os.path.join(sampleoutdir, 'outs')
        kwargs["gexjson"] = os.path.join(sampleoutdir, 'analysis', kwargs["gexname"], f"{kwargs['gexname']}_summary.json")
        kwargs["atacjson"] = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], f"{kwargs['atacname']}_summary.json")
        from .report_arc import report
        report(**kwargs)

        # summary file
        astep3Dir = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3")
        peakfile = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3", f"{kwargs['atacname']}_peaks.bed")
        rds_file = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step4", f"{kwargs['samplename']}.rds")

        cmd3 = f"cd {astep3Dir}; mv {kwargs['atacname']}_fragments.tsv.gz {kwargs['atacname']}_fragments.tsv.gz.tbi ../../../outs/; ln -sf ../../../outs/{kwargs['atacname']}_fragments.tsv.gz {kwargs['atacname']}_fragments.tsv.gz; ln -sf ../../../outs/{kwargs['atacname']}_fragments.tsv.gz.tbi {kwargs['atacname']}_fragments.tsv.gz.tbi"
        cmd_execute(cmd3, check=True)
        cmd5 = f"cp -rf {kwargs['atac_matrix']} {kwargs['outdir']} ; cp -rf {atac_rawmatrix} {kwargs['outdir']}"
        cmd_execute(cmd5, check=True)
        cmd6 = f"cp {peakfile} {kwargs['outdir']} ; mv {rds_file} {kwargs['outdir']}"
        cmd_execute(cmd6, check=True)

