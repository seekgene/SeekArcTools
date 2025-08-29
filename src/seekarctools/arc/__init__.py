import os
import sys
import json
import click
from ..utils.chemistry import CHEMISTRY
from ..utils.helper import logger, include_introns_callback, check_path
from ..utils.wrappers import cmd_execute
from ..utils.helper import check_all


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
@click.option("--chemistry", type=click.Choice(["DD-AG", "DD5-AG", "custom"]), help="DD-AG, DD5-AG.")
@click.pass_obj
def estep1(obj, **kwargs):
    # if kwargs["barcode"]:
    #     kwargs["chemistry"] = "custom"
    # else:
    #     kwargs["chemistry"] = "DD-Q"
    # if kwargs["sc5p"]:
    #     kwargs["chemistry"] = "DD_5G"
    os.makedirs(kwargs["outdir"], exist_ok=True)
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
@click.option("--sc5p",is_flag=True,default=False,show_default=True,help="If set, the single cell data is considered as 5' data.")
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
@click.pass_obj
def estep3(obj, **kwargs):
    from .estep3 import count, cell_calling
    count(**kwargs)

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
@click.option("--chemistry", type=click.Choice(["DD-AA", "custom"]), help="DD-AA.")
@click.pass_obj
def astep1(obj, **kwargs):
    # if kwargs["barcode"]:
    #     print(f'atac barcode path:{kwargs["barcode"]}')
    #     kwargs["chemistry"] = "custom"
    # else:
    #     kwargs["chemistry"] = "DD_AG"
    os.makedirs(kwargs["outdir"], exist_ok=True)
    from ..utils.atacbarcode import check_atac_options
    chemistry_kwargs = check_atac_options(**kwargs)
    kwargs.update(chemistry_kwargs)
    with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
        json.dump(kwargs, fh, indent=4)
    from ..utils.atacbarcode import barcode_main
    barcode_main(**kwargs)

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
@click.option('--gexoutdir', required=True,help='gex. step3/filtered_feature_bc_matrix/barcodes.tsv.gz')
@click.option("--refpath", "refpath", required=True, type=click.Path(exists=True, resolve_path=True), help="Path to reference.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--qvalue", default=0.05, show_default=True, help="Minimum FDR (q-value) cutoff for peak detection.")
@click.option("--nolambda", is_flag=True, default=False, show_default=True, help="If True, macs3 will use the background lambda as local lambda. This means macs3 will not consider the local bias at peak candidate regions.")
@click.option("--snapshift", default=0, show_default=True, help="The shift size in MACS.")
@click.option("--extsize", default=400, show_default=True, help="The extension size in MACS.")
@click.option("--min_len", default=400, show_default=True, help="The minimum length of a called peak. If None, it is set to 'extsize'.")
@click.option("--broad", is_flag=True, default=False, show_default=True, help="This option facilitates broad peak calling, producing results in the UCSC gappedPeak format which encapsulates a nested structure of peaks.")
@click.option("--broad_cutoff", default=0.1, show_default=True, help="Cutoff for the broad region.")
@click.option("--min_atac_count", default=None, show_default=True, help="Cell caller override: define the minimum number of ATAC transposition events in peaks (ATAC counts) for a cell barcode. Need to be used together with --min_gex_comnt.")
@click.option("--min_gex_count", default=None, show_default=True, help="Cell caller override: define the minimum number of GEX UMI counts for a cell barcode. Need to be used together with --min_atac_count.")
@click.option("--retry", is_flag=True, default=False, show_default=True, hidden=True, help="Skip the alignment step and rerun astep3.")
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
@click.option("--refpath", "refpath", required=True, type=click.Path(exists=True, resolve_path=True), help="Path to reference.")
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
@click.option("--chemistry", type=click.Choice(["DD-AG", "DD5-AG", "custom"]), 
              help="DD-AG, DD5-AG.")
@click.option("--core", default=4, show_default=True,
              help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--include-introns", "region", is_flag=True, default=False, callback=include_introns_callback, show_default=True,
              help="include introns or not.")
@click.option("--refpath", "refpath", required=True, type=click.Path(exists=True, resolve_path=True),
              help="Path to reference.")
@click.option("--star_path", "star_path", default="STAR",
              help="Path to executable STAR aligner.")
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
@click.option("--broad", is_flag=True, default=False, show_default=True, 
              help="This option facilitates broad peak calling, producing results in the UCSC gappedPeak format which encapsulates a nested structure of peaks.")
@click.option("--broad_cutoff", default=0.1, show_default=True, 
              help="Cutoff for the broad region.")
@click.option("--min_atac_count", default=None, show_default=True, 
              help="Cell caller override: define the minimum number of ATAC transposition events in peaks (ATAC counts) for a cell barcode.")
@click.option("--min_gex_count", default=None, show_default=True, 
              help="Cell caller override: define the minimum number of GEX UMI counts for a cell barcode.")
@click.option("--retry", is_flag=True, default=False, show_default=True, hidden=True, 
              help="Skip the alignment step and rerun astep3.")

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
    sampleoutdir = kwargs["outdir"]
    # check_all()

    # GEX
    kwargs["gexname"] = f'{kwargs["samplename"]}_E'
    kwargs["outdir"] = os.path.join(sampleoutdir, 'analysis', kwargs["gexname"])
    os.makedirs(kwargs["outdir"], exist_ok=True)

    if "estep1" in obj["steps"]:
        # if kwargs["barcode"]:
        #     kwargs["chemistry"] = "custom"
        # else:
        #     kwargs["chemistry"] = "DD-Q"
        # if kwargs["sc5p"]:
        #     kwargs["chemistry"] = "DD_5G"
        from ..utils.barcode import check_rna_options
        chemistry_kwargs = check_rna_options(fq1=kwargs["rnafq1"], fq2=kwargs["rnafq2"], **kwargs)
        kwargs.update(chemistry_kwargs)
        with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
            json.dump(kwargs, fh, indent=4)
        from ..utils.barcode import barcode_main
        barcode_main(fq1=kwargs["rnafq1"], fq2=kwargs["rnafq2"], **kwargs)

    # estep1_fq1 = os.path.join(kwargs["outdir"], "step1", f"{kwargs['gexname']}_1.fq.gz")
    estep1_fq2 = os.path.join(kwargs["outdir"], "step1", f"{kwargs['gexname']}_2.fq.gz")
    # paired: [r1.fq.gz, r2.fq.gz]
    kwargs["fq"] = [estep1_fq2, ]

    if "estep2" in obj["steps"]:
        from .estep2 import align
        align(stpes=_steps["estep2"], **kwargs)
    gexbam = os.path.join(kwargs["outdir"], "step2", "featureCounts",  f"{kwargs['gexname']}_SortedByName.bam")
   
    if "estep3" in obj["steps"]:
        from .estep3 import count
        count(bam=gexbam, **kwargs)

    raw_matrix = os.path.join(kwargs["outdir"], "step3", "raw_feature_bc_matrix")
    kwargs["raw_matrix"] = raw_matrix
       
    matrix = os.path.join(kwargs["outdir"], "step3", "filtered_feature_bc_matrix")
    kwargs["matrix"] = matrix
    kwargs["gexjson"] = os.path.join(kwargs["outdir"], f"{kwargs['gexname']}_summary.json")

    # ATAC
    kwargs["atacname"] = f'{kwargs["samplename"]}_A'
    kwargs["outdir"] = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"])
    os.makedirs(kwargs["outdir"], exist_ok=True)

    if "astep1" in obj["steps"]:
        kwargs["chemistry"] = "DD-AA"
        from ..utils.atacbarcode import check_atac_options
        chemistry_kwargs = check_atac_options(fq1=kwargs["atacfq1"], fq2=kwargs["atacfq2"], **kwargs)
        kwargs.update(chemistry_kwargs)
        with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
            json.dump(kwargs, fh, indent=4)
        from ..utils.atacbarcode import barcode_main
        barcode_main(fq1=kwargs["atacfq1"], fq2=kwargs["atacfq2"],**kwargs)
        kwargs["atacjson"] = os.path.join(kwargs["outdir"], f"{kwargs['atacname']}_summary.json")

    step1afq1 = os.path.join(kwargs["outdir"], "step1", f"{kwargs['atacname']}_1.fq.gz")
    step1afq2 = os.path.join(kwargs["outdir"], "step1", f"{kwargs['atacname']}_2.fq.gz")
    kwargs['afq'] = [step1afq1, step1afq2]

    if "astep2" in obj["steps"]:
        from .astep2 import align
        align(**kwargs)

    atacbam = os.path.join(kwargs["outdir"], "step2", "bwa_pe",  f"{kwargs['atacname']}_mem_pe_Sort.bam")
    kwargs["gexoutdir"] = os.path.join(sampleoutdir, 'analysis', kwargs["gexname"])

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
    gex_bam = os.path.join(sampleoutdir, 'analysis', kwargs["gexname"], "step2", "featureCounts", f"{kwargs['gexname']}_SortedByName.bam")
    atac_bam = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step2", "bwa_pe", f"{kwargs['atacname']}_mem_pe_Sort.bam")
    astep3Dir = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3")
    peakfile = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3", f"{kwargs['atacname']}_peaks.bed")
    rds_file = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step4", f"{kwargs['samplename']}.rds")
    cmd1 = f"mv {gex_bam} {kwargs['outdir']}"
    cmd_execute(cmd1, check=True)
    cmd2 = f"mv {atac_bam} {kwargs['outdir']}"
    cmd_execute(cmd2, check=True)
    cmd3 = f"cd {astep3Dir}; mv {kwargs['atacname']}_fragments.tsv.gz {kwargs['atacname']}_fragments.tsv.gz.tbi {kwargs['outdir']}"
    cmd_execute(cmd3, check=True)
    cmd4 = f"cp -rf {kwargs['gex_matrix']} {kwargs['outdir']}; cp -rf {kwargs['raw_matrix']} {kwargs['outdir']}"
    cmd_execute(cmd4, check=True)
    cmd5 = f"cp -rf {kwargs['atac_matrix']} {kwargs['outdir']} ; cp -rf {atac_rawmatrix} {kwargs['outdir']}"
    cmd_execute(cmd5, check=True)
    cmd6 = f"cp {peakfile} {kwargs['outdir']}"
    cmd_execute(cmd6, check=True)
    cmd7 = f"mv {rds_file} {kwargs['outdir']}"
    cmd_execute(cmd7, check=True)

@arc.command(help="Skip the alignment step and rerun astep3. Retry call peaks or cells.")
@click.pass_obj
@click.option("--samplename", required=True,
              help="Sample name.")
@click.option("--rawname", default="rawname", hidden=True,
              help="raw name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(resolve_path=True),
              help="Output dir.")
@click.option("--core", default=4, show_default=True,
              help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--refpath", "refpath", required=True, type=click.Path(exists=True, resolve_path=True),
              help="Path to reference.")
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
@click.option("--broad", is_flag=True, default=False, show_default=True, 
              help="This option facilitates broad peak calling, producing results in the UCSC gappedPeak format which encapsulates a nested structure of peaks.")
@click.option("--broad_cutoff", default=0.1, show_default=True, 
              help="Cutoff for the broad region.")
@click.option("--min_atac_count", default=None, show_default=True, 
              help="Cell caller override: define the minimum number of ATAC transposition events in peaks (ATAC counts) for a cell barcode.")
@click.option("--min_gex_count", default=None, show_default=True, 
              help="Cell caller override: define the minimum number of GEX UMI counts for a cell barcode.")

def retry(obj, **kwargs):
    kwargs["retry"] = True
    logger.info("Check the genomeDir path...")
    kwargs["genomeDir"] = os.path.join(kwargs["refpath"], "star")
    logger.info(check_path(kwargs["genomeDir"]))
    logger.info("Check the gtf path...")
    kwargs["gtf"] = os.path.join(kwargs["refpath"], "genes", "genes.gtf")
    logger.info(check_path(kwargs["gtf"]))
    logger.info("Check the genomefa path...")
    kwargs["genomefa"] = os.path.join(kwargs["refpath"], "fasta", "genome.fa")
    logger.info(check_path(kwargs["genomefa"]))
    sampleoutdir = kwargs["outdir"]
    kwargs["gexname"] = f'{kwargs["samplename"]}_E'
    kwargs["atacname"] = f'{kwargs["samplename"]}_A'

    if kwargs["min_atac_count"] != None and kwargs["min_gex_count"] != None:
        logger.info(f'retry force call cell with --min_atac_count:{kwargs["min_atac_count"]} and --min_gex_count:{kwargs["min_gex_count"]}')

    # astep3
    kwargs["outdir"] = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"])
    atacbam = os.path.join(kwargs["outdir"], "step2", "bwa_pe",  f"{kwargs['atacname']}_mem_pe_Sort.bam")
    kwargs["gexoutdir"] = os.path.join(sampleoutdir, 'analysis', kwargs["gexname"])
    from .astep3 import runpipe
    runpipe(bam=atacbam, **kwargs)

    # astep4
    kwargs["gex_matrix"] = os.path.join(kwargs["gexoutdir"], "step3", "filtered_feature_bc_matrix")
    gex_rawmatrix = os.path.join(kwargs["gexoutdir"], "step3", "raw_feature_bc_matrix")
    kwargs["atac_matrix"] = os.path.join(kwargs["outdir"], "step3", "filtered_peaks_bc_matrix")
    atac_rawmatrix = os.path.join(kwargs["outdir"], "step3", "raw_peaks_bc_matrix")
    kwargs["fragpath"] = os.path.join(sampleoutdir, "outs", f"{kwargs['atacname']}_fragments.tsv.gz")
    kwargs["gexjson"] = os.path.join(kwargs["gexoutdir"], f"{kwargs['gexname']}_summary.json")
    kwargs["atacjson"] = os.path.join(kwargs["outdir"], f"{kwargs['atacname']}_summary.json")

    from .astep4 import do_signac
    do_signac(**kwargs)

    # report
    kwargs["outdir"] = os.path.join(sampleoutdir, 'outs')
    from .report_arc import report
    report(**kwargs)

    # summary file
    astep3Dir = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3")
    peakfile = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step3", f"{kwargs['atacname']}_peaks.bed")
    rds_file = os.path.join(sampleoutdir, 'analysis', kwargs["atacname"], "step4", f"{kwargs['samplename']}.rds")
    
    cmd4 = f"cp -rf {kwargs['gex_matrix']} {kwargs['outdir']}; cp -rf {gex_rawmatrix} {kwargs['outdir']}"
    cmd_execute(cmd4, check=True)
    cmd5 = f"cp -rf {kwargs['atac_matrix']} {kwargs['outdir']} ; cp -rf {atac_rawmatrix} {kwargs['outdir']}"
    cmd_execute(cmd5, check=True)
    cmd6 = f"cp {peakfile} {kwargs['outdir']}"
    cmd_execute(cmd6, check=True)
    cmd7 = f"mv {rds_file} {kwargs['outdir']}"
    cmd_execute(cmd7, check=True)