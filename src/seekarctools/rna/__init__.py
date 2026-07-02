import json
import os

import click

from ..utils.helper import include_introns_callback, logger, check_path


@click.group(help="pipeline can analyze Gene Expression data only.")
def rna():
    pass


@rna.command(help="Gene Expression extract barcode and umi.")
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
@click.option("--chemistry", type=click.Choice(["DD_AG_RNA", "DD_AG", "DD5_AG", "custom"]), default="DD_AG_RNA", show_default=True, help="RNA chemistry.")
def step1(**kwargs):
    os.makedirs(kwargs["outdir"], exist_ok=True)

    from ..utils.barcode import check_rna_options, barcode_main

    chemistry_kwargs = check_rna_options(**kwargs)
    kwargs.update(chemistry_kwargs)

    with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
        json.dump(kwargs, fh, indent=4)

    barcode_main(**kwargs)


@rna.command(help="Gene Expression align reads to genome.")
@click.option("--fq", required=True, multiple=True, help="Expression Step1 Read2 fq file.")
@click.option("--genomeDir", "genomeDir", required=True, type=click.Path(), help="Path to dir which store the genome indices.")
@click.option("--gtf", required=True, type=click.Path(), help="Path to GTF file.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--star_path", "star_path", default="STAR", help="Path to executable STAR aligner.")
@click.option("--scoremin", default=None, type=float, help="STAR --outFilterScoreMinOverLread. Example: 0.33 for loose alignment.")
@click.option("--matchnmin", default=None, type=float, help="STAR --outFilterMatchNminOverLread. Example: 0.33 for loose alignment.")
@click.option("--pairendAlignment", "pairend_alignment", is_flag=True, default=False, show_default=True, help="Use paired-end STAR alignment with step1 R1 and R2.")
@click.option("--pairendMinR1Len", "pairend_min_r1_len", default=20, show_default=True, type=int, help="Minimum R1 length kept for paired-end STAR alignment.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--sc5p", is_flag=True, default=False, show_default=True, help="If set, the single cell data is considered as 5' data.")
@click.option("--include-introns", "region", is_flag=True, default=False, callback=include_introns_callback, show_default=True, help="include introns or not.")
def step2(**kwargs):
    from ..arc.estep2 import align

    align(**kwargs)


@rna.command(help="Gene Expression quantifies and calls cells.")
@click.option("--bam", required=True, help="Bam file which contain alignment info.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(), help="Output dir.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--gtf", required=True, type=click.Path(), help="Path to GTF file.")
@click.option("--umi_correct_method", type=click.Choice(["cluster", "adjacency", "directional"]), default="adjacency", show_default=True, help="cluster, adjacency, directional.")
@click.option("--expect_cells", default=3000, show_default=True, type=int, help="Expected recovered cells for RNA cell calling.")
@click.option("--pvalue", default=0.01, show_default=True, type=float, help="P-value cutoff for ambient RNA test.")
@click.option("--force_cells", default=0, show_default=True, type=int, help="Force top N barcodes as cells. 0 means auto cell calling.")
@click.option("--pairendAlignment", "pairend_alignment", is_flag=True, default=False, show_default=True, help="Treat input BAM as paired-end for RNA UMI counting.")
@click.option("--rscript_path", default="Rscript", show_default=True, help="Path to Rscript.")
def step3(**kwargs):
    from .step3 import count_and_call

    count_and_call(**kwargs)


@rna.command(help="Create Seurat RDS from RNA filtered matrix.")
@click.option("--gex_matrix", required=True, type=click.Path(), help="RNA filtered_feature_bc_matrix directory.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", required=True, type=click.Path(), help="Output directory for RDS and Seurat tables.")
@click.option("--core", default=4, show_default=True, type=int, help="Parallel cores for Seurat.")
@click.option("--memory", default=32, show_default=True, type=int, help="Memory limit in GB for Seurat future globals.")
@click.option("--rscript_path", default="Rscript", show_default=True, help="Path to Rscript.")
def rds(**kwargs):
    from ..utils.wrappers import cmd_execute

    rds_script = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "utils",
        "rna_seurat.R"
    )

    if not os.path.exists(rds_script):
        raise FileNotFoundError(f"rna_seurat.R not found: {rds_script}")

    os.makedirs(kwargs["outdir"], exist_ok=True)

    cmd = [
        kwargs["rscript_path"],
        rds_script,
        "--gex_matrix", kwargs["gex_matrix"],
        "--samplename", kwargs["samplename"],
        "--outdir", kwargs["outdir"],
        "--core", str(kwargs["core"]),
        "--memory", str(kwargs["memory"]),
    ]

    cmd_execute(cmd)

@rna.command(help="Create RNA QC summary CSV.")
@click.option("--summary_json", required=True, type=click.Path(), help="RNA summary json from step1/step3.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--outdir", required=True, type=click.Path(), help="Output directory for summary csv.")
@click.option("--tsne", default="", type=click.Path(), help="Optional gex_tsne_umi.xls for mito median.")
def report(**kwargs):
    from .report import make_summary_csv

    make_summary_csv(**kwargs)

@rna.command(help="run all RNA-only steps.")
@click.option("--rnafq1", "rnafq1", required=True, type=click.Path(exists=True, resolve_path=True), multiple=True, help="Expression Read1 fq file, can specify multiple times.")
@click.option("--rnafq2", "rnafq2", required=True, type=click.Path(exists=True, resolve_path=True), multiple=True, help="Expression Read2 fq file, can specify multiple times.")
@click.option("--samplename", required=True, help="Sample name.")
@click.option("--rawname", default="rawname", hidden=True, help="raw name.")
@click.option("--outdir", default="./", show_default=True, type=click.Path(resolve_path=True), help="Output dir.")
@click.option("--shift", is_flag=True, default=False, hidden=True, help="Shift, used to describe read1 structure.")
@click.option("--pattern", "shift_pattern", default="A", hidden=True, help="Anchor sequence, used to describe read1 structure.")
@click.option("--barcode", multiple=True, hidden=True, help="Barcode white list file, can specify multiple times, used to describe read1 structure.")
@click.option("--structure", hidden=True, help="Used to describe read1 structure.")
@click.option("--linker", multiple=True, hidden=True, help="Linker white list file, can specify multiple times, used to describe read1 structure.")
@click.option("--skip_misB", "do_B_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of barcode.")
@click.option("--skip_misL", "do_L_correction", is_flag=True, default=True, show_default=True, help="Not allow one base err correction in each part of linker.")
@click.option("--skip_multi", "use_multi", is_flag=True, default=True, show_default=True, help="Do not rescue barcode match multi when do correction.")
@click.option("--skip_len", "use_short_read", is_flag=True, default=False, show_default=True, help="Skip filtering short reads after adapter filter, short reads will be used.")
@click.option("--chemistry", type=click.Choice(["DD_AG_RNA", "DD_AG", "DD5_AG", "custom"]), default="DD_AG_RNA", show_default=True, help="RNA chemistry.")
@click.option("--core", default=4, show_default=True, help="Set max number of cpus that pipeline might request at the same time.")
@click.option("--memory", default=32, show_default=True, type=int, help="Memory limit in GB for Seurat future globals.")
@click.option("--include-introns", "region", is_flag=True, default=False, callback=include_introns_callback, show_default=True, help="include introns or not.")
@click.option("--refpath", "refpath", required=True, type=click.Path(exists=True, resolve_path=True), help="Path to reference.")
@click.option("--star_path", "star_path", default="STAR", help="Path to executable STAR aligner.")
@click.option("--scoremin", default=None, type=float, help="STAR --outFilterScoreMinOverLread. Example: 0.33 for loose alignment.")
@click.option("--matchnmin", default=None, type=float, help="STAR --outFilterMatchNminOverLread. Example: 0.33 for loose alignment.")
@click.option("--pairendAlignment", "pairend_alignment", is_flag=True, default=False, show_default=True, help="Use paired-end STAR alignment with step1 R1 and R2.")
@click.option("--pairendMinR1Len", "pairend_min_r1_len", default=20, show_default=True, type=int, help="Minimum R1 length kept for paired-end STAR alignment.")
@click.option("--umi_correct_method", type=click.Choice(["cluster", "adjacency", "directional"]), default="adjacency", show_default=True, hidden=True, help="cluster, adjacency, directional.")
@click.option("--expect_cells", default=3000, show_default=True, type=int, help="Expected recovered cells for RNA cell calling.")
@click.option("--pvalue", default=0.01, show_default=True, type=float, help="P-value cutoff for ambient RNA test.")
@click.option("--force_cells", default=0, show_default=True, type=int, help="Force top N barcodes as cells. 0 means auto cell calling.")
@click.option("--rscript_path", default="Rscript", show_default=True, help="Path to Rscript.")
def run(**kwargs):
    logger.info("Check the genomeDir path...")
    kwargs["genomeDir"] = os.path.join(kwargs["refpath"], "star")
    logger.info(check_path(kwargs["genomeDir"]))

    logger.info("Check the gtf path...")
    kwargs["gtf"] = os.path.join(kwargs["refpath"], "genes", "genes.gtf")
    logger.info(check_path(kwargs["gtf"]))

    sampleoutdir = kwargs["outdir"]

    kwargs["gexname"] = f'{kwargs["samplename"]}_E'
    kwargs["outdir"] = os.path.join(sampleoutdir, "analysis", kwargs["gexname"])
    os.makedirs(kwargs["outdir"], exist_ok=True)

    from ..utils.barcode import check_rna_options, barcode_main

    chemistry_kwargs = check_rna_options(fq1=kwargs["rnafq1"], fq2=kwargs["rnafq2"], **kwargs)
    kwargs.update(chemistry_kwargs)

    with open(f"{kwargs['outdir']}/.params.json", "w") as fh:
        json.dump(kwargs, fh, indent=4)

    barcode_main(fq1=kwargs["rnafq1"], fq2=kwargs["rnafq2"], **kwargs)

    step1_fq1 = os.path.join(kwargs["outdir"], "step1", f"{kwargs['gexname']}_1.fq.gz")
    step1_fq2 = os.path.join(kwargs["outdir"], "step1", f"{kwargs['gexname']}_2.fq.gz")

    if kwargs.get("pairend_alignment", False):
        kwargs["fq"] = [step1_fq1, step1_fq2]
        logger.info("RNA STAR alignment mode: paired-end")
    else:
        kwargs["fq"] = [step1_fq2]
        logger.info("RNA STAR alignment mode: single-end R2")

    from ..arc.estep2 import align

    align(**kwargs)

    gexbam = os.path.join(
        kwargs["outdir"],
        "step2",
        "featureCounts",
        f"{kwargs['gexname']}_SortedByName.bam",
    )

    from .step3 import count_and_call

    count_and_call(bam=gexbam, **kwargs)

    # Generate RNA-only Seurat RDS
    filtered_matrix = os.path.join(
        kwargs["outdir"],
        "step3",
        "filtered_feature_bc_matrix",
    )
    step4_outdir = os.path.join(kwargs["outdir"], "step4")
    os.makedirs(step4_outdir, exist_ok=True)

    logger.info("RNA-only RDS generation started.")

    from ..utils.wrappers import cmd_execute

    rds_script = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "utils",
        "rna_seurat.R",
    )

    if not os.path.exists(rds_script):
        raise FileNotFoundError(f"rna_seurat.R not found: {rds_script}")

    cmd_execute([
        kwargs["rscript_path"],
        rds_script,
        "--gex_matrix", filtered_matrix,
        "--samplename", kwargs["samplename"],
        "--outdir", step4_outdir,
        "--core", str(kwargs["core"]),
        "--memory", str(kwargs["memory"]),
    ])

    logger.info("RNA-only RDS generation done.")

    # Generate RNA-only QC summary CSV
    logger.info("RNA-only summary CSV generation started.")

    summary_json = os.path.join(kwargs["outdir"], f"{kwargs['gexname']}_summary.json")
    tsne_file = os.path.join(step4_outdir, "gex_tsne_umi.xls")
    outs_dir = os.path.join(sampleoutdir, "outs")
    os.makedirs(outs_dir, exist_ok=True)

    from .report import make_summary_csv

    make_summary_csv(
        summary_json=summary_json,
        samplename=kwargs["samplename"],
        outdir=outs_dir,
        tsne=tsne_file,
    )

    logger.info("RNA-only summary CSV generation done.")