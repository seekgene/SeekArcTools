import os
import json
import sys
import gzip
from collections import defaultdict
from ..utils.helper import logger
from ..utils.wrappers import cmd_execute

def mapping_summary(STARLog, RnaSeqMetrics):
    summary = defaultdict()
    with open(STARLog, "r") as fh:
        for line in fh:
            if "Number of input reads" in line:
                summary["Number of input reads"] = int(
                    line.strip().split("\t")[-1])
            if "Uniquely mapped reads number" in line:
                summary["Uniquely mapped reads number"] = int(
                    line.strip().split("\t")[-1])
            if "Number of reads mapped to multiple loci" in line:
                summary["Number of reads mapped to multiple loci"] = int(
                    line.strip().split("\t")[-1])
            if "Number of reads mapped to too many loci" in line:
                summary["Number of reads mapped to too many loci"] = int(
                    line.strip().split("\t")[-1])

    with open(RnaSeqMetrics, "r") as fh:
        while True:
            line = fh.readline().strip()
            if line.startswith("total alignments"):
                summary["total alignments"] = int(line.split()[-1].replace(",", ""))
            if line.startswith("reads aligned"):
                summary["reads aligned"] = int(line.split()[-1].replace(",", ""))
            if line.startswith("aligned to genes"):
                summary["aligned to genes"] = int(line.split()[-1].replace(",", ""))
            if line.startswith("no feature assigned"):
                summary["no feature assigned"] = int(line.split()[-1].replace(",", ""))
            if line.startswith("exonic"):
                summary["exonic"] = int(line.split()[-2].replace(",", ""))
            if line.startswith("intronic"):
                summary["intronic"] = int(line.split()[-2].replace(",", ""))
            if line.startswith("intergenic"):
                summary["intergenic"] = int(line.split()[-2].replace(",", ""))
                break
    return summary

def check_bam(bampath):
    if os.path.exists(bampath):
        cmd = f"rm {bampath}"
        cmd_execute(cmd, check=True)

def _read_fastq_record(fh):
    header = fh.readline()
    if not header:
        return None

    seq = fh.readline()
    plus = fh.readline()
    qual = fh.readline()

    if not seq or not plus or not qual:
        raise ValueError("Malformed FASTQ record found during paired-end filtering.")

    return header, seq, plus, qual

def _filter_pairend_fastq(fq1, fq2, outdir, prefix, min_r1_len=20):
    """
    Filter paired-end FASTQ files before STAR paired-end alignment.

    R1 from ARC RNA data may become empty or too short after barcode/UMI trimming.
    STAR paired-end alignment fails when any mate has an empty sequence line.
    This function filters read pairs synchronously:
      - keep pair if len(R1) >= min_r1_len and len(R2) > 0
      - drop both R1 and R2 otherwise
    """
    os.makedirs(outdir, exist_ok=True)

    out_fq1 = os.path.join(outdir, f"{prefix}_pairend_1.fq.gz")
    out_fq2 = os.path.join(outdir, f"{prefix}_pairend_2.fq.gz")

    total = 0
    kept = 0
    drop_short_r1 = 0
    drop_empty_r2 = 0

    with gzip.open(fq1, "rt") as r1_fh, gzip.open(fq2, "rt") as r2_fh, \
            gzip.open(out_fq1, "wt") as out1_fh, gzip.open(out_fq2, "wt") as out2_fh:

        while True:
            r1 = _read_fastq_record(r1_fh)
            r2 = _read_fastq_record(r2_fh)

            if r1 is None and r2 is None:
                break

            if r1 is None or r2 is None:
                raise ValueError("R1 and R2 FASTQ files have different record numbers.")

            total += 1

            r1_len = len(r1[1].rstrip("\r\n"))
            r2_len = len(r2[1].rstrip("\r\n"))

            if r1_len < min_r1_len:
                drop_short_r1 += 1
                continue

            if r2_len == 0:
                drop_empty_r2 += 1
                continue

            out1_fh.writelines(r1)
            out2_fh.writelines(r2)
            kept += 1

    logger.info(
        "paired-end FASTQ filter done. "
        f"total={total}, kept={kept}, "
        f"drop_short_r1={drop_short_r1}, drop_empty_r2={drop_empty_r2}, "
        f"min_r1_len={min_r1_len}"
    )

    if kept == 0:
        logger.error("paired-end FASTQ filter failed: no read pairs left.")
        sys.exit(1)

    return out_fq1, out_fq2

def align(
    fq:list, genomeDir:str, gtf:str, samplename:str, outdir:str, region:str, sc5p=False,
    core:int=4, star_path:str="STAR", **kwargs):

    if ("steps" not in kwargs) or (not kwargs["steps"]):
        kwargs["steps"] = ["STAR", "SortByPos", "qualimap", "FeatureCounts", "SortByName"]

    gexname = f"{samplename}_E"
    basedir = os.path.join(outdir, "step2")
    STAR_dir = os.path.join(basedir, "STAR")
    pairend_filter_dir = os.path.join(basedir, "pairend_filter")

    fq = list(fq)
    pairend_alignment = kwargs.get("pairend_alignment", False)
    if isinstance(pairend_alignment, str):
        pairend_alignment = pairend_alignment.lower() in ("true", "1", "yes", "y")
    pairend_min_r1_len = int(kwargs.get("pairend_min_r1_len", 20) or 20)

    if pairend_alignment:
        if len(fq) != 2:
            logger.error("pairendAlignment requires exactly two fq files.")
            sys.exit(1)

        logger.info(f"STAR paired-end raw input: {fq[0]}, {fq[1]}")

        fq1_filtered, fq2_filtered = _filter_pairend_fastq(
            fq1=fq[0],
            fq2=fq[1],
            outdir=pairend_filter_dir,
            prefix=gexname,
            min_r1_len=pairend_min_r1_len,
        )

        fq = [fq1_filtered, fq2_filtered]
        logger.info(f"STAR paired-end filtered input: {fq[0]}, {fq[1]}")
    else:
        if len(fq) != 1:
            logger.warning(f"single-end alignment expected one fq, got {len(fq)}.")
        logger.info(f"STAR single-end alignment input: {fq}")
    os.makedirs(STAR_dir, exist_ok=True)
    prefix = os.path.join(STAR_dir, gexname + "_")

    if "STAR" not in kwargs["steps"]:
        logger.info("STAR skiped!")
    else:
        from ..utils.wrappers import STAR_wrapper
        logger.info("STAR started!")
        bam, STARLog = STAR_wrapper(
            fq=fq,
            core=core,
            genomeDir=genomeDir,
            prefix=prefix,
            star_path=star_path,
            scoremin=kwargs.get("scoremin", None),
            matchnmin=kwargs.get("matchnmin", None),
        )
        logger.info("STAR done!")
    bam, STARLog =f"{prefix}Aligned.out.bam", f"{prefix}Log.final.out"

    # sort by pos
    if "SortByPos" not in kwargs["steps"]:
        logger.info("SortByPos skiped!")
    else:
        logger.info("SortByPos started!")
        from ..utils.wrappers import samtools_sort_wrapper
        bam = samtools_sort_wrapper(
            bam,
            f"{prefix}SortedByCoordinate.bam",
            core=core,
            clean=True,
        )
        logger.info("SortByPos done!")
    bam = f"{prefix}SortedByCoordinate.bam"

    logger.info("run_qualimap started!")
    from ..utils.wrappers import qualimap_wrapper
    RnaSeqMetrics = qualimap_wrapper(bam=bam, gtf=gtf, outdir=STAR_dir, SC5P=sc5p, pairend_alignment=pairend_alignment,)
    logger.info("run_qualimap done!")

    with open(os.path.join(outdir, f"{gexname}_summary.json")) as fh:
        refpath=os.path.dirname(genomeDir.rstrip("/"))
        reffile=os.path.join(refpath, "reference.json")
        if os.path.exists(reffile):
                with open(reffile) as refjson:
                    refj=json.load(refjson)
                    genome=refj["genomes"][0]
        else:
                genome=genomeDir
        summary = json.load(fh)
        summary["reference"] = genome
        Total = summary["stat"]["total"]
    
    if region=="exon":
        summary["include_introns"] = "False"
    else:
        summary["include_introns"] = "True"

    summary_tmp = defaultdict()
    tmp = mapping_summary(STARLog, RnaSeqMetrics)
    Total = tmp["Number of input reads"]
    mapped_genome_ratio = (tmp["Uniquely mapped reads number"] +
                           tmp["Number of reads mapped to multiple loci"] +
                           tmp["Number of reads mapped to too many loci"])/Total
    summary_tmp["Reads Mapped to Genome"] = mapped_genome_ratio

    mapped_confident_ratio = (tmp["aligned to genes"] + tmp["no feature assigned"]) / Total
    summary_tmp["Reads Mapped Confidently to Genome"] = mapped_confident_ratio

    mapped_intergenic_ratio = tmp["intergenic"] / Total
    summary_tmp["Reads Mapped to Intergenic Regions"] = mapped_intergenic_ratio

    mapped_intronic_ratio = tmp["intronic"] / Total
    summary_tmp["Reads Mapped to Intronic Regions"] = mapped_intronic_ratio

    mapped_exonic_ratio = tmp["exonic"] / Total
    summary_tmp["Reads Mapped to Exonic Regions"] = mapped_exonic_ratio

    with open(os.path.join(outdir, f"{gexname}_summary.json"), "w") as fh:
        summary["mapping"] = summary_tmp
        json.dump(summary, fh, indent=4)


    featureCounts_dir = os.path.join(basedir, "featureCounts")
    # FeatureCounts
    if "FeatureCounts" not in kwargs["steps"]:
        logger.info("run_featureCounts done!")
    else:
        os.makedirs(featureCounts_dir, exist_ok=True)
        logger.info("run_featureCounts started!")
        from ..utils.wrappers import featureCounts_wrapper
        bam = featureCounts_wrapper(
            bam=bam,
            gexname=gexname,
            outdir=featureCounts_dir,
            gtf=gtf,
            SC5P=sc5p,
            region=region,
            core=core,
            pairend_alignment=pairend_alignment,
        )
        logger.info("run_featureCounts done!")
    bam = os.path.join(featureCounts_dir, f"{gexname}_SortedByCoordinate.bam.featureCounts.bam")
    # sort by name
    if "SortByName" not in kwargs["steps"]:
        logger.info("SortByName skiped!")
    else:
        logger.info("SortByName started!")
        from ..utils.wrappers import samtools_sort_wrapper
        bampath = os.path.join(featureCounts_dir, f"{gexname}_SortedByName.bam")
        check_bam(bampath)
        bam = samtools_sort_wrapper(
            bam,
            os.path.join(featureCounts_dir, f"{gexname}_SortedByName.bam"),
            core=core,
            byname=True,
            clean=True,
        )
        logger.info("SortByName done!")
