import os
import re
import json
from collections import defaultdict
from subprocess import run
from ..utils.helper import logger
from ..utils.wrappers import bwamem_wrapper, cmd_execute

def check_bam(bampath):
    if os.path.exists(bampath):
        cmd = f"rm {bampath}"
        cmd_execute(cmd, check=True)

def parse_bam_stats(filename):
    total_reads = None
    map_pairs_reads = None
    unmapped_reads = None

    if not os.path.isfile(filename):
        raise FileNotFoundError(f"{filename} file not found!")

    with open(filename, "rt") as fh:
        for line in fh:
            if not line.startswith("SN"):
                continue
            match = re.search(r"raw total sequences:\s+(\d+)", line)
            if match:
                total_reads = int(match.group(1))
            match = re.search(r"reads mapped and paired:\s+(\d+)", line)
            if match:
                map_pairs_reads = int(match.group(1))
            match = re.search(r"reads unmapped:\s+(\d+)", line)
            if match:
                unmapped_reads = int(match.group(1))

    missing = []
    if total_reads is None:
        missing.append("raw total sequences")
    if map_pairs_reads is None:
        missing.append("reads mapped and paired")
    if unmapped_reads is None:
        missing.append("reads unmapped")

    if missing:
        missing_str = ", ".join(missing)
        logger.info(f"Warning: {missing_str} was None. Setting Confidently mapped read pairs and Unmapped read pairs to 0.")

    return total_reads, map_pairs_reads, unmapped_reads

def align(
    afq:list, genomefa:str, samplename:str, outdir:str, core:int=4, bwa_path:str="bwa", samtoolspath:str="samtools", **kwargs):
    atacname = f"{samplename}_A"
    basedir = os.path.join(outdir, "step2")
    bwa_dir = os.path.join(basedir, "bwa_pe")
    os.makedirs(bwa_dir, exist_ok=True)
    bamname = f'{atacname}_mem_pe_Sort.bam'
    bampath = os.path.join(bwa_dir, bamname)
    prefix = os.path.join(bwa_dir, atacname + "_")
    check_bam(bampath)

    logger.info("BWA mem started!")
    bam = bwamem_wrapper(
        afq=afq,
        genomefa=genomefa,
        atacname=atacname,
        prefix=prefix,
        core=core,
        bwa_path=bwa_path
    )
    logger.info("BWA mem done!")

    logger.info("count bam start...")
    cmd = ("cd {bwa_dir}; "
            "{samtoolspath} stats -@ {core} {atacname}_mem_pe_Sort.bam > {atacname}_mem_pe_Sort.bam.stats; "
            "{samtoolspath} view -@ {core} -b -q 30 {atacname}_mem_pe_Sort.bam > {atacname}_mem_pe_Sort_q30.bam; "
            "{samtoolspath} stats -@ {core} {atacname}_mem_pe_Sort_q30.bam > {atacname}_mem_pe_Sort_q30.bam.stats && rm -f {atacname}_mem_pe_Sort_q30.bam"
            ).format(bwa_dir=bwa_dir, samtoolspath=samtoolspath, atacname=atacname, core=core)
    logger.info(cmd)
    run(cmd, shell=True)

    bam_stats = os.path.join(bwa_dir, f'{atacname}_mem_pe_Sort.bam.stats')
    q30bam_stats = os.path.join(bwa_dir, f'{atacname}_mem_pe_Sort_q30.bam.stats')
    total_reads, map_pairs_reads, unmapped_reads = parse_bam_stats(bam_stats)

    q30total_reads, q30map_pairs_reads, q30unmapped_reads = parse_bam_stats(q30bam_stats)

    confidently_mapped_read_pairs = q30map_pairs_reads / total_reads
    unmapped_read_pairs = 1 - (map_pairs_reads / total_reads)

    mapd = {}
    mapd["confidently_mapped_read_pairs"] = confidently_mapped_read_pairs
    mapd["unmapped_read_pairs"] = unmapped_read_pairs

    atacjson = os.path.join(outdir, atacname+"_summary.json")
    with open(atacjson, "r") as fh:
        atac_summary = json.load(fh)

    with open(atacjson, "w") as fh:
        atac_summary["bam_count"] = mapd
        json.dump(
            atac_summary,
            fh,
            indent=4,
            default=lambda o: int(o) if isinstance(o, np.int64) else o
        )
    logger.info("count bam end!")