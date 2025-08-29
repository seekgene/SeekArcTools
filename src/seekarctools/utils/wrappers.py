import sys
import os
import shutil
from subprocess import run
import subprocess
import random
import psutil
from .helper import logger

def cmd_execute(
    args, check:bool=False, text:bool=True,
    capture_output:bool=True, env:os.environ=None
):
    """
    Parameters:
    args: list or str
        The command or arguments to execute.
    check: bool, optional
        Whether to check the command execution result, default is False.
    text: bool, optional
        Whether to execute the command in text mode, default is True.
    capture_output: bool, optional
        Whether to capture the command output, default is True.
    env: os.environ, optional
        The environment variables for executing the command, default is None.

    Returns:
    _call
        The result returned by executing the command.
    """

    if isinstance(args, list):
        args = [str(_) for _ in args]
        logger.info(" ".join(args))
        _call = run(args, check=False, text=text, capture_output=capture_output, env=env)
    elif isinstance(args, str):
        logger.info(args)
        _call = run(args, shell=True, check=False, text=text, capture_output=capture_output, env=env)
    
    if check:
        if _call.returncode!=0:
            logger.info(f"stderr: {_call.stderr}")
            logger.info(f"stdout: {_call.stdout}")
            sys.exit(1)
    return _call

def STAR_wrapper(
    fq:list, genomeDir:str, prefix:str, scoremin=None, matchnmin=None,
    core:int=4, star_path:str="STAR", readMapNumber:int=-1,
    ):
    """wrapper for STAR"""
    args = [
        star_path, '--runThreadN', core, '--limitOutSJcollapsed', 5000000, '--readMapNumber', readMapNumber,
        '--genomeDir', genomeDir, '--readFilesCommand', 'zcat', '--outFileNamePrefix',
        prefix, '--outSAMtype', 'BAM', 'Unsorted', 
    ]
    if scoremin:
        args += ['--outFilterScoreMinOverLread', scoremin,]
    if matchnmin:
        args += [ '--outFilterMatchNminOverLread', matchnmin,]

    args = args + ['--readFilesIn', ] + list(fq)
    _ = cmd_execute(args)
    return f'{prefix}Aligned.out.bam', f'{prefix}Log.final.out'

def bwamem_wrapper(
    afq:list, genomefa:str, atacname:str, prefix:str, core:int=4, samtools_path:str="samtools", bwa_path:str="bwa", 
    ):
    """wrapper for bwa-mem"""

    args = ('{bwa_path} mem -t {core} -M -I 250[150,1000,36] -R "@RG\\tID:{atacname}\\tLB:WGS\\tPL:Illumina\\tPU:{atacname}\\tSM:{atacname}" {genomefa} '
            '{astep1_fq1} {astep1_fq2} | {samtools_path} sort -@ {core} -o {prefix}mem_pe_Sort.bam').format(
                bwa_path=bwa_path, core=core, atacname=atacname, genomefa=genomefa, astep1_fq1=afq[0], astep1_fq2=afq[-1], samtools_path=samtools_path, prefix=prefix)
    _ = cmd_execute(args)
    return f'{prefix}mem_pe_Sort.bam'


def samtools_index_wrapper(inbam:str, core:int=4, samtools_path:str="samtools"):
    args = [samtools_path, "index", "-@", core, inbam]
    _call = cmd_execute(args, check=True)
    return 

def samtools_sort_wrapper(
    inbam:str, outbam:str, byname:bool=False, clean:bool=False,
    core:int=4, samtools_path:str="samtools"
)->str:
    """wrapper for samtools sort"""
    args = [
        samtools_path, "sort", "-O", "BAM", "-@",  core, "-o", outbam, inbam
    ]
    if byname:
        args.insert(2, '-n')

    _call = cmd_execute(args, check=True)

    if _call.returncode == 0:
        # index
        if not byname:
            samtools_index_wrapper(outbam,samtools_path=samtools_path)
        if clean:
            os.remove(inbam)
    return outbam

def unzip_wrapper(gzfile:str, outfile:str):
    args = [
        "gzip", "-dc", gzfile, ">", outfile
    ]
    cmd_execute(args, check=True)
    return outfile

def qualimap_wrapper(bam:str, gtf:str, outdir:str, SC5P:bool=None, qualimap_path:str="qualimap", **kwargs):
    strand = {
        'f': 'strand-specific-forward',
        'r': 'strand-specific-reverse',
        'non': 'non-strand-specific'
    }
    s = "non"
    if SC5P:
        s="r"
    # gtf.gz?
    if gtf.endswith('.gz'):
        gtf_file = os.path.join(outdir, 'tmp.gtf')
        unzip_wrapper(gtf, gtf_file)
    else:
        gtf_file = gtf
    available_memory = psutil.virtual_memory().available / (1024 * 1024 * 1024)
    args = [
        qualimap_path, 'rnaseq', '-outformat', 'PDF', '-outdir', outdir, '-bam',
        bam, '-gtf', gtf, '-p', strand[s], f'--java-mem-size={int(available_memory)}G'
    ]

    my_env = os.environ.copy()
    if 'DISPLAY' in my_env:
        del my_env['DISPLAY']
    
    cmd_execute(args, check=True, env=my_env)
    return os.path.join(outdir, "rnaseq_qc_results.txt")

def check_gtf_region(gtf, region):
    region_set = set()
    logger.info("check gtf feature.")
    with open(gtf) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            tmp = line.strip().split("\t")
            region_set.add(tmp[2])
            if tmp[2] == region:
                return region
    if region == "transcript" and "gene" in region_set:
        logger.warning("No transcript feature in gtf, use gene feature instead.")
        return "gene"
    logger.warning(f"No {region} feature in gtf")
    os.exit(1)

def featureCounts_wrapper(
    bam:str, gtf:str, gexname:str, outdir:str, region:str, SC5P:bool=None, 
    core:int=4, featureCounts_path="featureCounts", **kwargs
    ):
    outcounts = os.path.join(outdir, 'counts.txt')
    region = check_gtf_region(gtf, region)
    s = 0
    if SC5P:
        s = 2

    args = [
        featureCounts_path, "-T", core, "-t", region , "-s", s, "-M", "-O", "-g", "gene_id",
        "--fracOverlap", 0.5, "-a", gtf, "-o", outcounts, 
    ]
    
    args = args + ["-R", "BAM", bam]
    cmd_execute(args, check=True)
    return os.path.join(outdir, f"{outdir}/{os.path.basename(bam)}.featureCounts.bam")
