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
    执行命令，输出命令以进行调试。

    参数：
    args: list or str
        要执行的命令或参数。
    check: bool, 可选
        是否检查命令执行结果，默认为 False。
    text: bool, 可选
        是否以文本形式执行命令，默认为 True。
    capture_output: bool, 可选
        是否捕获命令输出，默认为 True。
    env: os.environ, 可选
        执行命令的环境变量，默认为 None。

    返回值：
    _call
        执行命令的返回结果。

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

    args = ('{bwa_path} mem -t {core} -M -R "@RG\\tID:{atacname}\\tLB:WGS\\tPL:Illumina\\tPU:{atacname}\\tSM:{atacname}" {genomefa} '
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

def qualimap_wrapper(bam:str, gtf:str, outdir:str, qualimap_path:str="qualimap", **kwargs):
    """ -pe,--paired? -s?"""
    strand = {
        'f': 'strand-specific-forward',
        'r': 'strand-specific-reverse',
        'non': 'non-strand-specific'
    }
    s="f"

    if '-pe' in kwargs :
       s = 'non'
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

    for k, v in kwargs.items():
        args += [f"{k}", f"{v}"]

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
    bam:str, gtf:str, gexname:str, outdir:str, region:str, 
    core:int=4, featureCounts_path="featureCounts", **kwargs
    ):
    outcounts = os.path.join(outdir, 'counts.txt')
    region = check_gtf_region(gtf, region)
    s = "1"

    if '-p' in kwargs:
       s = 0
    args = [
        featureCounts_path, "-T", core, "-t", region , "-s", s, "-M", "-O", "-g", "gene_id",
        "--fracOverlap", 0.5, "-a", gtf, "-o", outcounts, 
    ]
    
    for k, v in kwargs.items():
        args += [f"{k}", f"{v}"]
    
    args = args + ["-R", "BAM", bam]
    cmd_execute(args, check=True)
    return os.path.join(outdir, f"{outdir}/{os.path.basename(bam)}.featureCounts.bam")

def bowtie2_wrapper(
    fq:str, ref:str, bam:str, core:int=1, local_mode=True,
    bowtie2_path:str="bowtie2", samtools_path:str="samtools"
)->str:
    """使用bowtie2比对，并用samtools转bam"""
    local_option = ""
    if local_mode:
        local_option = "--local"
    args = (
        f"{bowtie2_path} -p {core} -x {ref} -U {fq} {local_option}|"
        f"samtools view -b > {bam}"
    )
    cmd_execute(args, check=True)
    return bam

def bowtie2_build_wrapper(
    ref:str,
    outdir:str,
    bowtie2_path:str="bowtie2",
)->str:
    os.makedirs(outdir, exist_ok=True)
    if os.path.dirname(os.path.abspath(ref)) != os.path.dirname(os.path.abspath(outdir)):
        ref_l = shutil.copy(ref, outdir)
    args = [f"{bowtie2_path}-build", ref_l, ref_l]
    cmd_execute(args, check=True)
    return ref_l

def igblastn_wrapper(
    input:str,
    output:str,
    organism:str,
    chain:str,
    core:int=4,
    igblastn_path:str="igblastn",
    auxiliary_data:str="",
    internal_data:str="",
):
    outdir = os.path.dirname(output)

    cmd = f"cd {outdir}; ln -s {internal_data} ./{organism}; "
    if chain=="TR":
        cmd += (
            f"{igblastn_path} -query {input} -auxiliary_data {auxiliary_data} -outfmt 19 -ig_seqtype TCR "
            f"-germline_db_J {organism}/TR_J-REGION -germline_db_D {organism}/TR_D-REGION "
            f"-germline_db_V '{organism}/TR_L-REGION+V-REGION' -c_region_db {organism}/TR_C-REGION "
            f"-num_threads {core} -organism {organism} -out {output}"
        )
    elif chain=="IG":
        cmd += (
            f"{igblastn_path} -query {input} -auxiliary_data {auxiliary_data} -outfmt 19 -ig_seqtype Ig "
            f"-germline_db_J {organism}/IG_J-REGION -germline_db_D {organism}/IG_D-REGION "
            f"-germline_db_V {organism}/IG_L-REGION+V-REGION -c_region_db {organism}/IG_C-REGION  "
            f"-num_threads {core} -organism {organism} -out {output}"
        )
    else:
        pass

    cmd_execute(cmd, check=True)
    return output

def picard_DownsampleSam_wrapper(bam, fra, downsambam):
    cmd = "picard DownsampleSam I=%s O=%s P=%s" %(bam, downsambam, fra)
    cmd_execute(cmd)
    samtools_index_wrapper(downsambam)
    return downsambam

def geneBody_coverage_wrapper(downsambam, downdir, gexname, gtf):
    resultbed = os.path.join(downdir, gexname+".bed")
    reductionbed = os.path.join(downdir, gexname+".reduction.bed")
    pel=os.path.join(os.path.abspath(os.path.dirname(__file__)),'gtf2bed.pl')
    with open(resultbed, 'wb', 0) as bedfile:
        subprocess.call(['perl', pel, gtf], stdout=bedfile)
    with open(resultbed) as infile:
        with open(reductionbed,'w') as out:
            allbed=infile.readlines()
            lines = random.sample(allbed,20000)
            for line in lines:
                out.write(line)
    cmd = f"cd {downdir}; geneBody_coverage.py -r {downdir}/{gexname}.reduction.bed -i {downsambam} -o {downdir}/{gexname}"
    cmd_execute(cmd)
    with open(f'{downdir}/{gexname}.geneBodyCoverage.txt','r') as infile:
        infile.readline()
        txt=infile.readline().strip().split('\t')[1:]
        numberlist=map(float,txt)
        all=sum(numberlist)
        mosttxt=txt[24:76]
        mostnumberlist=map(float,mosttxt)
        most=sum(mostnumberlist)
    return float(most)/float(all)
def cr_vdj_wrapper(
    projdir:str,
    fastqs:str,
    sample:str,
    cr_path:str,
    reference:str,
    chain:str="auto",
    core:int=4,
    mem:int=30,
):
    cmd = (f"mkdir -p {projdir}; cd {projdir}; "
        f"{cr_path} vdj --id={sample} --fastqs={fastqs} --sample={sample} --chain {chain} "
        f" --reference={reference} --localcores={core} --localmem={mem}")
    cmd_execute(cmd)
    return projdir
