import os
import psutil
from ..utils.helper import logger, get_file_path
from ..utils.wrappers import cmd_execute

def do_signac(
    gex_matrix:str, atac_matrix:str, fragpath:str, samplename:str, refpath:str, outdir:str, 
    rscriptpath:str="Rscript", rawname:str="rawname", core:int=4, **kwargs):

    logger.info("signac started!")
    outdir = os.path.join(outdir, "step4")
    ref_gtf = os.path.join(refpath, "genes", "genes.gtf")
    ref_fa = os.path.join(refpath, "fasta", "genome.fa")
    os.makedirs(outdir, exist_ok=True)
    if rawname == "rawname":
        rawname = samplename

    available_memory = psutil.virtual_memory().available / (1024 * 1024 * 1024)
    memory = int(available_memory * 0.9)

    Rapp = os.path.abspath(os.path.join(get_file_path(__file__),"../utils/count_link.R"))
    args = [
        rscriptpath, Rapp, "--gex_matrix", gex_matrix, "--atac_matrix", atac_matrix, "--fragpath", fragpath, "--rawname", rawname, "--samplename", samplename, "--outdir", outdir,
        "--core", core, "--memory", memory, "--ref_gtf", ref_gtf, "--ref_fa", ref_fa
    ]
    cmd_execute(args, check=True)

    logger.info("signac done!")
