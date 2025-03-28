import os
import psutil
from ..utils.helper import logger, get_file_path
from ..utils.wrappers import cmd_execute

def do_signac(
    gex_matrix:str, atac_matrix:str, fragpath:str, samplename:str, organism:str, outdir:str, 
    rscriptpath:str="Rscript", rawname:str="rawname", anno_rds:str="anno.rds", core:int=4, **kwargs):

    logger.info("signac started!")
    outdir = os.path.join(outdir, "step4")
    os.makedirs(outdir, exist_ok=True)
    if rawname == "rawname":
        rawname = samplename

    available_memory = psutil.virtual_memory().available / (1024 * 1024 * 1024)
    memory = int(available_memory * 0.9)

    Rapp = os.path.abspath(os.path.join(get_file_path(__file__),"../utils/count_link.R"))
    args = [
        rscriptpath, Rapp, "--gex_matrix", gex_matrix, "--atac_matrix", atac_matrix, "--fragpath", fragpath, "--rawname", rawname, "--samplename", samplename, "--outdir", outdir,
        "--core", core, "--memory", memory, "--organism", organism, "--anno_rds", anno_rds
    ]
    cmd_execute(args, check=True)

    logger.info("signac done!")
