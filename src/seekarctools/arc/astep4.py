import os
from ..utils.helper import logger, get_file_path
from ..utils.wrappers import cmd_execute

def do_signac(
    gex_matrix:str, atac_matrix:str, fragpath:str, species:str, anno_rds:str, outdir:str, rscriptpath:str="Rscript", core:int=4, memory:int=60, **kwargs):

    logger.info("signac started!")
    outdir = os.path.join(outdir, "step4")
    os.makedirs(outdir, exist_ok=True)

    Rapp = os.path.abspath(os.path.join(get_file_path(__file__),"../utils/count_link.R"))
    args = [
        rscriptpath, Rapp, "--gex_matrix", gex_matrix, "--atac_matrix", atac_matrix, "--fragpath", fragpath, "--outdir", outdir,
        "--core", core, "--species", species, "--anno_rds", anno_rds, "--memory", memory
    ]
    cmd_execute(args, check=True)

    logger.info("signac done!")
