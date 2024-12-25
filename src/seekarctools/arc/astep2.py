import os
import json
from collections import defaultdict
from ..utils.helper import logger
from ..utils.wrappers import bwamem_wrapper


def align(
    afq:list, genomefa:str, atacname:str, outdir:str, core:int=4, bwa_path:str="bwa", **kwargs):
    basedir = os.path.join(outdir, "step2")
    bwa_dir = os.path.join(basedir, "bwa_pe")
    os.makedirs(bwa_dir, exist_ok=True)
    prefix = os.path.join(bwa_dir, atacname + "_")

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
