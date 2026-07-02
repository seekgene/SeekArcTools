import os
import json
import gzip
import sys
import numpy as np
from collections import defaultdict
from ..utils.helper import logger, get_utils_path
from ..utils.wrappers import cmd_execute
from ..utils import countUtil 

__srcdir = get_utils_path(__file__)

def count(bam, outdir, gtf, umi_correct_method, **kwargs):
    basedir = os.path.join(outdir, "step3")
    os.makedirs(basedir, exist_ok=True)
    countUtil.count(bam, basedir, gtf, umi_correct_method, **kwargs)

def _count_barcodes(barcodes_file):
    n = 0
    with gzip.open(barcodes_file, "rt") as fh:
        for line in fh:
            if line.strip():
                n += 1
    return n


def cell_calling(
    outdir,
    samplename,
    expect_cells=3000,
    pvalue=0.01,
    force_cells=0,
    rscript_path="Rscript",
    **kwargs
):
    basedir = os.path.join(outdir, "step3")
    raw_matrix_dir = os.path.join(basedir, "raw_feature_bc_matrix")
    filtered_matrix_dir = os.path.join(basedir, "filtered_feature_bc_matrix")

    if not os.path.exists(raw_matrix_dir):
        logger.error(f"raw matrix not found: {raw_matrix_dir}")
        sys.exit(1)

    cell_identify_R = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "utils",
        "cell_identify.R"
    )

    cmd = [
        rscript_path,
        cell_identify_R,
        "-i", raw_matrix_dir,
        "-o", filtered_matrix_dir,
        "-e", str(expect_cells),
        "-p", str(pvalue),
    ]

    if force_cells and int(force_cells) > 0:
        cmd += ["-f", str(int(force_cells))]

    logger.info("ARC GEX cell calling started.")
    ret = cmd_execute(cmd, check=False)

    if ret.returncode != 0:
        logger.error("cell_identify.R returned non-zero exit code.")
        if ret.stdout:
            logger.error(f"cell_identify.R stdout: {ret.stdout}")
        if ret.stderr:
            logger.error(f"cell_identify.R stderr: {ret.stderr}")
        sys.exit(ret.returncode)

    filtered_barcodes_file = os.path.join(filtered_matrix_dir, "barcodes.tsv.gz")

    if not os.path.exists(filtered_barcodes_file):
        logger.error(f"filtered barcodes not found: {filtered_barcodes_file}")
        sys.exit(1)

    cell_num = _count_barcodes(filtered_barcodes_file)
    if cell_num == 0:
        logger.error("cell calling failed: filtered barcodes number is 0.")
        sys.exit(1)

    logger.info(f"ARC GEX cell calling done. Estimated cells: {cell_num}")
