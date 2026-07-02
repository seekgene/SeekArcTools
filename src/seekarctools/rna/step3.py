import gzip
import json
import os
import sys

import numpy as np

from ..utils import countUtil
from ..utils.helper import logger
from ..utils.wrappers import cmd_execute


def _count_barcodes(barcodes_file):
    n = 0
    with gzip.open(barcodes_file, "rt") as fh:
        for line in fh:
            if line.strip():
                n += 1
    return n


def _json_default(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def count_and_call(
    bam,
    outdir,
    samplename,
    gtf,
    umi_correct_method="adjacency",
    expect_cells=3000,
    pvalue=0.01,
    force_cells=0,
    rscript_path="Rscript",
    **kwargs
):
    """
    RNA-only step3:
    1. count UMI and write raw_feature_bc_matrix
    2. call cells from raw_feature_bc_matrix by cell_identify.R
    3. write filtered_feature_bc_matrix
    4. calculate RNA cell metrics and update summary json
    """

    gexname = kwargs.get("gexname", f"{samplename}_E")
    basedir = os.path.join(outdir, "step3")
    os.makedirs(basedir, exist_ok=True)

    pairend_alignment = kwargs.get("pairend_alignment", False)

    logger.info("RNA count started.")
    logger.info(f"RNA UMI counting mode: {'paired-end' if pairend_alignment else 'single-end'}")
    countUtil.count(
        bam=bam,
        outdir=basedir,
        gtf=gtf,
        umi_correct_method=umi_correct_method,
        ispair=pairend_alignment,
        **kwargs
    )
    logger.info("RNA count done.")

    raw_matrix_dir = os.path.join(basedir, "raw_feature_bc_matrix")
    filtered_matrix_dir = os.path.join(basedir, "filtered_feature_bc_matrix")

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
        "-e", expect_cells,
        "-p", pvalue,
    ]

    if force_cells and int(force_cells) > 0:
        cmd += ["-f", int(force_cells)]

    logger.info("RNA cell calling started.")
    ret = cmd_execute(cmd, check=False)

    if ret.returncode != 0:
        logger.warning("cell_identify.R returned non-zero exit code.")
        if ret.stdout:
            logger.warning(f"cell_identify.R stdout: {ret.stdout}")
        if ret.stderr:
            logger.warning(f"cell_identify.R stderr: {ret.stderr}")

    filtered_barcodes_file = os.path.join(filtered_matrix_dir, "barcodes.tsv.gz")
    filtered_features_file = os.path.join(filtered_matrix_dir, "features.tsv.gz")

    if not os.path.exists(filtered_barcodes_file):
        logger.error(f"RNA cell calling failed: {filtered_barcodes_file} not found.")
        sys.exit(1)

    cell_num = _count_barcodes(filtered_barcodes_file)
    if cell_num == 0:
        logger.error("RNA cell calling failed: filtered barcodes number is 0.")
        sys.exit(1)

    logger.info(f"RNA cell calling done. Estimated cells: {cell_num}")

    logger.info("RNA metrics calculation started.")
    summary_tmp, downsample = countUtil.calculate_metrics(
        counts_file=os.path.join(basedir, "counts.xls"),
        detail_file=os.path.join(basedir, "detail.xls"),
        filterd_barcodes_file=filtered_barcodes_file,
        filterd_features_file=filtered_features_file,
        gtf=gtf,
        basedir=basedir
    )

    summary_json = os.path.join(outdir, f"{gexname}_summary.json")

    if os.path.exists(summary_json):
        with open(summary_json, "r") as fh:
            summary = json.load(fh)
    else:
        logger.warning(f"{summary_json} not found. A new summary json will be created.")
        summary = {}

    total_reads = summary.get("stat", {}).get("total", 0)
    estimated_cell_num = summary_tmp["Estimated Number of Cells"]

    if estimated_cell_num > 0 and total_reads > 0:
        mean_reads_per_cell = int(total_reads / estimated_cell_num)
        summary_tmp["Mean Reads per Cell"] = mean_reads_per_cell

    if "Median lnc Genes per Cell" in summary_tmp:
        del summary_tmp["Median lnc Genes per Cell"]

    summary["cells"] = summary_tmp

    if "Mean Reads per Cell" in summary_tmp:
        downsample["Reads"] = [
            int(summary_tmp["Mean Reads per Cell"] * p)
            for p in downsample["percentage"]
        ]

    summary["downsample"] = downsample

    with open(summary_json, "w") as fh:
        json.dump(summary, fh, indent=4, default=_json_default)

    logger.info("RNA metrics calculation done.")