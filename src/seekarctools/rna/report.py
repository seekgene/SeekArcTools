import json
import os

import pandas as pd

from ..utils.helper import logger


HEADER = (
    "Samplename,"
    "Pipeline_version,"
    "Chemistry,"
    "Include_introns,"
    "Estimated_number_of_cells,"
    "GEX_Sequenced_read_pairs,"
    "GEX_Valid_barcodes,"
    "GEX_Too_Short,"
    "GEX_Q30_bases_in_barcode,"
    "GEX_Q30_bases_in_UMI,"
    "GEX_Fraction_of_reads_in_cells,"
    "GEX_Mean_raw_reads_per_cell,"
    "GEX_Median_UMI_counts_per_cell,"
    "GEX_Median_genes_per_cell,"
    "GEX_Total_Genes_Detected,"
    "GEX_Sequencing_Saturation,"
    "GEX_Mito_Median,"
    "GEX_Reads_mapped_to_genome,"
    "GEX_Reads_mapped_confidently_to_genome,"
    "GEX_Reads_mapped_confidently_to_intergenic_regions,"
    "GEX_Reads_mapped_confidently_to_intronic_regions,"
    "GEX_Reads_mapped_confidently_to_exonic_regions,"
    "Contamination_rate"
)


def _safe_get(d, *keys, default=""):
    cur = d
    for key in keys:
        if not isinstance(cur, dict):
            return default
        if key not in cur:
            return default
        cur = cur[key]
    return cur


def _fmt_float(value, ndigits=4, default=""):
    try:
        if value in ["", None]:
            return default
        return float(f"{float(value):.{ndigits}f}")
    except Exception:
        return default


def _fmt_int(value, default=""):
    try:
        if value in ["", None]:
            return default
        return int(float(value))
    except Exception:
        return default


def _safe_divide(a, b, ndigits=4, default=""):
    try:
        if a in ["", None] or b in ["", None]:
            return default
        if float(b) == 0:
            return default
        return float(f"{float(a) / float(b):.{ndigits}f}")
    except Exception:
        return default


def _q30_rate(qdict):
    """
    Calculate Q30 base ratio from summary["barcode_q"] or summary["umi_q"].
    """
    if not isinstance(qdict, dict):
        return ""

    total_base = 0
    q30_base = 0

    try:
        for values in qdict.values():
            if not isinstance(values, list):
                continue
            total_base += sum(values)
            q30_base += sum(values[30:])
    except Exception:
        return ""

    if total_base == 0:
        return ""

    return float(f"{q30_base / total_base:.4f}")


def _get_mito_median(tsne_file):
    """
    gex_tsne_umi.xls from rna_seurat.R contains mito / percent.mito information.
    Return median mito as a percent-like string, same display style as ARC report.
    """
    if not tsne_file:
        return "0.00%"

    if not os.path.exists(tsne_file):
        logger.info(f"Warning: The path of '{tsne_file}' is not exists.")
        return "0.00%"

    try:
        df = pd.read_table(tsne_file, index_col=0)
    except Exception as e:
        logger.info(f"Warning: failed to read '{tsne_file}': {e}")
        return "0.00%"

    mito_col = None
    for col in ["mito", "percent.mito", "percent_mito"]:
        if col in df.columns:
            mito_col = col
            break

    if mito_col is None:
        return "0.00%"

    try:
        median_mito = df[mito_col].median()
        if pd.isna(median_mito):
            return "0.00%"
        return f"{float(median_mito):.2f}%"
    except Exception:
        return "0.00%"


def _contamination_rate(summary):
    """
    Keep the same idea as parse_json_arc():
    if seq_17L19ME exists, contamination = seq_17L19ME / total.
    """
    total = _safe_get(summary, "stat", "total", default="")
    contam = _safe_get(summary, "stat", "seq_17L19ME", default="")

    try:
        if total in ["", None] or float(total) == 0:
            return ""
        if contam in ["", None]:
            return ""
        return f"{float(contam) / float(total):.2%}"
    except Exception:
        return ""


def make_summary_csv(summary_json, outdir, samplename, tsne="", rawname=None, **kwargs):
    """
    RNA-only QC summary csv.

    Input:
        *_E_summary.json
        optional step4/gex_tsne_umi.xls

    Output:
        {samplename}_summary.csv

    This file intentionally contains only RNA/GEX metrics.
    No ATAC / joint / peak / fragment / linkage fields are generated here.
    """
    os.makedirs(outdir, exist_ok=True)

    if not os.path.exists(summary_json):
        raise FileNotFoundError(f"summary_json not found: {summary_json}")

    with open(summary_json) as fh:
        gex_summary = json.load(fh)

    pipeline_version = gex_summary.get("__version__", "")
    if pipeline_version:
        pipeline_version = "seekarctools " + str(pipeline_version)

    total_reads = _safe_get(gex_summary, "stat", "total", default="")
    valid_reads = _safe_get(gex_summary, "stat", "valid", default="")

    summary_data = [
        samplename,
        pipeline_version,
        _safe_get(gex_summary, "stat", "chemistry", default=""),
        _safe_get(gex_summary, "include_introns", default=""),
        _fmt_int(_safe_get(gex_summary, "cells", "Estimated Number of Cells")),

        _fmt_int(total_reads),
        _safe_divide(valid_reads, total_reads),
        _fmt_int(_safe_get(gex_summary, "stat", "too_short")),
        _q30_rate(gex_summary.get("barcode_q", {})),
        _q30_rate(gex_summary.get("umi_q", {})),

        _fmt_float(_safe_get(gex_summary, "cells", "Fraction Reads in Cells")),
        _fmt_float(_safe_get(gex_summary, "cells", "Mean Reads per Cell")),
        _fmt_int(_safe_get(gex_summary, "cells", "Median UMI Counts per Cell")),
        _fmt_int(_safe_get(gex_summary, "cells", "Median Genes per Cell")),
        _fmt_int(_safe_get(gex_summary, "cells", "Total Genes Detected")),
        _fmt_float(_safe_get(gex_summary, "cells", "Sequencing Saturation")),
        _get_mito_median(tsne),

        _fmt_float(_safe_get(gex_summary, "mapping", "Reads Mapped to Genome")),
        _fmt_float(_safe_get(gex_summary, "mapping", "Reads Mapped Confidently to Genome")),
        _fmt_float(_safe_get(gex_summary, "mapping", "Reads Mapped to Intergenic Regions")),
        _fmt_float(_safe_get(gex_summary, "mapping", "Reads Mapped to Intronic Regions")),
        _fmt_float(_safe_get(gex_summary, "mapping", "Reads Mapped to Exonic Regions")),

        _contamination_rate(gex_summary),
    ]

    out_csv = os.path.join(outdir, f"{samplename}_summary.csv")

    with open(out_csv, "w") as fh:
        fh.write(HEADER + "\n")
        fh.write(",".join(str(x).replace(",", "") for x in summary_data) + "\n")

    logger.info(f"RNA-only summary csv saved: {out_csv}")
    return out_csv