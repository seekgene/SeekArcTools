import json
import os
from jinja2 import Environment, FileSystemLoader
import numpy as np
import pandas as pd
from ..utils.helper import logger

def pre_diff_data(diff_table, n=50):
    df = pd.read_table(diff_table)
    tmp = df.groupby('cluster').apply(lambda x: x.sort_values(by='avg_log2FC', ascending=False).head(n))
    diff_data = tmp.loc[:,['Ensembl', 'gene', 'avg_log2FC', 'p_val_adj','cluster']].to_dict('records')
    return diff_data

def get_gex_tsne(tsnefile):
    df = pd.read_table(tsnefile, index_col=0)
    # df = df.drop(['orig.ident', 'nFeature_RNA', 'percent.mito'], axis=1)
    tsne_d = df.to_dict(orient='list')
    tsne1 = tsne_d['tSNE_1']
    tsne2 = tsne_d['tSNE_2']
    depth = tsne_d['nCount_RNA']
    cluster = tsne_d['seurat_clusters']
    return tsne1, tsne2, depth, cluster

def get_atac_tsne(tsnefile):
    df = pd.read_table(tsnefile, index_col=0)
    # df = df.drop(['frac_dup', 'frac_mito', 'tsse', 'n_frag_overlap_peak', 'percent.mito', ], axis=1)
    tsne_d = df.to_dict(orient='list')
    tsne1 = tsne_d['atactsne_1']
    tsne2 = tsne_d['atactsne_2']
    depth = tsne_d['nCount_ATAC']
    cluster = tsne_d['seurat_clusters']
    return tsne1, tsne2, depth, cluster

def count_link(linkfile):
    linkdf = pd.read_table(linkfile)
    total_rows = linkdf.shape[0]
    unique_genes = linkdf['gene'].nunique()
    unique_peaks = linkdf['peak'].nunique()
    return total_rows, unique_genes, unique_peaks


def report(gexjson, atacjson, outdir, samplename, rawname, **kwargs):
    os.makedirs(outdir, exist_ok=True)
    datajson = os.path.join(os.path.dirname(__file__), '../utils/reportarc/sgarc.json')
    assert os.path.exists(datajson), f'{datajson} not found!'

    with open(gexjson) as gfh:
        gex_summary = json.load(gfh)
    with open(atacjson) as afh:
        atac_summary = json.load(afh)
    with open(datajson) as dfh:
        data_summary = json.load(dfh)

    atacdir = os.path.dirname(atacjson)
    linkfile = os.path.join(atacdir, 'step4', 'linked_feature.xls')
    if not os.path.exists(linkfile):
        logger.info(f"Warning : The path of '{linkfile}' is not exists. This may be due to the reference genome may be in NCBI format")
        feature_link = 0
        link_gene = 0
        link_peak = 0
    else:
        feature_link, link_gene, link_peak = count_link(linkfile)
    if rawname == "rawname":
        rawname = samplename

    # joint: title
    data_summary["Joint"][0]["left"][0]["data"]["Estimated number of cells"] = f'{atac_summary["atac"]["Cells"]["Estimated number of cells"]:,}'
    data_summary["Joint"][0]["right"][0]["data"]["GEX Median genes per cell"] = f'{gex_summary["cells"]["Median Genes per Cell"]:,}'
    data_summary["Joint"][0]["right"][0]["data"]["ATAC Median high-quality fragments overlapping peaks per cell"] = f'{atac_summary["atac"]["Cells"]["Median high-quality fragments in peaks per cell"]:,}'
    # joint: left_Sample riget_Joint Metrics
    data_summary["Joint"][1]["left"][0]["data"]["Sample ID"] = samplename
    data_summary["Joint"][1]["left"][0]["data"]["Sample description"] = rawname
    data_summary["Joint"][1]["left"][0]["data"]["Pipeline version"] = "seekarctools " + gex_summary["__version__"]
    data_summary["Joint"][1]["left"][0]["data"]["Reference path"] = atac_summary["atac"]["refpath"].split('/')[-1]
    data_summary["Joint"][1]["left"][0]["data"]["Chemistry"] = atac_summary["stat"]["chemistry"]
    data_summary["Joint"][1]["left"][0]["data"]["Include introns"] = gex_summary["include_introns"]
    data_summary["Joint"][1]["right"][0]["data"]["Feature linkages detected"] = f'{feature_link:,}'
    data_summary["Joint"][1]["right"][0]["data"]["Linked genes"] = f'{link_gene:,}'
    data_summary["Joint"][1]["right"][0]["data"]["Linked peaks"] = f'{link_peak:,}'
    # joint: Joint Cell Calling
    data_summary["Joint"][2]["data"]["data"]["x1"] = atac_summary["atac"]["joint_cell"]["nocell_events"]
    data_summary["Joint"][2]["data"]["data"]["y1"] = atac_summary["atac"]["joint_cell"]["nocell_umi"]
    data_summary["Joint"][2]["data"]["data"]["x2"] = atac_summary["atac"]["joint_cell"]["cell_events"]
    data_summary["Joint"][2]["data"]["data"]["y2"] = atac_summary["atac"]["joint_cell"]["cell_umi"]
    data_summary["Joint"][2]["data"]["data"]["x1num"] = atac_summary["atac"]["joint_cell"]["nocell_num"]
    data_summary["Joint"][2]["data"]["data"]["x2num"] = atac_summary["atac"]["joint_cell"]["cell_num"]
    # joint: left_atac cluster riget_gex cluster
    atacdir = os.path.dirname(atacjson)
    gex_tsnefile = os.path.join(atacdir, 'step4', 'gex_tsne_umi.xls')
    # assert os.path.exists(gex_tsnefile), f'{gex_tsnefile} not found!'
    if not os.path.exists(gex_tsnefile):
        logger.info(f"Warning : The path of '{gex_tsnefile}' is not exists.")
        gex_tsne1 = []
        gex_tsne2 = []
        gex_depth = []
        gex_cluster = []
    else:
        gex_tsne1, gex_tsne2, gex_depth, gex_cluster = get_gex_tsne(gex_tsnefile)

    atac_tsnefile = os.path.join(atacdir, 'step4', 'atac_tsne_umi.xls')
    # assert os.path.exists(atac_tsnefile), f'{atac_tsnefile} not found!'
    if not os.path.exists(atac_tsnefile):
        logger.info(f"Warning : The path of '{atac_tsnefile}' is not exists.")
        atac_tsne1 = []
        atac_tsne2 = []
        atac_depth = []
        atac_cluster = []
    else:
        atac_tsne1, atac_tsne2, atac_depth, atac_cluster = get_atac_tsne(atac_tsnefile)

    data_summary["Joint"][3]["left"][0]["data"]["coordinate"]["tSNE1"] = atac_tsne1
    data_summary["Joint"][3]["left"][0]["data"]["coordinate"]["tSNE2"] = atac_tsne2
    data_summary["Joint"][3]["left"][0]["data"]["cluster"] = atac_cluster
    data_summary["Joint"][3]["right"][0]["data"]["coordinate"]["tSNE1"] = gex_tsne1
    data_summary["Joint"][3]["right"][0]["data"]["coordinate"]["tSNE2"] = gex_tsne2
    data_summary["Joint"][3]["right"][0]["data"]["cluster"] = gex_cluster
    # joint: left_atac depth riget_gex depth
    data_summary["Joint"][4]["left"][0]["data"]["coordinate"]["tSNE1"] = atac_tsne1
    data_summary["Joint"][4]["left"][0]["data"]["coordinate"]["tSNE2"] = atac_tsne2
    data_summary["Joint"][4]["left"][0]["data"]["depth"] = atac_depth
    data_summary["Joint"][4]["right"][0]["data"]["coordinate"]["tSNE1"] = gex_tsne1
    data_summary["Joint"][4]["right"][0]["data"]["coordinate"]["tSNE2"] = gex_tsne2
    data_summary["Joint"][4]["right"][0]["data"]["depth"] = gex_depth

    # rna: title
    data_summary["Rna"][0]["left"][0]["data"]["Estimated number of cells"] = f'{atac_summary["atac"]["Cells"]["Estimated number of cells"]:,}'
    data_summary["Rna"][0]["right"][0]["data"]["GEX Median genes per cell"] = f'{gex_summary["cells"]["Median Genes per Cell"]:,}'
    data_summary["Rna"][0]["right"][0]["data"]["ATAC Median high-quality fragments overlapping peaks per cell"] = f'{atac_summary["atac"]["Cells"]["Median high-quality fragments in peaks per cell"]:,}'
    # rna: left_Sequencing riget_Saturation tu
    data_summary["Rna"][1]["left"][0]["data"]["Number of Reads"] = f'{gex_summary["stat"]["total"]:,}'
    data_summary["Rna"][1]["left"][0]["data"]["Valid Barcode"] = f'{gex_summary["stat"]["valid"]/gex_summary["stat"]["total"]:.2%}'
    data_summary["Rna"][1]["left"][0]["data"]["Sequencing Saturation"] = f'{gex_summary["cells"]["Sequencing Saturation"]:.2%}'
    data_summary["Rna"][1]["left"][0]["data"]["Too Short"] = f'{gex_summary["stat"]["too_short"]:,}'
    b_total_base = sum([sum(v) for v in gex_summary["barcode_q"].values()])
    b30_base = sum([sum(v[30:]) for v in gex_summary["barcode_q"].values()])
    data_summary["Rna"][1]["left"][0]["data"]["Q30 Bases in Barcode"] = f'{b30_base/b_total_base:.2%}'
    u_total_base = sum([sum(v) for v in gex_summary["umi_q"].values()])
    u30_base = sum([sum(v[30:]) for v in gex_summary["umi_q"].values()])
    data_summary["Rna"][1]["left"][0]["data"]["Q30 Bases in UMI"] = f'{u30_base/u_total_base:.2%}'
    # data_summary["Rna"][1]["left"][0]["data"]["Reads with Linker"] = f'{gex_summary["stat"]["seq_17L19ME"]/gex_summary["stat"]["total"]:.2%}'
    data_summary["Rna"][1]["right"][0]["data"]["x"] = [0, ] + gex_summary["downsample"]["Reads"]
    data_summary["Rna"][1]["right"][0]["data"]["y"] = [0, ] + gex_summary["downsample"]["saturation"]
    # rna: left_cells riget_Median Genes tu
    data_summary["Rna"][2]["left"][0]["data"]["Estimated Number of Cells"] = f'{atac_summary["atac"]["Cells"]["Estimated number of cells"]:,}'
    data_summary["Rna"][2]["left"][0]["data"]["Fraction Reads in Cells"] = f'{gex_summary["cells"]["Fraction Reads in Cells"]:.2%}'
    data_summary["Rna"][2]["left"][0]["data"]["Mean Reads per Cell"] = f'{gex_summary["cells"]["Mean Reads per Cell"]:,}'
    data_summary["Rna"][2]["left"][0]["data"]["Median Genes per Cell"] = f'{gex_summary["cells"]["Median Genes per Cell"]:,}'
    data_summary["Rna"][2]["left"][0]["data"]["Median UMI Counts per Cell"] = f'{gex_summary["cells"]["Median UMI Counts per Cell"]:,}'
    data_summary["Rna"][2]["left"][0]["data"]["Total Genes Detected"] = f'{gex_summary["cells"]["Total Genes Detected"]:,}'
    data_summary["Rna"][2]["right"][0]["data"]["x"] = [0, ] + gex_summary["downsample"]["Reads"]
    data_summary["Rna"][2]["right"][0]["data"]["y"] = [0, ] + gex_summary["downsample"]["median"]
    # rna: left_mapping 
    data_summary["Rna"][3]["left"][0]["data"]["Reads Mapped to Genome"] = f'{gex_summary["mapping"]["Reads Mapped to Genome"]:.2%}'
    data_summary["Rna"][3]["left"][0]["data"]["Reads Mapped Confidently to Genome"] = f'{gex_summary["mapping"]["Reads Mapped Confidently to Genome"]:.2%}'
    data_summary["Rna"][3]["left"][0]["data"]["Reads Mapped to Intergenic Regions"] = f'{gex_summary["mapping"]["Reads Mapped to Intergenic Regions"]:.2%}'
    data_summary["Rna"][3]["left"][0]["data"]["Reads Mapped to Intronic Regions"] = f'{gex_summary["mapping"]["Reads Mapped to Intronic Regions"]:.2%}'
    data_summary["Rna"][3]["left"][0]["data"]["Reads Mapped to Exonic Regions"] = f'{gex_summary["mapping"]["Reads Mapped to Exonic Regions"]:.2%}'

    # diff
    diff_table = os.path.join(atacdir, 'step4', 'gex_FindAllMarkers.xls')
    # assert os.path.exists(diff_table), f'{diff_table} not found!'
    if not os.path.exists(diff_table):
        logger.info(f"Warning : The path of '{diff_table}' is not exists.")
    else:
        diff_data = pre_diff_data(diff_table)
        data_summary["diff"] = diff_data

    # atac: title
    data_summary["ATAC"][0]["left"][0]["data"]["Estimated number of cells"] = f'{atac_summary["atac"]["Cells"]["Estimated number of cells"]:,}'
    data_summary["ATAC"][0]["right"][0]["data"]["GEX Median genes per cell"] = f'{gex_summary["cells"]["Median Genes per Cell"]:,}'
    data_summary["ATAC"][0]["right"][0]["data"]["ATAC Median high-quality fragments overlapping peaks per cell"] = f'{atac_summary["atac"]["Cells"]["Median high-quality fragments in peaks per cell"]:,}'
    # atac: left_Sequencing riget_median_fragments tu
    data_summary["ATAC"][1]["left"][0]["data"]["Sequenced read pairs"] = f'{atac_summary["atac"]["Sequencing"]["Sequenced read pairs"]:,}'
    data_summary["ATAC"][1]["left"][0]["data"]["Valid barcodes"] = f'{atac_summary["atac"]["Sequencing"]["Valid barcodes"]:.2%}'
    data_summary["ATAC"][1]["left"][0]["data"]["Too short"] = f'{atac_summary["atac"]["Sequencing"]["Too short"]:,}'
    data_summary["ATAC"][1]["left"][0]["data"]["Q30 bases in barcode"] = f'{atac_summary["atac"]["Sequencing"]["Q30 bases in barcode"]:.2%}'
    data_summary["ATAC"][1]["left"][0]["data"]["Q30 bases in read 1"] = f'{atac_summary["atac"]["Sequencing"]["Q30 bases in read 1"]:.2%}'
    data_summary["ATAC"][1]["left"][0]["data"]["Q30 bases in read 2"] = f'{atac_summary["atac"]["Sequencing"]["Q30 bases in read 2"]:.2%}'
    data_summary["ATAC"][1]["left"][0]["data"]["Percent duplicates"] = f'{atac_summary["atac"]["Sequencing"]["Percent duplicates"]:.2%}'
    data_summary["ATAC"][1]["right"][0]["data"]["x"] = atac_summary["atac"]["median"]["mean_reads"]
    data_summary["ATAC"][1]["right"][0]["data"]["y"] = atac_summary["atac"]["median"]["median_fragments"]
    # atac: left_cells riget_peaks targeting tu
    data_summary["ATAC"][2]["left"][0]["data"]["Estimated number of cells"] = f'{atac_summary["atac"]["Cells"]["Estimated number of cells"]:,}'
    data_summary["ATAC"][2]["left"][0]["data"]["Mean raw read pairs per cell"] = f'{atac_summary["atac"]["Cells"]["Mean raw read pairs per cell"]:,}'
    data_summary["ATAC"][2]["left"][0]["data"]["Fraction of high-quality fragments in cells"] = f'{atac_summary["atac"]["Cells"]["Fraction of high-quality fragments in cells"]:.2%}'
    data_summary["ATAC"][2]["left"][0]["data"]["Fraction of transposition events in peaks in cells"] = f'{atac_summary["atac"]["Cells"]["Fraction of transposition events in peaks in cells"]:.2%}'
    data_summary["ATAC"][2]["left"][0]["data"]["Median high-quality fragments per cell"] = f'{atac_summary["atac"]["Cells"]["Median high-quality fragments per cell"]:,}'
    data_summary["ATAC"][2]["left"][0]["data"]["Median high-quality fragments overlapping peaks per cell"] = f'{atac_summary["atac"]["Cells"]["Median high-quality fragments in peaks per cell"]:,}'
    data_summary["ATAC"][2]["right"][0]["data"]["x1"] = atac_summary["atac"]["peaks_target"]["nocell_fragments"]
    data_summary["ATAC"][2]["right"][0]["data"]["y1"] = atac_summary["atac"]["peaks_target"]["nocell_frac"]
    data_summary["ATAC"][2]["right"][0]["data"]["x2"] = atac_summary["atac"]["peaks_target"]["cell_fragments"]
    data_summary["ATAC"][2]["right"][0]["data"]["y2"] = atac_summary["atac"]["peaks_target"]["cell_frac"]
    # atac: left_targeting riget_tss tu
    data_summary["ATAC"][3]["left"][0]["data"]["Number of peaks"] = f'{atac_summary["atac"]["Targeting"]["Number of peaks"]:,}'
    data_summary["ATAC"][3]["left"][0]["data"]["Fraction of genome in peaks"] = f'{atac_summary["atac"]["Targeting"]["Fraction of genome in peaks"]:.2%}'
    data_summary["ATAC"][3]["left"][0]["data"]["TSS enrichment score"] = f'{atac_summary["atac"]["Targeting"]["TSS enrichment score"]:.2f}'
    data_summary["ATAC"][3]["left"][0]["data"]["Fraction of high-quality fragments overlapping TSS"] = f'{atac_summary["atac"]["Targeting"]["Fraction of high-quality fragments overlapping TSS"]:.2%}'
    data_summary["ATAC"][3]["left"][0]["data"]["Fraction of high-quality fragments overlapping peaks"] = f'{atac_summary["atac"]["Targeting"]["Fraction of high-quality fragments overlapping peaks"]:.2%}'
    data_summary["ATAC"][3]["right"][0]["data"]["x"] = atac_summary["atac"]["tss"]["position"]
    data_summary["ATAC"][3]["right"][0]["data"]["y"] = atac_summary["atac"]["tss"]["score"]
    # atac: left_mapping riget_insert tu
    data_summary["ATAC"][4]["left"][0]["data"]["Confidently mapped read pairs"] = f'{atac_summary["atac"]["Mapping"]["Confidently mapped read pairs"]:.2%}'
    data_summary["ATAC"][4]["left"][0]["data"]["Unmapped read pairs"] = f'{atac_summary["atac"]["Mapping"]["Unmapped read pairs"]:.2%}'
    data_summary["ATAC"][4]["left"][0]["data"]["Non-nuclear read pairs"] = f'{atac_summary["atac"]["Mapping"]["Non-nuclear read pairs"]:.2%}'
    data_summary["ATAC"][4]["right"][0]["data"]["x"] = atac_summary["atac"]["insert"]["size"]
    data_summary["ATAC"][4]["right"][0]["data"]["y"] = atac_summary["atac"]["insert"]["count"]

    # report
    template_dir_new = os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/reportarc'))
    env = Environment(loader=FileSystemLoader(template_dir_new))
    template = env.get_template('base.html')
    with open(os.path.join(outdir, f'{samplename}_report.html'), 'w') as fh:
        fh.write(template.render(websummary_json_data=data_summary))

    # summary csv
    header=('Samplename,'
            'Pipeline_version,'
            'Estimated_number_of_cells,'
            'Feature_linkages_detected,'
            'Linked_genes,'
            'Linked_peaks,'
            'ATAC_Confidently_mapped_read_pairs,'
            'ATAC_Fraction_of_genome_in_peaks,'
            'ATAC_Fraction_of_high-quality_fragments_in_cells,'
            'ATAC_Fraction_of_high-quality_fragments_overlapping_TSS,'
            'ATAC_Fraction_of_high-quality_fragments_overlapping_peaks,'
            'ATAC_Fraction_of_transposition_events_in_peaks_in_cells,'
            'ATAC_Mean_raw_read_pairs_per_cell,'
            'ATAC_Median_high-quality_fragments_per_cell,'
            'ATAC_Non-nuclear_read_pairs,'
            'ATAC_Number_of_peaks,'
            'ATAC_Percent_duplicates,'
            'ATAC_Q30_bases_in_barcode,'
            'ATAC_Q30_bases_in_read_1,'
            'ATAC_Q30_bases_in_read_2,'
            'ATAC_Sequenced_read_pairs,'
            'ATAC_TSS_enrichment_score,'
            'ATAC_Unmapped_read_pairs,'
            'ATAC_Too_Short,'
            'ATAC_Valid_barcodes,'
            'GEX_Fraction_of_reads_in_cells,'
            'GEX_Mean_raw_reads_per_cell,'
            'GEX_Median_UMI_counts_per_cell,'
            'GEX_Median_genes_per_cell,'
            'GEX_Sequencing_Saturation,'
            'GEX_Q30_bases_in_UMI,'
            'GEX_Q30_bases_in_barcode,'
            'GEX_Reads_mapped_confidently_to_exonic_regions,'
            'GEX_Reads_mapped_confidently_to_genome,'
            'GEX_Reads_mapped_confidently_to_intergenic_regions,'
            'GEX_Reads_mapped_confidently_to_intronic_regions,'
            'GEX_Reads_mapped_to_genome,'
            'GEX_Sequenced_read_pairs,'
            'GEX_Too_Short,'
            'GEX_Total_Genes_Detected,'
            'GEX_Valid_barcodes'
            )
    summary_data = [
            samplename,
            "seekarctools " + gex_summary["__version__"],
            atac_summary["atac"]["Cells"]["Estimated number of cells"],
            feature_link,
            link_gene,
            link_peak,
            float(f'{atac_summary["atac"]["Mapping"]["Confidently mapped read pairs"]:.4f}'),
            float(f'{atac_summary["atac"]["Targeting"]["Fraction of genome in peaks"]:.4f}'),
            float(f'{atac_summary["atac"]["Cells"]["Fraction of high-quality fragments in cells"]:.4f}'),
            float(f'{atac_summary["atac"]["Targeting"]["Fraction of high-quality fragments overlapping TSS"]:.4f}'),
            float(f'{atac_summary["atac"]["Targeting"]["Fraction of high-quality fragments overlapping peaks"]:.4f}'),
            float(f'{atac_summary["atac"]["Cells"]["Fraction of transposition events in peaks in cells"]:.4f}'),
            float(f'{atac_summary["atac"]["Cells"]["Mean raw read pairs per cell"]:.4f}'),
            atac_summary["atac"]["Cells"]["Median high-quality fragments per cell"],
            float(f'{atac_summary["atac"]["Mapping"]["Non-nuclear read pairs"]:.4f}'),
            atac_summary["atac"]["Targeting"]["Number of peaks"],
            float(f'{atac_summary["atac"]["Sequencing"]["Percent duplicates"]:.4f}'),
            float(f'{atac_summary["atac"]["Sequencing"]["Q30 bases in barcode"]:.4f}'),
            float(f'{atac_summary["atac"]["Sequencing"]["Q30 bases in read 1"]:.4f}'),
            float(f'{atac_summary["atac"]["Sequencing"]["Q30 bases in read 2"]:.4f}'),
            atac_summary["atac"]["Sequencing"]["Sequenced read pairs"],
            float(f'{atac_summary["atac"]["Targeting"]["TSS enrichment score"]:.4f}'),
            float(f'{atac_summary["atac"]["Mapping"]["Unmapped read pairs"]:.4f}'),
            atac_summary["atac"]["Sequencing"]["Too short"],
            float(f'{atac_summary["atac"]["Sequencing"]["Valid barcodes"]:.4f}'),
            float(f'{gex_summary["cells"]["Fraction Reads in Cells"]:.4f}'),
            float(f'{gex_summary["cells"]["Mean Reads per Cell"]:.4f}'),
            gex_summary["cells"]["Median UMI Counts per Cell"],
            gex_summary["cells"]["Median Genes per Cell"],
            float(f'{gex_summary["cells"]["Sequencing Saturation"]:.4f}'),
            float(f'{u30_base/u_total_base:.4f}'),
            float(f'{b30_base/b_total_base:.4f}'),
            float(f'{gex_summary["mapping"]["Reads Mapped to Exonic Regions"]:.4f}'),
            float(f'{gex_summary["mapping"]["Reads Mapped Confidently to Genome"]:.4f}'),
            float(f'{gex_summary["mapping"]["Reads Mapped to Intergenic Regions"]:.4f}'),
            float(f'{gex_summary["mapping"]["Reads Mapped to Intronic Regions"]:.4f}'),
            float(f'{gex_summary["mapping"]["Reads Mapped to Genome"]:.4f}'),
            gex_summary["stat"]["total"],
            gex_summary["stat"]["too_short"],
            gex_summary["cells"]["Total Genes Detected"],
            float(f'{gex_summary["stat"]["valid"]/gex_summary["stat"]["total"]:.4f}')
        ]
    with open(os.path.join(outdir, f'{samplename}_summary.csv'), 'w') as fh:
        fh.write(header + '\n')
        fh.write(','.join(str(_).replace(',', '') for _ in summary_data)+ '\n')
