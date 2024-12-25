import snapatac2 as snap
import snapatac2._snapatac2
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
from subprocess import run
import os
import json
import gzip
import scipy.io
from ..utils.helper import logger


def makedict(chrlenfile):
    d = {}
    size_list = []
    with open(chrlenfile,'r') as fh:
        for line in fh:
            chrom = str(line.strip().split('\t')[0])
            lenth = int(line.strip().split('\t')[1])
            d[chrom] = lenth
            size_list.append(f'{chrom}:1-{lenth}')
    return d, size_list

def calulate_chunck_size(detail_file):
    chunck_size = 5000000
    with open(detail_file) as f:
        line_count = sum(1 for _ in f)
    n_split = int(line_count/chunck_size)
    if n_split == 0:
        n_split = 1
    chunck_size = int(line_count/n_split) + 1
    return chunck_size


def runpipe(bam:str, atacjson:str, outdir:str, atacname:str, filtercb:str, countxls:str, species:str, refpath:str, 
            bedtoolspath:str="bedtools", gunzippath:str="gunzip", bgzippath:str="bgzip", tabixpath:str="tabix", 
            core:int=4, qvalue:float=0.05, nolambda=False, snapshift:int=0, extsize:int=400, min_len:int=400, blacklist=None, **kwargs):
    step3dir = os.path.join(outdir, "step3")
    os.makedirs(step3dir, exist_ok=True)

    chrNameLength = os.path.join(refpath, 'star/chrNameLength.txt')
    assert os.path.exists(chrNameLength), f'{chrNameLength} not found!'
    gtf = os.path.join(refpath, 'genes/genes.gtf')
    assert os.path.exists(gtf), f'{gtf} not found!'

    qc = {
        "Sequencing": {},
        "Cells": {},
        "Mapping": {},
        "Targeting": {},
        "median": {},
        "peaks_target": {},
        "tss": {},
        "insert": {},
        "joint_cell": {}
    }
    qc["refpath"] = refpath
    qc["Organism"] = species

    # make fragments and QC
    logger.info("make fragments and QC started!")
    fragments_file = os.path.join(step3dir, atacname+"_fragments.tsv.gz")
    bam_qc = snap.pp.make_fragment_file(
        bam_file=bam,
        output_file=fragments_file,
        barcode_regex = "^(.*?)_",
        compression = 'gzip',
        compression_level = 6)

    logger.info("make fragments and QC end!")


    with open(atacjson, "r") as fh:
        atac_summary = json.load(fh)

    qc["Sequencing"]["Sequenced read pairs"] = int(atac_summary["stat"]["total"])
    qc["Sequencing"]["Valid barcodes"] = atac_summary["stat"]["valid"]/atac_summary["stat"]["total"]
    qc["Sequencing"]["Too short"] = atac_summary["stat"]["valid"] - atac_summary["stat"]["step1_available"]
    b_total_base = sum([sum(v) for v in atac_summary["barcode_q"].values()])
    b30_base = sum([sum(v[30:]) for v in atac_summary["barcode_q"].values()])
    qc["Sequencing"]["Q30 bases in barcode"] = b30_base/b_total_base
    qc["Sequencing"]["Q30 bases in read 1"] = bam_qc["frac_q30_bases_read1"]
    qc["Sequencing"]["Q30 bases in read 2"] = bam_qc["frac_q30_bases_read2"]
    qc["Sequencing"]["Percent duplicates"] = bam_qc["frac_duplicates"]
    qc["Mapping"]["Confidently mapped read pairs"] = bam_qc["frac_confidently_mapped"]
    qc["Mapping"]["Unmapped read pairs"] = bam_qc["frac_unmapped"]
    qc["Mapping"]["Non-nuclear read pairs"] = bam_qc["frac_nonnuclear"]

    logger.info(f'Sequencing QC : {qc["Sequencing"]}')
    logger.info(f'Mapping QC : {qc["Mapping"]}')

    # make Anndata obj
    logger.info("make Anndata obj started!")
    size_dict, size_list = makedict(chrNameLength)
    atac = snap.pp.import_data(fragment_file=fragments_file, 
                               chrom_sizes=size_dict,
                               min_num_fragments=0)

    logger.info("make Anndata obj end!")


    # plot fragments size fenbu
    logger.info("plot fragments size started!")
    snap.pl.frag_size_distr(atac, out_file=os.path.join(step3dir, atacname+"_fragments_size.png"))
    # fragments size fenbu data to json
    frag_size_list = list(range(1001))
    frag_count_list = list(atac.uns['frag_size_distr'])
    qc["insert"]["size"] = frag_size_list
    qc["insert"]["count"] = frag_count_list

    logger.info("plot fragments size end!")


    # TSS
    logger.info("plot TSS score started!")
    snap.metrics.tsse(atac, gene_anno=gtf)
    max_value = max(list(atac.uns['TSS_profile'])[999:3000])
    min_value = min(list(atac.uns['TSS_profile'])[999:3000])
    listtss=list(atac.uns['TSS_profile'])[999:3000]
    listtssbz = [x / min_value for x in listtss]
    plt.plot(range(-1000, 1001), listtssbz)
    plt.xlabel('Relative Position (bp from TSS)')
    plt.ylabel('Relative Enrichment')
    plt.savefig(os.path.join(step3dir, atacname+"_tss_enrichment.png"))
    # TSS score fenbu data to json
    tss_position_list = list(range(-1000,1001))
    tss_score_list = listtssbz
    qc["tss"]["position"] = tss_position_list
    qc["tss"]["score"] = tss_score_list

    logger.info("plot TSS score end!")


    # call peak
    logger.info("call peak started!")
    genomelen=atac.uns['reference_sequences']['reference_seq_length'].sum()
    logger.info(f'Parameter : qvalue={qvalue}, nolambda={nolambda}, shift={snapshift}, extsize={extsize}, min_len={min_len}, blacklist={blacklist}, n_jobs={core}')

    snap.tl.macs3(atac, qvalue=qvalue, nolambda=nolambda, shift=snapshift, extsize=extsize, min_len=min_len, blacklist=blacklist, n_jobs=core)
    logger.info(f'Call Peak Done !!!')
    peakfile=os.path.join(step3dir, atacname+"_rawpeaks.bed")
    peaksortu=os.path.join(step3dir, atacname+"_peaksort-u.bed")
    peakuniq=os.path.join(step3dir, atacname+"_peaks.bed")
    with open(peakfile, 'w') as fhout:
        for index, row in atac.uns['macs3_pseudobulk'].iterrows():
            chrom = row['chrom']
            start = row['start']
            end = row['end']
            fhout.write(f'{chrom}\t{start}\t{end}\n')
    # quchong
    cmd = ("less {peakfile} |sort -u > {peaksortu}; "
           "{bedtoolspath} sort -i {peaksortu} > {peakuniq} && rm {peaksortu} {peakfile}"
        ).format(bedtoolspath=bedtoolspath, peakfile=peakfile, peaksortu=peaksortu, peakuniq=peakuniq)
    run(cmd, shell=True)
    # count peaks
    peakuniqlen = 0
    peaknum=0
    peaks = []
    with open(peakuniq, 'r') as fh:
        for line in fh:
            if line.startswith('#'): continue
            peaknum+=1
            tmp = line.strip().split('\t')
            chrom = tmp[0]
            start = int(tmp[1])
            end = int(tmp[2])
            peakuniqlen += end - start
            peaks.append(f'{chrom}:{start}-{end}')

    logger.info(f'snapatac2 call peaks number: {peaknum}')
    logger.info(f'snapatac2 call peaks length: {peakuniqlen}')
    logger.info(f'snapatac2 call peaks fraction: {peakuniqlen/genomelen:.2%}')

    logger.info("call peak end!")


    # ---------------count metrics per barcode---------------
    logger.info("count metrics per barcode started!")
    snapatac2.metrics.frip(atac, {"n_frag_overlap_peak": peaks}, normalized=False)
    snapatac2.metrics.frip(atac, {"events_overlap_peak": peaks}, normalized=False, count_as_insertion=True, inplace=True)
    snapatac2.metrics.frip(atac, {"events_all": size_list}, normalized=False, count_as_insertion=True, inplace=True)
    atac.write(os.path.join(step3dir, atacname+"_snapatac2_raw.h5ad"))
    raw_peak_mat = snap.pp.make_peak_matrix(atac, use_rep=peaks, counting_strategy='insertion')
    # raw_peak_mat.write(os.path.join(step3dir, "raw_peaks_bc_matrix.h5ad"))
    # output raw_peaks_bc_matrix dir
    os.makedirs(os.path.join(step3dir, "raw_peaks_bc_matrix"), exist_ok=True)
    scipy.io.mmwrite(os.path.join(step3dir, "raw_peaks_bc_matrix/matrix.mtx"), raw_peak_mat.X.T.astype(np.float32))
    with open(os.path.join(step3dir, "raw_peaks_bc_matrix/matrix.mtx"), 'rb') as f_in:
        with gzip.open(os.path.join(step3dir, "raw_peaks_bc_matrix/matrix.mtx.gz"), 'wb') as f_out:
            f_out.writelines(f_in)
    with gzip.open(os.path.join(step3dir, "raw_peaks_bc_matrix/features.tsv.gz"), 'wt') as f:
        for feature in raw_peak_mat.var_names:
            f.write(f"{feature}\t{feature}\tpeaks\n")
    with gzip.open(os.path.join(step3dir, "raw_peaks_bc_matrix/barcodes.tsv.gz"), 'wt') as f:
        for cellbarcode in raw_peak_mat.obs_names:
            f.write(f"{cellbarcode}\n")
    os.remove(os.path.join(step3dir, "raw_peaks_bc_matrix/matrix.mtx"))

    logger.info("output raw_peaks_bc_matrix done!")


    # snap fragments number and fragment
    frag_df = pd.DataFrame(atac.obs['n_fragment'])
    frag_df['barcode'] = atac.obs.index

    # fragments overlap peaks
    frag_overpeak_df = pd.DataFrame(atac.obs['n_frag_overlap_peak'])
    frag_overpeak_df['barcode'] = atac.obs.index
    merged_df = pd.merge(frag_df, frag_overpeak_df, on='barcode', how='inner')

    # events overlap peak
    events_overpeak_df = pd.DataFrame(atac.obs['events_overlap_peak'])
    events_overpeak_df['barcode'] = atac.obs.index
    merged_df = pd.merge(merged_df, events_overpeak_df, on='barcode', how='inner')
    # events per barcode
    events_all_df = pd.DataFrame(atac.obs['events_all'])
    events_all_df['barcode'] = atac.obs.index
    merged_df = pd.merge(merged_df, events_all_df, on='barcode', how='inner')

    # atac reads per barcode
    logger.info("read atac step3 fragments.tsv.gz ...")
    d = {}
    with gzip.open(fragments_file, 'rt') as fh, open(os.path.join(step3dir, "frag_counts.xls"), 'w') as fhout:
        fhout.write(f'barcode\tfragment\tnum\treads\n')
        for row in fh:
            tmp = row.strip().split('\t')
            cb = tmp[3]
            fragment = f'{tmp[0]}_{tmp[1]}-{tmp[2]}'
            readsnum = tmp[4]
            d[cb] = d.get(cb, {'fragments_num':0, 'atac_reads':0})
            d[cb]['fragments_num'] += 1
            d[cb]['atac_reads'] += int(readsnum)
            fhout.write(f'{cb}\t{fragment}\t{1}\t{int(readsnum)}\n')
    atac_reads_df = pd.DataFrame.from_dict(d, orient='index')
    atac_reads_df.reset_index(inplace=True)
    atac_reads_df.rename(columns={'index': 'barcode'}, inplace=True)
    merged_df = pd.merge(merged_df, atac_reads_df, on='barcode', how='inner')
    merged_df['fraction_frag_overlap_peak'] = merged_df['n_frag_overlap_peak'] / merged_df['fragments_num']
    logger.info("count atac fragments and reads completed.")


    # ---------------gex gene UMI and reads per barcode---------------
    logger.info("read gex step3 counts.xls...")
    d = {}
    with open(countxls,"rt") as gex_file:
        for line in gex_file:
            if line.startswith("cellID"): continue
            ls = line.strip().split("\t")
            barcode = ls[0]
            gene = ls[1]
            umi = ls[2]
            reads = ls[3]
            d[barcode] = d.get(barcode, {'gex_gene':0, 'gex_umi':0, 'gex_reads':0})
            d[barcode]['gex_gene'] += 1
            d[barcode]['gex_umi'] += int(umi)
            d[barcode]['gex_reads'] += int(reads)
    
    gexdf = pd.DataFrame.from_dict(d, orient='index')
    gexdf.reset_index(inplace=True)
    gexdf.rename(columns={'index': 'barcode'}, inplace=True)

    merged_df = pd.merge(gexdf, merged_df, on='barcode', how='outer')
    merged_df.fillna(0, inplace=True)
    merged_df['gex_gene'] = merged_df['gex_gene'].astype(int)
    merged_df['gex_umi'] = merged_df['gex_umi'].astype(int)
    merged_df['gex_reads'] = merged_df['gex_reads'].astype(int)
    merged_df['n_fragment'] = merged_df['n_fragment'].astype(int)
    merged_df['n_frag_overlap_peak'] = merged_df['n_frag_overlap_peak'].astype(int)
    merged_df['events_overlap_peak'] = merged_df['events_overlap_peak'].astype(int)
    merged_df['events_all'] = merged_df['events_all'].astype(int)
    merged_df['fragments_num'] = merged_df['fragments_num'].astype(int)
    merged_df['atac_reads'] = merged_df['atac_reads'].astype(int)

    logger.info("count metrics per barcode end!")

    # ---------------filter gex cell barcode---------------
    logger.info("joint cell barcode started!")
    joint_cb_list = []
    with gzip.open(filtercb, 'rt') as fh:
        for row in fh:
            joint_cb_list.append(row.strip())
    merged_df['is_cell'] = merged_df['barcode'].isin(joint_cb_list).astype(int)
    merged_df.to_csv(os.path.join(step3dir, "per_barcode_metrics.csv"), index=False)
    cell_merged_df = merged_df[merged_df['is_cell'] == 1]

    logger.info("joint cell barcode end!")


    # ---------------output filter matrix---------------
    filter_atac = atac[joint_cb_list, :]
    filter_atac.write(os.path.join(step3dir, atacname+"_snapatac2_filter.h5ad"))
    filter_peak_mat = snap.pp.make_peak_matrix(filter_atac, use_rep=peaks, counting_strategy='insertion')

    # output filter_peaks_bc_matrix dir
    os.makedirs(os.path.join(step3dir, "filter_peaks_bc_matrix"), exist_ok=True)
    scipy.io.mmwrite(os.path.join(step3dir, "filter_peaks_bc_matrix/matrix.mtx"), filter_peak_mat.X.T.astype(np.float32))
    with open(os.path.join(step3dir, "filter_peaks_bc_matrix/matrix.mtx"), 'rb') as f_in:
        with gzip.open(os.path.join(step3dir, "filter_peaks_bc_matrix/matrix.mtx.gz"), 'wb') as f_out:
            f_out.writelines(f_in)
    with gzip.open(os.path.join(step3dir, "filter_peaks_bc_matrix/features.tsv.gz"), 'wt') as f:
        for feature in filter_peak_mat.var_names:
            f.write(f"{feature}\t{feature}\tpeaks\n")
    with gzip.open(os.path.join(step3dir, "filter_peaks_bc_matrix/barcodes.tsv.gz"), 'wt') as f:
        for cellbarcode in filter_peak_mat.obs_names:
            f.write(f"{cellbarcode}\n")
    os.remove(os.path.join(step3dir, "filter_peaks_bc_matrix/matrix.mtx"))

    logger.info("output filter_peaks_bc_matrix done!")


    # ---------------cell count---------------
    logger.info("cell atac metrics count...")

    n_cells = len(joint_cb_list)
    logger.info(f"Estimated number of cells : {n_cells}")
    logger.info(f'Mean raw read pairs per cell : {int(qc["Sequencing"]["Sequenced read pairs"]/n_cells)}')
    logger.info(f"Fraction of high-quality fragments in cells : {cell_merged_df['fragments_num'].sum()/merged_df['fragments_num'].sum():.2%}")
    logger.info(f"Fraction of transposition events in peaks in cells : {cell_merged_df['events_overlap_peak'].sum() / cell_merged_df['events_all'].sum():.2%}")
    logger.info(f"Median high-quality fragments per cell : {int(cell_merged_df['fragments_num'].median())}")

    qc["Cells"]["Estimated number of cells"] = n_cells
    qc["Cells"]["Mean raw read pairs per cell"] = int(qc["Sequencing"]["Sequenced read pairs"]/n_cells)
    qc["Cells"]["Fraction of high-quality fragments in cells"] = cell_merged_df['fragments_num'].sum()/merged_df['fragments_num'].sum()
    qc["Cells"]["Fraction of transposition events in peaks in cells"] = cell_merged_df['events_overlap_peak'].sum() / cell_merged_df['events_all'].sum()
    qc["Cells"]["Median high-quality fragments per cell"] = int(cell_merged_df['fragments_num'].median())

    # ---------------Targeting Count---------------
    logger.info(f"Number of peaks : {len(peaks)}")
    logger.info(f"Fraction of genome in peaks : {peakuniqlen/genomelen}")
    logger.info(f"TSS enrichment score : {max_value/min_value}")
    logger.info(f"Fraction of high-quality fragments overlapping TSS : {atac.uns['frac_overlap_TSS']:.2%}")
    logger.info(f"Fraction of high-quality fragments overlapping peaks : {cell_merged_df['n_frag_overlap_peak'].sum()/cell_merged_df['fragments_num'].sum():.2%}")

    qc["Targeting"]["Number of peaks"] = len(peaks)
    qc["Targeting"]["Fraction of genome in peaks"] = peakuniqlen/genomelen
    qc["Targeting"]["TSS enrichment score"] = max_value/min_value
    qc["Targeting"]["Fraction of high-quality fragments overlapping TSS"] = atac.uns['frac_overlap_TSS']
    qc["Targeting"]["Fraction of high-quality fragments overlapping peaks"] = cell_merged_df['n_frag_overlap_peak'].sum()/cell_merged_df['fragments_num'].sum()

    logger.info("cell atac metrics count done!")


    # ---------------count medain fragments--------------- 
    logger.info("count medain fragments started!")
    chunck_size = calulate_chunck_size(os.path.join(step3dir, "frag_counts.xls"))
    csv_reader = pd.read_csv(os.path.join(step3dir, "frag_counts.xls"),
                            dtype={
                                "barcode": "category",
                                "fragment": "category",
                                "num": "category",
                                "reads": "int32"
                            },
                            sep="\t",
                            chunksize=chunck_size)

    median_tmp = defaultdict(list)
    cell_reads_total = 0
    for df in csv_reader:
        df = df.loc[df["barcode"].isin(joint_cb_list), :].reset_index(drop=True)
        cell_reads_total += df["reads"].sum()
        rep = df["reads"]
        df = df.drop(["reads"], axis=1)
        idx = df.index.repeat(rep)
        df = df.iloc[idx].reset_index(drop=True)
        del rep, idx
        # shuffle
        df = df.sample(frac=1.0).reset_index(drop=True)
        # downsample
        n_cols_key = [str((i+1)/ 10) for i in range(0,10)]
        for n, interval in enumerate(np.array_split(np.arange(df.shape[0]), 10)):
            idx = interval[-1]
            percentage = n_cols_key[n]
            sampled_df = df.iloc[:idx]
            sampled_df = sampled_df.assign(**{percentage: 1})
            # calculate median for each portion
            median = sampled_df.groupby([sampled_df["barcode"]],observed=True)["fragment"] \
                            .nunique() \
                            .reset_index(drop=True) \
                            .median()
            median_tmp[percentage].append(int(median))

    median_fragments_list = [0]
    for perc, medians in median_tmp.items():
        median = int(sum(medians)/len(medians))
        median_fragments_list.append(median)
    median_fragments_list
    mean_reads_list = [0] + [int(qc["Sequencing"]["Sequenced read pairs"]/n_cells * float(p))for p in n_cols_key]
    percentage_list = ['0'] + [str((i+1)/ 10) for i in range(0,10)]


    qc["median"]["percentage"] = percentage_list
    qc["median"]["mean_reads"] = mean_reads_list
    qc["median"]["median_fragments"] = median_fragments_list

    logger.info("count medain fragments end!")


    # ---------------joint call cell scatter plot, deduplication and density downsampling---------------
    logger.info("downsample started!")
    joint_df = merged_df[['events_overlap_peak', 'gex_umi', 'is_cell']]
    joint_df_unique = joint_df.drop_duplicates(subset=['events_overlap_peak', 'gex_umi'])
    joint_df_sorted = joint_df_unique.sort_values(by='events_overlap_peak')
    x_max = joint_df_sorted['events_overlap_peak'].max()
    bins = np.linspace(0, x_max, num=21)

    sampled_points = []

    for i in range(len(bins) - 1):
        bin_data = joint_df_sorted[(joint_df_sorted['events_overlap_peak'] > bins[i]) & (joint_df_sorted['events_overlap_peak'] <= bins[i + 1])]
        chuqu = int((len(bin_data) / len(joint_df_unique))*2000)
        if len(bin_data) > 50:
            sampled = bin_data.sample(n=chuqu, random_state=42)
        else:
            sampled = bin_data       
        sampled_points.append(sampled)

    final_joint_df = pd.concat(sampled_points)
    cell_umi = final_joint_df.loc[final_joint_df['is_cell'] == 1, 'gex_umi'].tolist()
    cell_events = final_joint_df.loc[final_joint_df['is_cell'] == 1, 'events_overlap_peak'].tolist()
    nocell_umi = final_joint_df.loc[final_joint_df['is_cell'] == 0, 'gex_umi'].tolist()
    nocell_events = final_joint_df.loc[final_joint_df['is_cell'] == 0, 'events_overlap_peak'].tolist()
    qc["joint_cell"]["cell_umi"] = cell_umi # json y2
    qc["joint_cell"]["cell_events"] = cell_events # json x2
    qc["joint_cell"]["nocell_umi"] = nocell_umi # json y1
    qc["joint_cell"]["nocell_events"] = nocell_events # json x1
    qc["joint_cell"]["cell_num"] = (merged_df['is_cell'] == 1).sum() # json x2num
    qc["joint_cell"]["nocell_num"] = (merged_df['is_cell'] == 0).sum() # json x1num


    # ---------------Fragments scatter plot, deduplication and density downsampling---------------
    frag_df = merged_df[['fragments_num', 'fraction_frag_overlap_peak', 'is_cell']]
    frag_df_unique = frag_df.drop_duplicates(subset=['fragments_num', 'fraction_frag_overlap_peak'])
    frag_df_sorted = frag_df_unique.sort_values(by='fragments_num')
    x_max = frag_df_sorted['fragments_num'].max()
    bins = np.linspace(0, x_max, num=21)
    sampled_points = []
    for i in range(len(bins) - 1):
        bin_data = frag_df_sorted[(frag_df_sorted['fragments_num'] > bins[i]) & (frag_df_sorted['fragments_num'] <= bins[i + 1])]
        chuqu = int((len(bin_data) / len(frag_df_unique))*2000)
        if len(bin_data) > 50:
            sampled = bin_data.sample(n=chuqu, random_state=42)
        else:
            sampled = bin_data       
        sampled_points.append(sampled)
    final_frag_df = pd.concat(sampled_points)
    cell_fragments = final_frag_df.loc[final_frag_df['is_cell'] == 1, 'fragments_num'].tolist()
    cell_frac = final_frag_df.loc[final_frag_df['is_cell'] == 1, 'fraction_frag_overlap_peak'].tolist()
    nocell_fragments = final_frag_df.loc[final_frag_df['is_cell'] == 0, 'fragments_num'].tolist()
    nocell_frac = final_frag_df.loc[final_frag_df['is_cell'] == 0, 'fraction_frag_overlap_peak'].tolist()
    qc["peaks_target"]["cell_fragments"] = cell_fragments
    qc["peaks_target"]["cell_frac"] = cell_frac
    qc["peaks_target"]["nocell_fragments"] = nocell_fragments
    qc["peaks_target"]["nocell_frac"] = nocell_frac

    logger.info("downsample end!")

    with open(atacjson, "w") as fh:
        atac_summary["atac"] = qc
        json.dump(
            atac_summary,
            fh,
            indent=4,
            default=lambda o: int(o) if isinstance(o, np.int64) else o
        )

    # sort & index
    logger.info("sort fragments & creat index started!")
    cmd = ("cd {step3dir}; "
           "{gunzippath} {atacname}_fragments.tsv.gz; "
           "{bedtoolspath} sort -i {atacname}_fragments.tsv > {atacname}_fragments_sort.tsv; "
           "{bgzippath} -c {atacname}_fragments_sort.tsv > {atacname}_fragments.tsv.gz; "
           "{tabixpath} -p bed {atacname}_fragments.tsv.gz && rm {atacname}_fragments.tsv {atacname}_fragments_sort.tsv"
        ).format(step3dir=step3dir, atacname=atacname, gunzippath=gunzippath, bedtoolspath=bedtoolspath, bgzippath=bgzippath, tabixpath=tabixpath)
    run(cmd, shell=True)

    logger.info("sort fragments & creat index end!")

