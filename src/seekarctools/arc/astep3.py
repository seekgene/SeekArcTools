import snapatac2 as snap
import snapatac2._snapatac2
import numpy as np
import pandas as pd
from collections import defaultdict
from subprocess import run
import psutil
import os
import json
import gzip
import scipy.io
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestCentroid
from ..utils.helper import logger
from ..utils.wrappers import cmd_execute
from ..utils import countUtil

pd.options.mode.chained_assignment = None 

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
    chunck_size = 10000000
    with open(detail_file) as f:
        line_count = sum(1 for _ in f)
    n_split = int(line_count/chunck_size)
    if n_split == 0:
        n_split = 1
    chunck_size = int(line_count/n_split) + 1
    return chunck_size

def check_fragments(step3dir, atacname):
    fragpath=os.path.join(step3dir, atacname) + "_fragments.tsv.gz"
    fragindex=os.path.join(step3dir, atacname) + "_fragments.tsv.gz.tbi"
    if os.path.exists(fragpath):
        cmd7 = f"rm {fragpath}"
        cmd_execute(cmd7, check=True)
    if os.path.exists(fragindex):
        cmd8 = f"rm {fragindex}"
        cmd_execute(cmd8, check=True)

def log_transform(x, alpha=0.1):
    return np.log(x + alpha)


def runpipe(bam:str, outdir:str, samplename:str, gexoutdir:str, refpath:str, 
            bedtoolspath:str="bedtools", macs3app:str="macs3", sortbedpath:str="sort-bed", gunzippath:str="gunzip", bgzippath:str="bgzip", tabixpath:str="tabix", sortpath:str="sort", 
            core:int=4, qvalue:float=0.05, nolambda=False, snapshift:int=0, extsize:int=400, min_len:int=400, broad=False, broad_cutoff:float=0.1, 
            min_atac_count:int=None, min_gex_count:int=None, retry=False, **kwargs):
    gexstep3dir = os.path.join(gexoutdir, "step3")
    gexname = f"{samplename}_E"
    gexjson = os.path.join(gexoutdir, gexname+"_summary.json")
    atacname = f"{samplename}_A"
    step3dir = os.path.join(outdir, "step3")
    os.makedirs(step3dir, exist_ok=True)

    chrNameLength = os.path.join(refpath, "star", 'chrNameLength.txt')
    assert os.path.exists(chrNameLength), f'{chrNameLength} not found!'
    gtf = os.path.join(refpath, "genes", "genes.gtf")
    assert os.path.exists(gtf), f'{gtf} not found!'

    atacjson = os.path.join(outdir, atacname+"_summary.json")
    with open(atacjson, "r") as fh:
        atac_summary = json.load(fh)

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

    if retry:
        size_dict, size_list = makedict(chrNameLength)
        qc["Sequencing"] = atac_summary["atac"]["Sequencing"]
        qc["Mapping"] = atac_summary["atac"]["Mapping"]

        logger.info("make retry Anndata obj started!")
        raw_snap_h5ad = os.path.join(step3dir, atacname+"_snapatac2_raw.h5ad")
        assert os.path.exists(raw_snap_h5ad), f'{raw_snap_h5ad} not found!'
        atac = sc.read_h5ad(raw_snap_h5ad)
        logger.info("make retry Anndata obj end!")

        # fragment.tsv.gz sort by barcode
        logger.info("fragment.tsv.gz sort by barcode started!")
        fragments_file = os.path.join(outdir, f"../../outs/{atacname}_fragments.tsv.gz")
        cmd = ("cd {step3dir}; "
            "{gunzippath} -c {fragments_file} > {atacname}_fragments.tsv; "
            "{sortpath} -k4,4 {atacname}_fragments.tsv > {atacname}_fragments_sorted_by_barcode.tsv && rm {atacname}_fragments.tsv"
            ).format(step3dir=step3dir, fragments_file=fragments_file, atacname=atacname, gunzippath=gunzippath, sortpath=sortpath)
        logger.info(cmd)
        run(cmd, shell=True)
        sort_fragments_file = os.path.join(step3dir, f"{atacname}_fragments_sorted_by_barcode.tsv")
        logger.info("fragment.tsv.gz sort by barcode end!")

        # atac reads per barcode
        logger.info("read atac step3 fragments.tsv.gz ...")
        d = {}
        d_insert = {}
        fcout = open(os.path.join(step3dir, "frag_counts.xls"), 'w')
        fcout.write(f'barcode\tfragment\tnum\treads\n')
        with open(sort_fragments_file, 'rt') as fh:
            for row in fh:
                tmp = row.strip().split('\t')
                cb = tmp[3]
                fragment = f'{tmp[0]}_{tmp[1]}-{tmp[2]}'
                insertsize = int(tmp[2]) - int(tmp[1])
                if insertsize <= 1000:
                    d_insert[cb] = d_insert.get(cb, {})
                    d_insert[cb][str(insertsize)] = d_insert[cb].get(str(insertsize), 0)
                    d_insert[cb][str(insertsize)] += 1
                readsnum = tmp[4]
                fcout.write(f'{cb}\t{fragment}\t{1}\t{int(readsnum)}\n')
                d[cb] = d.get(cb, {'fragments_num':0, 'atac_reads':0})
                d[cb]['fragments_num'] += 1
                d[cb]['atac_reads'] += int(readsnum)
        fcout.close()
        atac_reads_df = pd.DataFrame.from_dict(d, orient='index')
        atac_reads_df.reset_index(inplace=True)
        atac_reads_df.rename(columns={'index': 'barcode'}, inplace=True)
        os.remove(sort_fragments_file)
        logger.info("retry fragments count end!")
    else:
        # make fragments and QC
        logger.info("make fragments and QC started!")
        raw_fragments_file = os.path.join(step3dir, atacname+"_raw_fragments.tsv.gz")
        bam_qc = snap.pp.make_fragment_file(
            bam_file=bam,
            output_file=raw_fragments_file,
            barcode_regex = "^(.*?)_",
            compression = 'gzip',
            compression_level = 6)
        logger.info("make fragments and QC end!")

        qc["Sequencing"]["Sequenced read pairs"] = int(atac_summary["stat"]["total"])
        qc["Sequencing"]["Valid barcodes"] = atac_summary["stat"]["valid"]/atac_summary["stat"]["total"]
        qc["Sequencing"]["Too short"] = atac_summary["stat"]["too_short"]
        b_total_base = sum([sum(v) for v in atac_summary["barcode_q"].values()])
        b30_base = sum([sum(v[30:]) for v in atac_summary["barcode_q"].values()])
        qc["Sequencing"]["Q30 bases in barcode"] = b30_base/b_total_base
        qc["Sequencing"]["Q30 bases in read 1"] = bam_qc["frac_q30_bases_read1"]
        qc["Sequencing"]["Q30 bases in read 2"] = bam_qc["frac_q30_bases_read2"]
        qc["Sequencing"]["Percent duplicates"] = bam_qc["frac_duplicates"]
        qc["Mapping"]["Confidently mapped read pairs"] = atac_summary["bam_count"]["confidently_mapped_read_pairs"]
        qc["Mapping"]["Unmapped read pairs"] = atac_summary["bam_count"]["unmapped_read_pairs"]
        qc["Mapping"]["Non-nuclear read pairs"] = bam_qc["frac_nonnuclear"]

        logger.info("dedup raw fragments file & del chrM & count started!")
        fragments_file = os.path.join(step3dir, atacname+"_fragments.tsv.gz")
        d = {}
        d_insert = {}
        fcout = open(os.path.join(step3dir, "frag_counts.xls"), 'w')
        fcout.write(f'barcode\tfragment\tnum\treads\n')
        with gzip.open(raw_fragments_file, 'rt') as fh, gzip.open(fragments_file, "wt") as fhout:
            prev_key = None
            total_reads = 0
            for line in fh:
                fields = line.strip().split("\t")
                chrom, start, end, barcode, reads_count = fields
                if chrom in ['chrM','M']:
                    continue
                key = (chrom, start, end, barcode)            
                if key == prev_key:
                    total_reads += int(reads_count)
                else:
                    if prev_key is not None:
                        chrom_prev, start_prev, end_prev, barcode_prev = prev_key
                        fhout.write(f"{chrom_prev}\t{start_prev}\t{end_prev}\t{barcode_prev}\t{total_reads}\n")
                        cb = barcode_prev
                        fragment = f'{chrom_prev}_{start_prev}-{end_prev}'
                        insertsize = int(end_prev) - int(start_prev)
                        if insertsize <= 1000:
                            d_insert[cb] = d_insert.get(cb, {})
                            d_insert[cb][str(insertsize)] = d_insert[cb].get(str(insertsize), 0)
                            d_insert[cb][str(insertsize)] += 1
                        fcout.write(f'{cb}\t{fragment}\t{1}\t{int(total_reads)}\n')
                        d[cb] = d.get(cb, {'fragments_num':0, 'atac_reads':0})
                        d[cb]['fragments_num'] += 1
                        d[cb]['atac_reads'] += int(total_reads)
                    prev_key = key
                    total_reads = int(reads_count)
            
            if prev_key is not None:
                chrom_prev, start_prev, end_prev, barcode_prev = prev_key
                fhout.write(f"{chrom_prev}\t{start_prev}\t{end_prev}\t{barcode_prev}\t{total_reads}\n")
                cb = barcode_prev
                fragment = f'{chrom_prev}_{start_prev}-{end_prev}'
                insertsize = int(end_prev) - int(start_prev)
                if insertsize <= 1000:
                    d_insert[cb] = d_insert.get(cb, {})
                    d_insert[cb][str(insertsize)] = d_insert[cb].get(str(insertsize), 0)
                    d_insert[cb][str(insertsize)] += 1
                fcout.write(f'{cb}\t{fragment}\t{1}\t{int(total_reads)}\n')
                d[cb] = d.get(cb, {'fragments_num':0, 'atac_reads':0})
                d[cb]['fragments_num'] += 1
                d[cb]['atac_reads'] += int(total_reads)
        fcout.close()

        atac_reads_df = pd.DataFrame.from_dict(d, orient='index')
        atac_reads_df.reset_index(inplace=True)
        atac_reads_df.rename(columns={'index': 'barcode'}, inplace=True)
        logger.info("dedup raw fragments file & del chrM & count end!")


        # make Anndata obj
        logger.info("make Anndata obj started!")
        size_dict, size_list = makedict(chrNameLength)
        atac = snap.pp.import_data(fragment_file=fragments_file, 
                                chrom_sizes=size_dict,
                                min_num_fragments=0)
        # atac.write(os.path.join(step3dir, atacname+"_snapatac2_raw.h5ad"))
        os.makedirs('/tmp', exist_ok=True)
        atac.write(os.path.join('/tmp', atacname+"_snapatac2_raw.h5ad"))
        tmp_rawh5 = os.path.join('/tmp', atacname+"_snapatac2_raw.h5ad")
        rawh5 = os.path.join(step3dir, atacname+"_snapatac2_raw.h5ad")
        cmd = (f"mv {tmp_rawh5} {rawh5}")
        run(cmd, shell=True)

        logger.info("make Anndata obj end!")

    # TSS
    logger.info("plot TSS score started!")
    snap.metrics.tsse(atac, gene_anno=gtf)
    max_value = max(list(atac.uns['TSS_profile'])[999:3000])
    min_value = min(list(atac.uns['TSS_profile'])[999:3000])
    listtss=list(atac.uns['TSS_profile'])[999:3000]
    listtssbz = [x / min_value for x in listtss]

    # TSS score fenbu data to json
    tss_position_list = list(range(-1000,1001))
    tss_score_list = listtssbz
    qc["tss"]["position"] = tss_position_list
    qc["tss"]["score"] = tss_score_list

    logger.info("plot TSS score end!")


    # call peak
    logger.info("call peak started!")
    genomelen=atac.uns['reference_sequences']['reference_seq_length'].sum()
    fragments_bedgz = os.path.join(step3dir, atacname+"_fragments.bed.gz")
    fragments_bed = os.path.join(step3dir, atacname+"_fragments.bed")
    cmd = ("cp {fragments_file} {fragments_bedgz}; {gunzippath} {fragments_bedgz}"
        ).format(fragments_file=fragments_file, fragments_bedgz=fragments_bedgz, gunzippath=gunzippath)
    run(cmd, shell=True)
    logger.info(f'Parameter : qvalue={qvalue}, nolambda={nolambda}, shift={snapshift}, extsize={extsize}, min_len={min_len}, call_broad_peaks={broad}, broad_cutoff={broad_cutoff}, n_jobs={core}')
    if broad:
        cmd = ("{macs3app} callpeak -t {fragments_bed} -f BED -g {genomelen} --broad -n {atacname} --outdir {step3dir} -q {qvalue} --broad-cutoff {broad_cutoff} --shift {snapshift} --extsize {extsize} --min-length {min_len} --nomodel"
            ).format(macs3app=macs3app, fragments_bed=fragments_bed, genomelen=genomelen, atacname=atacname, step3dir=step3dir, qvalue=qvalue, broad_cutoff=broad_cutoff, snapshift=snapshift, extsize=extsize, min_len=min_len)
        run(cmd, shell=True)
        broad_peaks_file = os.path.join(step3dir, atacname+"_peaks.broadPeak")
        broad_peaks_df = pd.read_csv(
            broad_peaks_file, 
            sep='\t', 
            header=None, 
            names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'fold_change', 'pvalue', 'qvalue'])
        atac.uns['macs3_pseudobulk'] = broad_peaks_df
        os.remove(broad_peaks_file)
        gapped_peaks_file = os.path.join(step3dir, atacname+"_peaks.gappedPeak")
        os.remove(gapped_peaks_file)
    else:
        snap.tl.macs3(atac, qvalue=qvalue, nolambda=nolambda, shift=snapshift, extsize=extsize, min_len=min_len, blacklist=None, n_jobs=core)
    logger.info(f'Call Peak Done !!!')
    peakfile=os.path.join(step3dir, atacname+"_rawpeaks.bed")
    peakuniq=os.path.join(step3dir, atacname+"_peaks.bed")
    with open(peakfile, 'w') as fhout:
        for index, row in atac.uns['macs3_pseudobulk'].iterrows():
            chrom = row['chrom']
            start = row['start']
            end = row['end']
            if start == 0:
                start = 1
            if end > size_dict[chrom]:
                end = size_dict[chrom] - 1
            fhout.write(f'{chrom}\t{start}\t{end}\n')
    # Remove duplicate lines
    cmd = ("{sortbedpath} --max-mem 1G --unique {peakfile} > {peakuniq} && rm {peakfile}"
        ).format(sortbedpath=sortbedpath, peakfile=peakfile, peakuniq=peakuniq)
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
    genome_in_peaks = peakuniqlen/genomelen

    logger.info("call peak end!")


    fragments_inpeakbed = os.path.join(step3dir, atacname+"_fragments_inPeaks.bed")
    fragments_inpeak_uniq = os.path.join(step3dir, atacname+"_fragments_inPeaks_uniq.bed")
    cmd = ("{bedtoolspath} intersect -a {fragments_bed} -b {peakuniq} -wa > {fragments_inpeakbed}; "
           "{sortbedpath} --max-mem 4G --unique {fragments_inpeakbed} > {fragments_inpeak_uniq}"
        ).format(bedtoolspath=bedtoolspath, fragments_bed=fragments_bed, peakuniq=peakuniq, fragments_inpeakbed=fragments_inpeakbed, sortbedpath=sortbedpath, fragments_inpeak_uniq=fragments_inpeak_uniq)
    run(cmd, shell=True)

    dtmp1 = {}
    with open(fragments_inpeak_uniq, 'rt') as fh:
        for row in fh:
            tmp = row.strip().split('\t')
            cb = tmp[3]
            dtmp1[cb] = dtmp1.get(cb, 0)
            dtmp1[cb] += 1
            fragment = f'{tmp[0]}_{tmp[1]}-{tmp[2]}'
            readsnum = tmp[4]

    dftmp1 = pd.DataFrame(list(dtmp1.items()), columns=["barcode", "fragments_num_overlap_peaks"])


    # ---------------count metrics per barcode---------------
    logger.info("count metrics per barcode started!")
    snapatac2.metrics.frip(atac, {"n_frag_overlap_peak": peaks}, normalized=False)
    snapatac2.metrics.frip(atac, {"events_overlap_peak": peaks}, normalized=False, count_as_insertion=True, inplace=True)
    snapatac2.metrics.frip(atac, {"events_all": size_list}, normalized=False, count_as_insertion=True, inplace=True)

    raw_peak_mat = snap.pp.make_peak_matrix(atac, use_rep=peaks, counting_strategy='insertion')
    # output raw_peaks_bc_matrix dir
    os.makedirs(os.path.join(step3dir, "raw_peaks_bc_matrix"), exist_ok=True)
    scipy.io.mmwrite(os.path.join(step3dir, "raw_peaks_bc_matrix/matrix.mtx"), raw_peak_mat.X.T.astype(np.int32))
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

    merged_df = pd.merge(merged_df, atac_reads_df, on='barcode', how='inner')
    merged_df = pd.merge(merged_df, dftmp1, on='barcode', how='left')
    merged_df.fillna(0, inplace=True)
    merged_df['fraction_frag_overlap_peak'] = merged_df['fragments_num_overlap_peaks'] / merged_df['fragments_num']
    
    logger.info("count atac fragments and reads completed.")


    # ---------------gex gene UMI and reads per barcode---------------
    logger.info("read gex step3 counts.xls...")
    d = {}
    countxls = os.path.join(gexstep3dir, "counts.xls")
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
    logger.info("joint cell barcode start...")

    if min_atac_count is not None and min_gex_count is not None:
        logger.info("Call cells through min_atac_comunt and min_gex_comunt ...")
        logger.info(f"min_atac_comunt : {int(min_atac_count)}")
        logger.info(f"min_gex_comunt : {int(min_gex_count)}")
        joint_df_cell = merged_df[(merged_df['gex_umi'] >= int(min_gex_count)) & (merged_df['events_overlap_peak'] >= int(min_atac_count))]
        joint_cb_list = list(joint_df_cell['barcode'])

    else:
        df_filter3 = merged_df[(merged_df['gex_umi'] > 1) & (merged_df['events_overlap_peak'] > 1)]
        df_filter1 = df_filter3[df_filter3['fraction_frag_overlap_peak'] > genome_in_peaks]
        joint_df = df_filter1[['barcode', 'gex_umi', 'events_overlap_peak']]
        joint_df_unique = joint_df.drop_duplicates(subset=['gex_umi', 'events_overlap_peak'])
        points = log_transform(joint_df_unique[['gex_umi', 'events_overlap_peak']].values)

        kmeans = KMeans(n_clusters=2, init='k-means++', n_init=10, max_iter=300)
        kmeans.fit(points)
        labels = kmeans.labels_
        centroids = kmeans.cluster_centers_
        if centroids[0].mean() > centroids[1].mean():
            cell_label = 0
            nocell_label = 1
        else:
            cell_label = 1
            nocell_label = 0
        joint_df_unique['kmeans_labels'] = labels
        joint_df_unique['is_cell'] = (joint_df_unique['kmeans_labels'] == cell_label).astype(int)
        joint_df_uniquecell = joint_df_unique[joint_df_unique['is_cell'] == 1]
        joint_cb_list = list(joint_df_uniquecell['barcode'])

    merged_df['is_cell'] = np.where(merged_df['barcode'].isin(joint_cb_list), 1, 0)
    merged_df.drop('n_fragment', axis=1, inplace=True)
    merged_df.drop('n_frag_overlap_peak', axis=1, inplace=True)
    merged_df.to_csv(os.path.join(step3dir, "per_barcode_metrics.csv"), index=False)
    cell_merged_df = merged_df[merged_df['is_cell'] == 1]

    logger.info("joint cell barcode end!")


    # plot fragments size fenbu
    logger.info("count fragments size start...")
    max_insert = 1000
    joint_cb_list_set = set(joint_cb_list)
    insert_counts = {str(i): 0 for i in range(max_insert + 1)}
    for barcode, insert_dict in d_insert.items():
        if barcode in joint_cb_list_set:
            for insert_len, count in insert_dict.items():
                insert_counts[insert_len] += count
    insert_size_list = [i for i in range(max_insert + 1)]
    insert_counts_list = list(insert_counts.values())
    # fragments size fenbu data to json
    qc["insert"]["size"] = insert_size_list
    qc["insert"]["count"] = insert_counts_list

    logger.info("count fragments size end!")
    

    # ---------------output filtered_peaks_bc_matrix dir---------------
    logger.info("output filtered_feature_bc_matrix start...")
    gex_adata = sc.read_10x_mtx(os.path.join(gexstep3dir, "raw_feature_bc_matrix"), var_names='gene_ids')

    gex_cell_adata = gex_adata[joint_cb_list, :]
    dict_idgene = gex_cell_adata.var['gene_symbols'].to_dict()
    os.makedirs(os.path.join(gexstep3dir, "filtered_feature_bc_matrix"), exist_ok=True)

    scipy.io.mmwrite(os.path.join(gexstep3dir, "filtered_feature_bc_matrix/matrix.mtx"), gex_cell_adata.X.T.astype(np.int32))
    with open(os.path.join(gexstep3dir, "filtered_feature_bc_matrix/matrix.mtx"), 'rb') as f_in:
        with gzip.open(os.path.join(gexstep3dir, "filtered_feature_bc_matrix/matrix.mtx.gz"), 'wb') as f_out:
            f_out.writelines(f_in)
    with gzip.open(os.path.join(gexstep3dir, "filtered_feature_bc_matrix/features.tsv.gz"), 'wt') as f:
        for feature in gex_cell_adata.var_names:
            f.write(f"{feature}\t{dict_idgene[feature]}\tGene Expression\n")
    with gzip.open(os.path.join(gexstep3dir, "filtered_feature_bc_matrix/barcodes.tsv.gz"), 'wt') as f:
        for cellbarcode in gex_cell_adata.obs_names:
            f.write(f"{cellbarcode}\n")
    os.remove(os.path.join(gexstep3dir, "filtered_feature_bc_matrix/matrix.mtx"))
    logger.info("output filtered_feature_bc_matrix done!")


    # ---------------output filtered_peaks_bc_matrix---------------
    logger.info("output filtered_peaks_bc_matrix start...")
    filter_atac = atac[joint_cb_list, :]
    # filter_atac.write(os.path.join(step3dir, atacname+"_snapatac2_filter.h5ad"))
    filter_peak_mat = snap.pp.make_peak_matrix(filter_atac, use_rep=peaks, counting_strategy='insertion')

    # output filtered_peaks_bc_matrix dir
    os.makedirs(os.path.join(step3dir, "filtered_peaks_bc_matrix"), exist_ok=True)
    scipy.io.mmwrite(os.path.join(step3dir, "filtered_peaks_bc_matrix/matrix.mtx"), filter_peak_mat.X.T.astype(np.int32))
    with open(os.path.join(step3dir, "filtered_peaks_bc_matrix/matrix.mtx"), 'rb') as f_in:
        with gzip.open(os.path.join(step3dir, "filtered_peaks_bc_matrix/matrix.mtx.gz"), 'wb') as f_out:
            f_out.writelines(f_in)
    with gzip.open(os.path.join(step3dir, "filtered_peaks_bc_matrix/features.tsv.gz"), 'wt') as f:
        for feature in filter_peak_mat.var_names:
            f.write(f"{feature}\t{feature}\tpeaks\n")
    with gzip.open(os.path.join(step3dir, "filtered_peaks_bc_matrix/barcodes.tsv.gz"), 'wt') as f:
        for cellbarcode in filter_peak_mat.obs_names:
            f.write(f"{cellbarcode}\n")
    os.remove(os.path.join(step3dir, "filtered_peaks_bc_matrix/matrix.mtx"))

    logger.info("output filtered_peaks_bc_matrix done!")


    # ---------------atac cell count---------------
    logger.info("cell atac metrics count...")

    n_cells = len(joint_cb_list)
    qc["Cells"]["Estimated number of cells"] = n_cells
    qc["Cells"]["Mean raw read pairs per cell"] = int(qc["Sequencing"]["Sequenced read pairs"]/n_cells)
    qc["Cells"]["Fraction of high-quality fragments in cells"] = cell_merged_df['fragments_num'].sum()/merged_df['fragments_num'].sum()
    qc["Cells"]["Fraction of transposition events in peaks in cells"] = cell_merged_df['events_overlap_peak'].sum() / cell_merged_df['events_all'].sum()
    qc["Cells"]["Median high-quality fragments per cell"] = int(cell_merged_df['fragments_num'].median())
    qc["Cells"]["Median high-quality fragments in peaks per cell"] = int(cell_merged_df['fragments_num_overlap_peaks'].median())

    # ---------------atac Targeting Count---------------
    qc["Targeting"]["Number of peaks"] = len(peaks)
    qc["Targeting"]["Fraction of genome in peaks"] = peakuniqlen/genomelen
    qc["Targeting"]["TSS enrichment score"] = max_value/min_value
    qc["Targeting"]["Fraction of high-quality fragments overlapping TSS"] = atac.uns['frac_overlap_TSS']
    qc["Targeting"]["Fraction of high-quality fragments overlapping peaks"] = cell_merged_df['fragments_num_overlap_peaks'].sum()/cell_merged_df['fragments_num'].sum()

    logger.info("cell atac metrics count done!")


    # ---------------gex downsample---------------------------
    logger.info("gex calculate metrics start...")
    summary_tmp, downsample = countUtil.calculate_metrics(
        counts_file = os.path.join(gexstep3dir, "counts.xls"),
        detail_file = os.path.join(gexstep3dir, "detail.xls"),
        filterd_barcodes_file = os.path.join(gexstep3dir, "filtered_feature_bc_matrix/barcodes.tsv.gz"),
        filterd_features_file = os.path.join(gexstep3dir, "filtered_feature_bc_matrix/features.tsv.gz"),
        gtf = gtf,
        basedir = gexstep3dir
    )
    with open(gexjson, "r") as fh:
        summary = json.load(fh)
        Total = summary["stat"]["total"]
    with open(gexjson, "w") as fh:
        estimated_cell_num = summary_tmp["Estimated Number of Cells"]
        mean_reads_per_cell = Total / estimated_cell_num
        summary_tmp["Mean Reads per Cell"] = int(mean_reads_per_cell)
        del summary_tmp["Median lnc Genes per Cell"]
        summary["cells"] = summary_tmp
        downsample["Reads"] = [
            int(summary_tmp["Mean Reads per Cell"] * p)
            for p in downsample["percentage"]
        ]
        summary["downsample"] = downsample
        json.dump(
            summary,
            fh,
            indent=4,
            default=lambda o: int(o) if isinstance(o, np.int64) else o
        )

    logger.info("gex calculate metrics done!")

    # ---------------atac count medain fragments--------------- 
    logger.info("count medain fragments started!")
    chunck_size = calulate_chunck_size(os.path.join(step3dir, "frag_counts.xls"))
    csv_reader = pd.read_csv(os.path.join(step3dir, "frag_counts.xls"),
                            dtype={
                                "barcode": "category",
                                "fragment": "category",
                                "num": "int32",
                                "reads": "int32"
                            },
                            sep="\t",
                            chunksize=chunck_size)

    median_tmp = defaultdict(list)
    chunck_count = 0
    for df in csv_reader:
        df = df.loc[df["barcode"].isin(joint_cb_list), :].reset_index(drop=True)
        rep = df["reads"]
        df = df.drop(["reads"], axis=1)
        idx = df.index.repeat(rep)
        df = df.iloc[idx].reset_index(drop=True)
        del rep, idx
        # shuffle
        df = df.sample(frac=1.0).reset_index(drop=True)
        # downsample
        n_cols_key = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        for n, interval in enumerate(np.array_split(np.arange(df.shape[0]), 10)):
            idx = interval[-1]
            percentage = str(n_cols_key[n])
            sampled_df = df.iloc[:idx]
            sampled_df = sampled_df.assign(**{percentage: 1})
            median = sampled_df.groupby([sampled_df["barcode"]],observed=True)["fragment"] \
                            .nunique() \
                            .reset_index(drop=True) \
                            .median()
            median_tmp[percentage].append(median)
        chunck_count += 1

    median_sampling = []
    for perc, medians in median_tmp.items():
        median = int(np.mean(medians))
        median_sampling.append(median)

    frac_1_median = int(median_sampling[-1])
    median_frag = int(cell_merged_df['fragments_num'].median())
    k = median_frag / frac_1_median
    median_downsample = [int(x * k) for x in median_sampling]
    median_downsample.pop()
    median_fragments_list = [0] + median_downsample + [median_frag]
    mean_reads_list = [0] + [int(qc["Sequencing"]["Sequenced read pairs"]/n_cells * float(p))for p in n_cols_key]
    percentage_list = [0] + n_cols_key

    qc["median"]["percentage"] = percentage_list
    qc["median"]["mean_reads"] = mean_reads_list
    qc["median"]["median_fragments"] = median_fragments_list
    logger.info("count medain fragments end!")


    # ---------------joint call cell scatter plot, deduplication and density downsampling---------------
    logger.info("downsample started!")
    joint_df = merged_df[['barcode', 'events_overlap_peak', 'gex_umi', 'is_cell']]
    joint_df_unique = joint_df.drop_duplicates(subset=['events_overlap_peak', 'gex_umi'])

    cell_merged_df = joint_df_unique[joint_df_unique['is_cell'] == 1]
    if len(cell_merged_df['barcode']) <= 4000:
        cell_merged_df_4000 = cell_merged_df
    else:
        cell_merged_df_4000 = cell_merged_df.sample(n=4000, random_state=42)
        
    nocell_merged_df = joint_df_unique[joint_df_unique['is_cell'] == 0]
    nocell_merged_df_sorted = nocell_merged_df.sort_values(by='events_overlap_peak')
    x_max = nocell_merged_df_sorted['events_overlap_peak'].max()
    log_x_max = np.log10(x_max)
    logbins = np.linspace(0, log_x_max, num=21)
    bins = [10 ** v for v in logbins]

    sampled_points = []

    for i in range(len(bins) - 1):
        bin_data = nocell_merged_df_sorted[(nocell_merged_df_sorted['events_overlap_peak'] >= bins[i]) & (nocell_merged_df_sorted['events_overlap_peak'] < bins[i + 1])]
        chuqu = int((len(bin_data) / len(nocell_merged_df))*2000)
        if len(bin_data) > 100:
            if chuqu < len(bin_data):
                sampled = bin_data.sample(n=chuqu, random_state=42)
            else:
                sampled = bin_data
        else:
            sampled = bin_data   
        sampled_points.append(sampled)

    nocell_merged_df2000 = pd.concat(sampled_points)
    final_joint_df = pd.concat([cell_merged_df_4000, nocell_merged_df2000], ignore_index=True)
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
    frag_df = merged_df[['barcode', 'fragments_num', 'fraction_frag_overlap_peak', 'is_cell']]
    final_frag_df = frag_df[frag_df['barcode'].isin(final_joint_df['barcode'])]
    
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

    if retry:
        logger.info("do not need sort fragments & creat index!")
    else:
        # sort & index
        logger.info("sort fragments & creat index started!")
        available_memory = psutil.virtual_memory().available / (1024 * 1024 * 1024)
        memory = int(available_memory * 0.9)
        cmd = ("cd {step3dir}; "
            "{gunzippath} {atacname}_fragments.tsv.gz; "
            "{sortbedpath} --max-mem {memory}G {atacname}_fragments.tsv > {atacname}_fragments_sort.tsv; "
            "{bgzippath} -c {atacname}_fragments_sort.tsv > {atacname}_fragments.tsv.gz; "
            "{tabixpath} -p bed {atacname}_fragments.tsv.gz && rm {atacname}_fragments.tsv {atacname}_fragments_sort.tsv {atacname}_raw_fragments.tsv.gz"
            ).format(step3dir=step3dir, atacname=atacname, memory=memory, gunzippath=gunzippath, sortbedpath=sortbedpath, bgzippath=bgzippath, tabixpath=tabixpath)
        run(cmd, shell=True)
        logger.info("sort fragments & creat index end!")

    os.remove(fragments_bed)
    os.remove(fragments_inpeakbed)
    os.remove(fragments_inpeak_uniq)

