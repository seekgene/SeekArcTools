import os
import sys
import json
import dnaio
import itertools
from functools import partial
from .pipeline import Pipeline
from collections import defaultdict, Counter
from .helper import AdapterFilter, QcStat, parse_structure, read_file, get_new_bc, logger
from .chemistry import ATAC_R1_MINLEN, ATAC_R2_MINLEN, CHEMISTRY
from .wrappers import cmd_execute

def barcode_report(logfile):
    with open(logfile) as fh:
        summary = json.load(fh)
    d = summary["stat"]
    return float(d["valid"])/d["total"]

def cut_fq(fq:str, outdir:str, reads_num:int=400000):
    fq_tmp = os.path.join(outdir, os.path.basename(fq))
    cmd = f'zcat "{fq}"|head -n {reads_num}|gzip > {fq_tmp}'
    cmd_execute(cmd, check=True)
    return fq_tmp

def chemistry_auto(fq1, fq2, outdir, use_short_read=False, **kwargs):
    """"
    test with 1M reads
    """
    _outdir = f"{outdir}/.test"
    rawdata = f"{_outdir}/data"
    os.makedirs(rawdata, exist_ok=True)
    fq1_1M, fq2_1M = cut_fq(fq1[0], rawdata), cut_fq(fq2[0], rawdata)

    rate_dict = {}
    
    if kwargs["chemistry"] == "Auto":
        test_chems = ["MM", "DDV2"]
    elif kwargs["chemistry"] == "custom":
        CHEMISTRY["custom"] = {
            "shift": kwargs["shift"],
            "shift_pattern": kwargs["shift_pattern"],
            "structure": kwargs["structure"],
            "barcode": kwargs["barcode"],
            "linker": kwargs["linker"],
            "match_type": [1,],           
        }
        test_chems = ["custom",]
    else:
        test_chems = [kwargs["chemistry"],]

    for chem in test_chems:
        logger.info(f"test {chem}!")
        barcode_main([fq1_1M,], [fq2_1M,], chem, f"{_outdir}/{chem}", core=1, use_short_read=use_short_read, **CHEMISTRY[chem])
        rate = barcode_report(f"{_outdir}/{chem}/{chem}_A_summary.json")
        rate_dict[chem] = rate
    for k, v in rate_dict.items():
        logger.info(f"valid barcode rate of {k}: {v*100:.3f}%")
    if kwargs["chemistry"] == "Auto":
        if rate_dict["DDV2"] > rate_dict["MM"]:
            chemistry = "DDV2"
        else:
            chemistry = "MM"
    else:
        chemistry = kwargs["chemistry"]

    return CHEMISTRY[chemistry]

def check_atac_options(**kwargs):
    kwargs = chemistry_auto(**kwargs)
    return kwargs

def summary(seq, seq_q, seq_dict, qua_dict):
    
    for i, (base, q) in enumerate(zip(seq, seq_q)):
        seq_dict[(i,base)] += 1
        qua_dict[(i,q)] += 1
        
    return seq_dict, qua_dict

def process_barcode(fq1, fq2, fq_out, fqout_multi, r1_structure, shift, shift_pattern,
                    barcode_wl_dict, linker_wl_dict, match_type_dict, adapter1=[["AAAAAAAAAAAA", "3"],],
                    adapter2=[["AAAAAAAAAAAA", "3"],], do_B_correction=True, do_L_correction=True,
                    use_multi=True, use_short_read=False, paired_out=True):
    
    barcode_list_flag = False
    linker_list_flag = False
    if len(barcode_wl_dict)>0:
        barcode_list_flag = True

    if len(linker_wl_dict)>0:
        linker_list_flag = True

    stat_Dict = defaultdict(int)
    Barcode_Counter = Counter()
    
    Barcode_GC_Counter = Counter()
    UMI_GC_Counter = Counter()
    R2_GC_Counter = Counter()
    Barcode_Q_Counter = Counter()
    UMI_Q_Counter = Counter()
    R2_Q_Counter = Counter()
    
    adapter_filter = AdapterFilter(adapter1=adapter1, adapter2=adapter2)
    
    fh = dnaio.open(fq1, fq2, fileformat="fastq", mode="r")
    if paired_out:
        outfh = dnaio.open(fq_out[0], fq_out[1], fileformat="fastq", mode="w")
    else:
        outfh = dnaio.open(fq_out[0], fileformat="fastq", mode="w")

    if use_multi:
        if paired_out:
            outfh_multi = dnaio.open(fqout_multi[0], fqout_multi[1], fileformat="fastq", mode="w")
        else:
            outfh_multi = dnaio.open(fqout_multi[0], fileformat="fastq", mode="w")
    
    for r1, r2 in fh:
        stat_Dict["total"] += 1
        
        start_pos = 0
        end_pos = 0
        sequence = r1.sequence
        qualities = r1.qualities
        
        # deal with shift
        if shift:
            shift_pos = sequence[:7].find(shift_pattern)
            if shift_pos < 0:
                stat_Dict["no_anchor"] += 1
                logger.debug(f"{r1.name},{sequence},{sequence[:7]} no anchor!")
                continue
            else:
                start_pos = shift_pos + 1
        
        # get barcode/umi/quality sequence          
        old_seqs = defaultdict(list)
        new_seqs = defaultdict(list)
        seq_quals = defaultdict(list)
        B = 0
        L = 0
        is_valid = True
        is_multi = False
        is_correct = False
        is_B_no_correction = False
        is_L_no_correction = False

        for _, (code, n) in enumerate(r1_structure):
            end_pos = start_pos + n
            seq = sequence[start_pos:end_pos]
            quals = qualities[start_pos:end_pos]

            if code == "B":
                old_seqs["B"].append(seq)
                seq_quals["B"].append(quals)

                if barcode_list_flag: # match barcode in whitelist
                    if seq in barcode_wl_dict.get(B, barcode_wl_dict[0]):
                        new_seqs["B"].append({seq})
                    else:
                        if do_B_correction:
                            bc_set = get_new_bc(seq, barcode_wl_dict.get(B, barcode_wl_dict[0]), match_type_dict.get(B, match_type_dict[0]))

                            if len(bc_set) == 0:
                                is_valid = False
                                is_B_no_correction = True
                                logger.debug(f"{r1.name},{sequence[:start_pos]}[{seq}]{sequence[end_pos:]},{seq} no barcode!")
                                break
                            elif len(bc_set) == 1:
                                new_seqs["B"].append(bc_set)
                                is_correct = True
                                logger.debug(f"{r1.name},{sequence[:start_pos]}[{seq}]{sequence[end_pos:]},{seq} -> {list(bc_set)} do_B_correction!")
                            else:
                                new_seqs["B"].append(bc_set)
                                logger.debug(f"{r1.name},{sequence[:start_pos]}[{seq}]{sequence[end_pos:]},{seq} -> {list(bc_set)} do_B_correction!")
                                is_multi = True
                        else:
                            is_valid = False
                            break
                else:
                    new_seqs["B"].append({seq})
                B += 1

            elif code == "L":   
                if linker_list_flag: # linker correction step
                    if seq in linker_wl_dict.get(L, linker_wl_dict[0]):
                        pass
                    else:
                        if do_L_correction:
                            lk_set = get_new_bc(seq, linker_wl_dict.get(L, linker_wl_dict[0]))
                            if len(lk_set) == 0:
                                is_valid = False
                                is_L_no_correction = True
                                logger.debug(f"{r1.name},{sequence[:start_pos]}[{seq}]{sequence[end_pos:]},{seq} -> {list(lk_set)} no linker!")
                                break
                        else:
                            is_valid = False
                            break
                L += 1
                
            elif code == "U":
                new_seqs["U"].append(seq)
                seq_quals["U"].append(quals)
                
            start_pos = start_pos + n

        # check double instances
        if is_valid:

            barcode_old = "".join(old_seqs["B"])
            Barcode_Counter[barcode_old] += 1

            #get base summary for umi/r2
            umi = "".join(new_seqs["U"])
            umi_q = "".join(seq_quals["U"])
            barcode_q = "".join(seq_quals["B"])
            
            UMI_GC_Counter, UMI_Q_Counter = summary(umi, umi_q, UMI_GC_Counter, UMI_Q_Counter)
            R2_GC_Counter, R2_Q_Counter = summary(r2.sequence, r2.qualities, R2_GC_Counter, R2_Q_Counter)
            
            r1.sequence = sequence[start_pos:]
            r1.qualities = qualities[start_pos:]
                        
            if is_multi: #write r2 multi files
                if use_multi:         
                    #update barcode quality
                    Barcode_Q_Counter.update(enumerate(barcode_q))
                    bc_new_lst = []
                    for element in itertools.product(*new_seqs["B"]):          
                        barcode_new = "".join(element)
                        bc_new_lst.append(barcode_new)
                        
                    bc_new_all = ":".join(bc_new_lst)
                    r2.name = "_".join([barcode_old, bc_new_all, umi, r2.name])
                    r1.name = "_".join([barcode_old, bc_new_all, umi, r1.name])
                    r1.sequence = sequence[start_pos:]
                    r1.qualities = qualities[start_pos:]
                    outfh_multi.write(r1, r2)
            else:  #write r2 files
                stat_Dict["valid"] += 1
                flag, r1, r2 = adapter_filter.filter(r1, r2)
                if flag:
                    if (not use_short_read) or len(r1) == 0 or len(r2) == 0:
                        if len(r1) < ATAC_R1_MINLEN or len(r2) < ATAC_R2_MINLEN:
                            stat_Dict["too_short"] += 1
                            continue
                    else:
                        stat_Dict["trimmed"] += 1

                barcode_new = "".join([_.pop() for _ in new_seqs["B"]])
                Barcode_GC_Counter, Barcode_Q_Counter = summary(barcode_old, barcode_q, Barcode_GC_Counter, Barcode_Q_Counter)
                
                #find alterations
                if is_correct:
                    _alt = "".join([str(i)+o for i, (o,n) in enumerate(zip(barcode_old, barcode_new)) if o != n])
                else:
                    _alt = "M"

                r2.name = "_".join([barcode_new, umi, _alt, r2.name])
                r1.name = "_".join([barcode_new, umi, _alt, r1.name])
                outfh.write(r1, r2)
            if is_correct:
                stat_Dict["B_corrected"] += 1
        else:
            if is_B_no_correction:
                stat_Dict["B_no_correction"] += 1

            if is_L_no_correction:
                stat_Dict["L_no_correction"] += 1
    if use_multi:
        outfh_multi.close()
    outfh.close()

    return {
            "stat": Counter(stat_Dict),
            "barcode_count": Barcode_Counter,
            "barcode_gc": Barcode_GC_Counter,
            "umi_gc": UMI_GC_Counter,
            "r2_gc": R2_GC_Counter,
            "barcode_q": Barcode_Q_Counter,
            "umi_q": UMI_Q_Counter,
            "r2_q": R2_Q_Counter
        }

def barcode_main(fq1:list, fq2:list, samplename: str, outdir:str,
                 barcode:list=[], match_type:list=[], shift:str=True, shift_pattern:str="A",
                 structure:str="B8L8B8L10B8U12T15", linker: list=[],
                 core:int=4, do_B_correction=True, do_L_correction=True,
                 use_multi=True, use_short_read=False, adapter1=[["TTTTTTTTTTTT", "5"], ],
                 adapter2=[["AAAAAAAAAAAA", "3"], ], paired_out=True, **kwargs):
    logger.info("extract barcode start!")
    #parse r1 structure
    r1_structure = parse_structure(structure)
    
    #get wl dict for bc/linker
    # 当barcode或linker为空时，返回空字典
    barcode_wl_dict = read_file(barcode)
    linker_wl_dict = read_file(linker)
    match_type_dict = {ind: val for ind, val in enumerate(match_type)}

    if len(barcode_wl_dict)>0 and do_B_correction:
        logger.info("barcode one base mismatch allowed.")
    else:
        logger.info("barcode mismatch NOT allowed.")

    if "L" in structure:
        if len(linker_wl_dict)>0 and do_L_correction:
            logger.info("linker one base mismatch allowed.")
        else:
            logger.info("linker mismatch NOT allowed.")

    if use_multi:
        logger.info("rescue barcode match multi barcode in whitelist.")
    else:
        logger.info("ignore barcode match multi barcode in whitelist.")
    
    #worker function
    worker_func = partial(
        process_barcode,
        r1_structure=r1_structure,
        shift=shift,
        shift_pattern=shift_pattern,
        barcode_wl_dict=barcode_wl_dict,
        linker_wl_dict=linker_wl_dict,
        match_type_dict=match_type_dict,
        do_B_correction=do_B_correction,
        do_L_correction=do_L_correction,
        use_multi=use_multi,
        use_short_read=use_short_read,
        adapter1=adapter1,
        adapter2=adapter2,
    )
    
    stat = QcStat()
    
    atacname = f"{samplename}_A"
    os.makedirs(f"{outdir}/step1", exist_ok=True)
    fqout = os.path.join(outdir, f"step1/{atacname}")
    fqout_multi = os.path.join(outdir, f"step1/{atacname}_multi")
    json_multi = os.path.join(outdir, f"step1/{atacname}_multi.json")
    
    pipeline = Pipeline(
        func=worker_func,
        fq1=fq1,
        fq2=fq2,
        fqout=fqout,
        fqout_multi=fqout_multi,
        stat=stat,
        core=core,
        paired_out=paired_out
    )
    pipeline.run()

    fqout1 = f"{fqout}_1.fq.gz"
    fqout2 = f"{fqout}_2.fq.gz"
    if use_multi:
        # find the multiple barcodes
        logger.info("deal multi start!")
        fqout_multi1 = f"{fqout_multi}_1.fq.gz"
        fqout_multi2 = f"{fqout_multi}_2.fq.gz"
        adapter_filter = AdapterFilter(adapter1=adapter1, adapter2=adapter2)
        multi_stat = defaultdict(int)
        with dnaio.open(fqout1, fqout2, mode="a") as f:
            fh = dnaio.open(fqout_multi1, fqout_multi2, fileformat="fastq", mode="r")
            for r1, r2 in fh:
                multi_stat["total"] += 1
                final_barcode = None
                
                bc_old, r2_candidate, umi, r2_name = r2.name.split("_", 3)
                r2_candidate = r2_candidate.split(":")
                
                read_num = 0
                for _ in sorted(r2_candidate):
                    v = stat.data["barcode_count"].get(_, 0)
                    if v > read_num:
                        read_num = v
                        final_barcode = _

                if not final_barcode:
                    multi_stat["B_no_correction"] += 1
                    stat.data["stat"]["B_no_correction"] += 1
                    continue
                    
                multi_stat["valid"] += 1
                stat.data["stat"]["valid"] += 1
                # stat.data["barcode_count"][_] += 1

                flag, r1, r2 = adapter_filter.filter(r1, r2)
                if flag:
                    if (not use_short_read) or len(r1) == 0 or len(r2) == 0:
                        if len(r1) < ATAC_R1_MINLEN or len(r2) < ATAC_R2_MINLEN:
                            multi_stat["too_short"] += 1
                            stat.data["stat"]["too_short"] += 1
                            continue
                    else:
                        multi_stat["trimmed"] += 1
                        stat.data["stat"]["trimmed"] += 1

                alt_l = [str(i)+o for i, (o,n) in enumerate(zip(bc_old, final_barcode)) if o != n]
                _alt = "".join([alt for alt in alt_l])
                r2.name = "_".join([final_barcode, umi, _alt, r2_name])

                r1.name = r2.name
                f.write(r1, r2)

        with open(json_multi, "w") as fh:
            json.dump(multi_stat, fp=fh, indent=4)

    del stat.data["barcode_count"]
    logger.info("deal multi done!")
    stat.data["stat"]["chemistry"] = kwargs.get("chemistry", "custom")
    stat.data["stat"]["atacname"] = atacname
    stat.save(os.path.join(outdir, f"{atacname}_summary.json"))
    logger.info("extract barcode done!")
    return fqout1, fqout2
