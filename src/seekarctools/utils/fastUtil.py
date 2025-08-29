import re
from collections import defaultdict
from .helper import logger


def read_gtf(gtf):
    genetype_table = []
    mt_regex = re.compile("^(MT|mt|Mt)-")
    with open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            tmp = line.strip().split("\t")
            if tmp[2] == "gene":
                if "gene_id" in tmp[-1] and 'type' in  tmp[-1]:
                    gene_id=re.findall("gene_id \"([A-Za-z0-9_\.\-\:/\ ()]+)\";",tmp[-1])[0]
                    gene_names=re.findall("gene_name \"([A-Za-z0-9_\.\-\:/\ ()]+)\"",tmp[-1])
                    genetype=re.findall("gene_type \"([A-Za-z0-9_\.\-\:/\ ()+]+)\"",tmp[-1])
                    biotype=re.findall("gene_biotype \"([A-Za-z0-9_\.\-\:/\ ()+]+)\"",tmp[-1])
                    if len(genetype) !=0 :
                        gene_biotype = genetype[0]
                    elif len(biotype) !=0 :
                        gene_biotype = biotype[0]
                    else:
                        logger.info(gene_id,gene_names)
                        continue
                    if len(gene_names)==0:
                          gene_name = gene_id
                    else:
                          gene_name = gene_names[0]
                    
                    # gene_table.append([ gene_id, gene_name ])
                    if gene_biotype == 'lincRNA' or gene_biotype == 'antisense' or gene_biotype == "lncRNA":
                        biotype="lnc"
                    else:
                        biotype="coding"
                    genetype_table.append([gene_id, biotype])
    return  genetype_table
