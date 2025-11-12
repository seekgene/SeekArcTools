import os

__srcdir = os.path.dirname(os.path.abspath(__file__))

ADAPTERS = {
    "MM": [["AAAAAAAAAAAA", "3"],  ],
    "DD": [["AAGCAGTGGTATCAACGCAGAGTACATGGG", "5"], ],
}

R1_MINLEN = 20
R2_MINLEN = 50
ATAC_R1_MINLEN = 20
ATAC_R2_MINLEN = 50

CHEMISTRY = {
    "DD_AG":{
        'shift': False,
        'structure': 'B17U12X7',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
        'adapter1': [["TTGCTGT", "5"], ["AAGCAGTGGTATCAACGCAGAGTACATGG", "3"], ], ## 7F TSO
        'adapter2': [["CCATGTACTCTGCGTTGATACCACTGCTT", "5"], ["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", "3"], ], ## TSO-rev SP1-rev
        'sc5p': None,
        'match_type': (1,),
    },
    "DD_AA":{
        'shift': False,
        'structure': 'B17U12',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
        'adapter1': [["CTGTCTCTTATACACATCTCCGAGCCCACGAGAC", "3"], ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"]], ## SP2 SP1 17L19ME, ["CGTCCGTCGTTGCTCGTAGATGTGTATAAGAGACAG", "5"], ["AGATGTGTATAAGAGACAG", "5"]
        'adapter2': [["GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG", "5"], ["GGTGATTAACGATGCATTAGATGTGTATAAGAGACAG","5"], ["CTGTCTCTTATACACATCTACGAGCAACGACGGACG", "3"]], ## SP2-rev zhixiaoling 17L19ME-rev
        'match_type': (1,),
    },
    "DD5_AG":{
       'shift': False,
        'structure': 'B17U12X13',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
        'match_type': (1,),
        'sc5p': True,
        'adapter1': [["TTTATTATATGGG", "5"], ["AAGCAGTGGTATCAACGCAGAGTACATGG", "3"], ], 
        'adapter2': [["CCATGTACTCTGCGTTGATACCACTGCTT", "5"], ["CCCATATAATAAA", "3"], ],
    },
}