import os
import json
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
