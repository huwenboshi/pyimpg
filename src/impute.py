import tarfile
import numpy as np
import pandas as pd
import scipy, scipy.stats
from sumstats import *
from refpanel import *

eps = 10.0**-8

def get_ld(snpdata):
    """
    Compute the LD matrix
    """
    ld = np.corrcoef(snpdata)
    ld = np.nan_to_num(ld)
    return ld

def impute(refpanel_fnm, sumstats_fnm, window_size, buffer_size,
    chrom, min_maf, lam, ld_prune, out_fnm):

    # create the plink file reader
    refpanel = PlinkReader(refpanel_fnm)

    # load the summary stats and apply the filter
    sumstats = SumStats(sumstats_fnm, chrom)
    filtered = sumstats.filter_sumstats(refpanel.get_map())

    print filtered
