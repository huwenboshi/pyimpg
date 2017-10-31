import numpy as np
import pandas as pd
import sys, gzip, os, logging
import pysnptools
from pysnptools.snpreader import Bed

def load_partition(filename, chrom):
    """
    Load the partition in bed format
    """

    # load the partition file
    partition = pd.read_table(filename, delim_whitespace=True)
    partition = partition[partition['chr'] == ('chr'+str(chrom))]
    logging.info('Loaded {} partitions on chromosome {}'\
        .format(partition.shape[0], chrom))
    partition = partition.reset_index(drop=True)

    if partition.shape[0] == 0:
        logging.error('No partition found')
        sys.exit(1)

    # compute the average width of the window
    mean_width = int(np.mean(partition['stop'] - partition['start']))
    logging.info('Average window size is {}'.format(mean_width))

    # check if the window is large enough
    if(mean_width < 1000000):
        logging.warning('Window size seems a bit too small. This may lead'\
            ' to biased estimate.')

    # check if the window are continuous
    for i in xrange(partition.shape[0]-1):
        if partition['stop'][i] != partition['start'][i+1]:
            logging.warning('Partition of the genome is not continuous')
            break

    return partition

class PlinkReader(object):
    """
    Handles reference panel input
    """
    
    def __init__(self, filename):
        """
        Returns a RefPanel object
        """
        
        # check if all files exist
        if not os.path.exists('{}.bim'.format(filename)) or \
           not os.path.exists('{}.bed'.format(filename)):
            logging.error('Missing files for the reference panel {}'\
                .format(filename))
            sys.exit(1)
        
        # read in the data
        self.snpdata = Bed(filename, count_A1=False)
        self.snpmap = pd.read_table('{}.bim'.format(filename),
            delim_whitespace=True, usecols=[1, 3, 4, 5], header=None,
            names=['SNP', 'BP', 'A0', 'A1'])

        nsnp = self.snpmap.shape[0]
        logging.info('{} SNPs read from reference panel'.format(nsnp))

    def get_map(self):
        """
        Returns the legend / map of the reference panel
        """
        return self.snpmap

    def get_locus(self, start, stop, maf_thres):
        """
        Returns the legend and genotype matrix at a locus, specified by
        start (inclusive) and stop (exclusive)
        """

        # extract the legend from the locus
        snpmap_locus = self.snpmap[(self.snpmap['BP'] >= start) &
                                  (self.snpmap['BP'] < stop)]
        nsnp_locus, _ = snpmap_locus.shape
        
        # extract the genotype data from the locus
        start_idx = snpmap_locus.index.values[0]
        stop_idx = snpmap_locus.index.values[nsnp_locus-1]
        snpdata_locus = self.snpdata[:, start_idx:stop_idx+1].read().val.T
        _, nindv = snpdata_locus.shape

        # impute missing value with the average
        nanidx = np.where(np.isnan(snpdata_locus))
        mean_geno = np.nanmean(snpdata_locus, axis=1)
        snpdata_locus[nanidx] = mean_geno[nanidx[0]]

        # reset the index
        snpmap_locus = snpmap_locus.reset_index(drop=True)

        # get minor allele frequency
        maf = np.sum(snpdata_locus, axis=1) / (2.0 * nindv)
        maf[maf > 0.5] = 1.0 - maf[maf > 0.5]

        # filter based on maf_thres
        drop_idx = np.where(maf < maf_thres)[0]
        snpmap_locus = snpmap_locus.drop(drop_idx)
        snpdata_locus = np.delete(snpdata_locus, drop_idx, axis=0)

        # reset the index
        snpmap_locus = snpmap_locus.reset_index(drop=True)

        # remove duplicates
        if(np.sum(snpmap_locus['SNP'].duplicated()) > 0 or
           np.sum(snpmap_locus['BP'].duplicated()) > 0):
            logging.warning('Duplicate SNPs found in reference panel')

        snpmap_locus = snpmap_locus.drop_duplicates('SNP', keep=False)
        snpmap_locus = snpmap_locus.drop_duplicates('BP', keep=False)
        snpdata_locus = snpdata_locus[snpmap_locus.index.values,:]

        # reset the index
        snpmap_locus = snpmap_locus.reset_index(drop=True)
        
        return (snpmap_locus, snpdata_locus)
