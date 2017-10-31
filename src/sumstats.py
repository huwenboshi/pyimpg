import numpy as np
import pandas as pd
import os, sys, gzip, logging

# required columns for the summary stats file
required = ['SNP', 'CHR', 'BP', 'A1', 'A2', 'Z']

# magic bits to discern file type
magic_dict = {
    "\x1f\x8b\x08": "gz",
    "\x42\x5a\x68": "bz2",
    "\x50\x4b\x03\x04": "zip"
    }
max_len = max(len(x) for x in magic_dict)

# define equivalent alleles
equiv = dict()
equiv["AC"] = set(["TG", "AC", "TC", "AG"])
equiv["AG"] = set(["TC", "AG", "TG", "AC"])
equiv["CA"] = set(["GT", "CA", "GA", "CT"])
equiv["CT"] = set(["GA", "CT", "GT", "CA"])
equiv["TC"] = set(["AG", "TC", "AC", "TG"])
equiv["TG"] = set(["AC", "TG", "AG", "TC"])
equiv["GA"] = set(["CT", "GA", "CA", "GT"])
equiv["GT"] = set(["CA", "GT", "CT", "GA"])

# define reversed alleles
reverse = dict()
reverse["AC"] = set(["GT", "CA", "CT", "GA"])
reverse["AG"] = set(["CT", "GA", "GT", "CA"])
reverse["CA"] = set(["TG", "AC", "AG", "TC"])
reverse["CT"] = set(["AG", "TC", "TG", "AC"])
reverse["TC"] = set(["GA", "CT", "CA", "GT"])
reverse["TG"] = set(["CA", "GT", "GA", "CT"])
reverse["GA"] = set(["TC", "AG", "AC", "TG"])
reverse["GT"] = set(["AC", "TG", "TC", "AG"])

# define strand ambiguous alleles
ambiguous = set(["AT", "CG", "TA", "GC"])

def file_type(filename):
    """
    Check the type of the file

    Args:
        filename: The file for which the type will be checked

    Returns:
        The type of the file.
    """

    with open(filename) as f:
        file_start = f.read(max_len)
    
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            if filetype != 'gz':
                logging.error('{} compression type is not supported'\
                    .format(filetype))
                sys.exit(1)
            return filetype
    
    return "txt"

class SumStats(object):
    """
    Handles loading, filtering, and extraction of summary association data
    """

    def __init__(self, filename, chrom):
        """
        Load GWAS summary association data for the specified chromosome.
        Perform initial filtering, including removing SNPs without rs ID and
        SNPs with allele length greater than 1.

        Args:
            filename: File name of the summary association data.
            chrom: The chromosome for which the data will bed loaded.

        Returns:
            A data frame containing the GWAS summary association data.

        Raises:
            KeyError: Raises an exception.
        """
    
        # check the file type and create a file handle
        if not os.path.exists(filename):
            logging.error('{} does not exist.'.format(filename))
            sys.exit(1)
        ftype = file_type(filename)
        if ftype == 'gz':
            fhandle = gzip.open(filename, 'rb')
        else:
            fhandle = open(filename, 'r')

        # get index of the columns
        cols = fhandle.readline().strip().split()
        for req in required:
            if req not in cols:
                logging.error('Missing column {}'.format(req))
                sys.exit(1)
        idx_map = dict(zip(cols, [required.index(name) for name in required])) 

        # parse the summary stats file
        sumstats = dict(zip(required, [[] for name in required]))
        for line in fhandle:
            cols = line.strip().split()
            tmp = dict()
            for name in required:
                idx = idx_map[name]; val = cols[idx]
                if name == 'CHR' and val != chrom: break
                if name == 'SNP' and val[0:2] != 'rs': break
                if (name == 'A1' or name == 'A2') and len(val) != 1: break
                if name == 'BP' or name == 'Z': val = float(val)
                tmp[name] = val
            if len(tmp) == len(required):
                for key in tmp:
                    sumstats[key].append(tmp[key])

        # close the file handle
        fhandle.close()

        # create the data frame
        sumstats = pd.DataFrame(sumstats)

        # log the results
        logging.info('Loaded {} SNPs with rs IDs and single-letter alleles'
            ' on chromosome {} from the GWAS summary data file'\
            .format(sumstats.shape[0], chrom))

        self.sumstats = sumstats

    def filter_sumstats(self, refpanel_snpmap):
        """
        Filter out SNPs with ambiguous rs ID, BP, alleles. Flip signs
        accordingly based on the alleles.
        """
        
        # sort the data based on position
        self.sumstats = self.sumstats.sort_values(by='BP')       
        
        # remove duplicates based on both rs ID and position
        self.sumstats = self.sumstats.drop_duplicates('SNP', keep=False)
        self.sumstats = self.sumstats.drop_duplicates('BP', keep=False)

        # merge summstats and reference panel
        self.sumstats = self.sumstats.merge(refpanel_snpmap, on='SNP')
        self.sumstats = self.sumstats[['SNP', 'CHR', 'BP_y', 'A1_x',
            'A2', 'Z']]
        self.sumstats.columns = ['SNP', 'CHR', 'BP', 'A1', 'A2', 'Z']

        # re-index snps
        self.sumstats = self.sumstats.reset_index(drop=True)

        # remove ambiguous snps and flip sign accordingly (A1 is the set bit
        # in the reference panel)
        refpanel_snpmap['A1A0'] = refpanel_snpmap['A1'] + refpanel_snpmap['A0']
        refpanel_a1a0 = dict(zip(refpanel_snpmap.SNP, refpanel_snpmap.A1A0))
        flip = []; filt = []
        for i in xrange(self.sumstats.shape[0]):
            snp = self.sumstats['SNP'][i]
            a1a2 = self.sumstats['A1'][i] + self.sumstats['A2'][i]
            a1a0 = refpanel_a1a0[snp]
            if a1a0 in ambiguous or a1a2 in ambiguous: filt.append(i)
            elif a1a0 in reverse[a1a2]: flip.append(i)
            elif a1a0 in equiv[a1a2]: pass
            else: filt.append(i)
        
        self.sumstats.loc[flip, 'Z'] *= -1.0
        print filt
        filtered = self.sumstats.loc[filt,:]
        self.sumstats = self.sumstats.drop(filt)

        # re-index snps
        self.sumstats = self.sumstats.reset_index(drop=True)

        logging.info('{} SNPs left after filtering'\
            .format(self.sumstats.shape[0]))
        
        return filtered

    def get_locus(self, start, stop):
        """
        Extract the data in the region defined by start (inclusive) and
        stop (exclusive)
        """
        sumstats_locus = self.sumstats[(self.sumstats['BP'] >= start) &
                                    (self.sumstats['BP'] < stop)]
        sumstats_locus = sumstats_locus.reset_index(drop=True)

        return sumstats_locus

