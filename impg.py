#!/usr/bin/python
# (c) 2017-2022 Huwenbo Shi

import numpy as np
import pandas as pd
import argparse, math, sys, logging, time
from pysnptools.snpreader import Bed
from src.estimation import *

title_str = """\
@----------------------------------------------------------@
       |        IMPG       |      v2.0      |    30/October/2017  |
       |----------------------------------------------------------|
       |  (C) 2017 Huwenbo Shi, GNU General Public License, v3    |
       |----------------------------------------------------------|
       |  For documentation, citation & bug-report instructions:  |
       |   http://bogdan.bioinformatics.ucla.edu/software/hess/   |
       @----------------------------------------------------------@\
"""

# main function
def main():

    # get command line argument and initialize log
    args = get_command_line()
    init_log(args)
    argmap = check_command_line(args)
    
    print argmap

    # end the log
    end_log()

# initialize log
def init_log(args):

    # get log file name
    log_file_name = args.out
    if args.chrom != None:
        log_file_name = log_file_name + '_chr' + args.chrom
    log_file_name += '.log'

    # create the log file
    log_format = '[%(levelname)s] %(message)s'
    logging.basicConfig(filename=log_file_name, filemode="w",
        level=logging.DEBUG, format=log_format)

    # add stderr as a stream handler
    stderr_handler = logging.StreamHandler(sys.stderr)
    formatter = logging.Formatter(log_format)
    stderr_handler.setFormatter(formatter)
    logging.getLogger().addHandler(stderr_handler)

    # log time and command issued
    logging.info(title_str)
    specified = set([val for val in sys.argv if val[0] == '-'])
    cmd_str = sys.argv[0] + ' \\\n'
    for arg in vars(args):
        if ('--' + arg) in specified:
            param = getattr(args, arg)
            if type(param) == list:
                param = ' '.join([str(p) for p in param])
            elif type(param) == bool:
                param = ''
            cmd_str += '        --{} {} \\\n'.format(arg, param)
    cmd_str = cmd_str.strip()[0:len(cmd_str)-3]
    cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    logging.info('Command started at: %s' % cur_time)
    logging.info('Command issued:\n    {}'.format(cmd_str))

# end the log
def end_log():
    cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    logging.info('Command finished at: %s' % cur_time)

# get command line input
def get_command_line():
    
    parser = argparse.ArgumentParser(description='Impute GWAS summary '
        'association data under the multivariate normal distribution')

    parser.add_argument('--bfile', dest='bfile', type=str, required=True,
        default=None, help='Reference panel file in PLINK file format')

    parser.add_argument('--chrom', dest='chrom', type=str, required=True,
        help='Specifies the chromosome number')
   
    parser.add_argument('--sumstats', dest='sumstats', type=str,
        required=True, help='Specifies the summary statistics file')
 
    parser.add_argument('--window-size', dest='window-size', type=int,
        help='Size of the window (default 1,000,000)', default=1000000,
        required=False)
   
    parser.add_argument('--buffer-size', dest='buffer-size', type=int,
        help='Size of the buffer (default 250,000)', default=250000,
        required=False)

    parser.add_argument('--lambda', dest='lambda', type=float,
        help='Regularization constant', default=0.1, required=False)

    parser.add_argument('--out', dest='out', type=str,
        help='Output file name', required=True)
    
    return parser.parse_args()

# check command line
def check_command_line(args):

    argmap = dict()
    for arg in vars(args):
        argmap[arg] = getattr(args, arg)

    return argmap

if(__name__ == '__main__'):
    main()
