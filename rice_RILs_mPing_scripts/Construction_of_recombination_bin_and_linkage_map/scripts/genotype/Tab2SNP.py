#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO
import subprocess
import multiprocessing as mp
import gzip

def usage():
    test="name"
    message='''
python Tab2SNP.py --input RILs_ALL_bam_correct

Convert genotype.tab to Maq.p1.map.pileup.SNP. 

    '''
    print message

#/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/bin/RILs_ALL_275line_core/NB.RILs.dbSNP.SNPs.Markers
#SNP_id  Allele
#0100021547A     A
def read_parents(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'SNP_id'): 
                unit = re.split(r'\t',line)
                data[unit[0]] = unit[1]
    return data

##CHROM  POS     REF     RIL103_0_GAGTGG_FC1213L5
#Chr1    31071   A       A/A
##0100031071A     GN278   G
def convert_tab2SNP(infile, markers, outfile):
    #outfile = re.sub(r'.genotype.tab', r'.Maq.p1.map.pileup.SNP', infile)
    ofile = open (outfile, 'w')
    with gzip.open (infile, 'r') as filehd:
        headers = re.split(r'\t', filehd.readline())
        rils    = re.split(r'_', headers[-1])
        ril     = re.sub(r'RIL', r'GN', rils[0])
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'): 
                unit = re.split(r'\t',line)
                #pos  = int(unit[0][2:-1])
                chrs  = re.sub(r'Chr', r'', unit[0]) 
                pos   = '%02d%08d%s' %(int(chrs), int(unit[1]), unit[2])
                if unit[3][0] == unit[3][2] and not unit[3][0] == '.':
                    if markers.has_key(pos):
                        print >> ofile, '%s\t%s\t%s' %(pos, ril, unit[3][0])
    ofile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
 

    markers = read_parents('NB.RILs.dbSNP.SNPs.Markers')
    snp_files = glob.glob('%s/*.genotype.tab.gz' %(args.input))
    #convert_tab2SNP('RILs_ALL_bam_correct/GN87.genotype.tab', markers)

    for tab in sorted(snp_files):
        snp = re.sub(r'.genotype.tab.gz', r'.Maq.p1.map.pileup.SNP', tab)
        if not os.path.exists(snp):
            print '%s to %s' %(tab, snp)
            convert_tab2SNP(tab, markers, snp)

if __name__ == '__main__':
    main()

