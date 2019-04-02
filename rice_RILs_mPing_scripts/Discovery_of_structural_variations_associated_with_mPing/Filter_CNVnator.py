#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python Filter_CNVnator.py --input RIL275_correction

Read all bed files in RIL275_correction, filter bed files with normalized depth, high/low coverage region in sequencing, calls present in parent and RILs.
The results will put in RIL275_correction_filtered
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#1       1581901 1582600 0.0176711       1       1581804 1582688 699
def overlap2bed(infile, outfile, rate):
    data = defaultdict(lambda : int())
    ofile = open(outfile, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                index = '%s_%s_%s_%s' %(unit[0], unit[1], unit[2], unit[3])
                data[index] += int(unit[-1])
    for idx in sorted(data.keys()):
        chrs, start, end, depth = re.split(r'_', idx)
        if float(data[idx])/(int(end)-int(start)) < rate:
            print >> ofile, '%s\t%s\t%s\t%s' %(chrs, start, end, depth)
    ofile.close()
    return data

#filtered by normalized depth
def filter_depth(infile, outfile):
    r = re.compile(r'Sy|Un')
    ofile = open(outfile, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                #if not line.startswith(r'Chr'):
                #    unit[0] = 'Chr%s' %(unit[0])
                unit[0] = re.sub(r'Chr', r'', unit[0])
                if r.search(unit[0]):
                    continue
                #if float(unit[3]) <= 0.05 or float(unit[3]) >= 2:
                if float(unit[3]) <= 0.05:
                    print >> ofile, '\t'.join(unit)
    ofile.close()

def filter_high_low(ril, filter_depth_bed, MSU_black, filter_high_low_bed, high_low_NB_bed, high_low_HEG4_bed):
    cmd = []
    fraction = 0.50
    cmd.append('bedtools intersect -wao -a %s -b %s > %s.2.MSU7_black.overlap' %(filter_depth_bed, MSU_black, ril))
    os.system('bedtools intersect -wao -a %s -b %s > %s.2.MSU7_black.overlap' %(filter_depth_bed, MSU_black, ril))  
    overlap2bed('%s.2.MSU7_black.overlap' %(ril), '%s.2.MSU7_black.bed' %(ril), fraction)
    cmd.append('bedtools intersect -wao -a %s.2.MSU7_black.bed -b %s > %s.3.MSU7_NB_black.overlap' %(ril, high_low_NB_bed, ril))
    os.system('bedtools intersect -wao -a %s.2.MSU7_black.bed -b %s > %s.3.MSU7_NB_black.overlap' %(ril, high_low_NB_bed, ril))
    overlap2bed('%s.3.MSU7_NB_black.overlap' %(ril), '%s.3.MSU7_NB_black.bed' %(ril), fraction)
    cmd.append('bedtools intersect -wao -a %s.3.MSU7_NB_black.bed -b %s > %s' %(ril, high_low_HEG4_bed, re.sub(r'.bed', r'.overlap',filter_high_low_bed)))
    os.system('bedtools intersect -wao -a %s.3.MSU7_NB_black.bed -b %s > %s' %(ril, high_low_HEG4_bed, re.sub(r'.bed', r'.overlap',filter_high_low_bed)))
    overlap2bed(re.sub(r'.bed', r'.overlap',filter_high_low_bed), filter_high_low_bed, fraction)
    for c in cmd:
        print c
        #os.system(c)

def filter_parent(ril, filter_high_low_bed, filter_parent_bed, parent_NB_bed, parent_HEG4_bed):
    cmd = []
    fraction = 0.50
    cmd.append('bedtools intersect -wao -a %s -b %s > %s.5.filter_parent_NB.overlap' %(filter_high_low_bed, parent_NB_bed, ril))
    os.system('bedtools intersect -wao -a %s -b %s > %s.5.filter_parent_NB.overlap' %(filter_high_low_bed, parent_NB_bed, ril))
    overlap2bed('%s.5.filter_parent_NB.overlap' %(ril), '%s.5.filter_parent_NB.bed' %(ril), fraction)
    cmd.append('bedtools intersect -wao -a %s.5.filter_parent_NB.bed -b %s > %s' %(ril, parent_HEG4_bed, re.sub(r'.bed', r'.overlap',filter_parent_bed)))
    os.system('bedtools intersect -wao -a %s.5.filter_parent_NB.bed -b %s > %s' %(ril, parent_HEG4_bed, re.sub(r'.bed', r'.overlap',filter_parent_bed)))
    overlap2bed(re.sub(r'.bed', r'.overlap',filter_parent_bed), filter_parent_bed, fraction)
    for c in cmd:
        print c
        #os.system(c)

def get_line_number(infile):
    num_lines = sum(1 for line in open(infile))
    return num_lines

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    directory='/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/cnvnator/input'
    MSU_black       = '%s/%s' %(directory, 'MSU7.black.bed')
    high_low_NB_bed = '%s/%s' %(directory, 'NB_P.high_low_500bp.bedgraph')
    high_low_HEG4_bed = '%s/%s' %(directory, 'HEG4_P.high_low_500bp.bedgraph')
    parent_NB_bed   = '%s/%s' %(directory, 'NB_P.bam.readdepth.bed')
    parent_HEG4_bed = '%s/%s' %(directory, 'HEG4_P.bam.readdepth.bed')
    beds = glob.glob('%s/*.bed' %(args.input))
    os.system('mkdir %s_filtered' %(args.input))
 
    data = defaultdict(lambda : list())
    sumfile = open('%s_filtered.sum' %(args.input), 'w')
    print >> sumfile, 'RILs\tRaw\tFiltered_depth\tFiltered_black\tFiltered_parents'
    for bed in sorted(beds):
        ril = re.split(r'\.', os.path.split(bed)[1])[0]
        print ril
        #normalized depth
        filter_depth_bed = '%s.1.readdepth_filter_depth.bed' %(ril)
        filter_depth(bed, filter_depth_bed)
        #high and low coverage region in illumina sequence, chromosome start, end and centromere
        filter_high_low_bed = '%s.4.readdepth_filter_black.bed' %(ril)
        filter_high_low(ril, filter_depth_bed, MSU_black, filter_high_low_bed, high_low_NB_bed, high_low_HEG4_bed)
        #calls present in parents strain
        filter_parent_bed  = '%s.6.readdepth_filter_parent.bed' %(ril)
        filter_parent(ril, filter_high_low_bed, filter_parent_bed, parent_NB_bed, parent_HEG4_bed)
        data[ril] = [get_line_number(bed), get_line_number(filter_depth_bed), get_line_number(filter_high_low_bed), get_line_number(filter_parent_bed)]
        print >> sumfile, '%s\t%s' %(ril, '\t'.join(map(str, data[ril])))
        os.system('mv %s.* %s_filtered' %(ril, args.input))
    sumfile.close()

if __name__ == '__main__':
    main()

