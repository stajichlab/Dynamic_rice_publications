#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import subprocess
sys.path.append('%s/lib' %(os.getcwd()))
from excision_PingA import bamcheck, bamcheck_simple_cover, convert_MAP, genotyping, convert_MAP_SNP, genotyping_SNP
import glob

def usage():
    test="name"
    message='''
python mPing_Boundary_Coverage_PingA.py --bam_ref 3k_stowaway_strains_Ref_bam --bam_pseudo 3k_stowaway_strains_Pseudo_bam

Check read coverage at mPing boundary from bam files of read2reference and read2mping_flanking. We will create matrix
of RILs and mPing, which gives information of each end of mPing.

    '''
    print message

#Chr1    PseudoGenome    Transposable_element    1132975 1133405 -       .       .       ID=Chr1_1132975_1133405;Original_ID=Chr1.1132977.spanners;TE=mping;TSD=TAA;
#Chr1    PseudoGenome    Transposable_element    2642232 2647573 +       .       .       ID=Chr1_2642232_2647573;Original_ID=Chr1.2640500;TE=ping;TSD=TAA;
def id_mapping(infile, mping2ID_0):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                start = int(unit[3]) 
                end   = int(unit[4])
                chro  = unit[0]
                strand= unit[6]
                temp  = defaultdict(str)
                attrs = re.split(r';', unit[8])
                for attr in attrs:
                    if not attr == '':
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value
                temp['Original_ID'] = re.sub(r'.spanners', r'', temp['Original_ID'])
                temp_ids   = re.split(r'_', temp['ID'])
                temp['ID'] = '%s:%s-%s' %(temp_ids[0], temp_ids[1], temp_ids[2])
                temp_oids  = re.split(r'\.', temp['Original_ID'])
                #print temp_oids
                temp['Original_ID'] = '%s:%s-%s' %(temp_oids[0], str(int(temp_oids[1])-2), temp_oids[1])
                mping2ID_0[temp['ID']]  = temp['Original_ID']
                #print '%s\t%s' %(temp['ID'], mping2ID_0[temp['ID']])
    return data

def decode(flag):
    if int(flag) == 0 or int(flag) == 4:
        return 'covered'
    elif int(flag) == 1:
        return 'clipped'
    elif int(flag) == 5:
        return 'covered/clipped'
    else:
        return 'unknown'

def decode_gt(gt):
    if gt == 'NA':
        return 'NA'
    elif str(gt) == '0':
        return 'NB'
    elif str(gt) == '1':
        return 'HEG4'
    else:
        return 'NA'

def sort_mping_chr(mpings):
    data = defaultdict(lambda : defaultdict(lambda : str()))
    for mping in mpings.values():
        reg     = re.compile(r'Chr(\d+):(\d+)\-(\d+)')
        match   = reg.search(mping)
        chro    = match.groups(0)[0]
        start   = match.groups(0)[1]
        end     = match.groups(0)[2]
        data[chro][start] = mping
    sorted_mping = []
    for c in sorted(data.keys(), key=int):
        for s in sorted(data[c].keys(), key=int):
            sorted_mping.append(data[c][s])
    return sorted_mping

def get_rils(bams):
    rils = []
    for bam in bams:
        ril = os.path.split(bam)[1]
        ril = re.sub(r'_Pseudo.NM2.bam', r'', ril)
        #ril = re.sub(r'RIL', r'', ril)
        rils.append(ril)
    return rils


def bamcheck_ref(bam, mping, bamck_file, ril):
    reg     = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match   = reg.search(mping)
    chro    = match.groups(0)[0]
    l_start = int(match.groups(0)[1]) - 2
    l_end   = int(match.groups(0)[1]) + 2
    l_mping = '%s:%s-%s' %(chro, l_start, l_end)
    r_start = int(match.groups(0)[2]) - 2
    r_end   = int(match.groups(0)[2]) + 2
    r_mping = '%s:%s-%s' %(chro, r_start, r_end)
    l_flag  = bamcheck_simple_cover(bam, l_mping, bamck_file)
    r_flag  = bamcheck_simple_cover(bam, r_mping, bamck_file)
    return (l_flag, r_flag)

#B001    ERS470219
def acc2name(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[1]] = unit[0]
    return data



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_ref')
    parser.add_argument('--bam_pseudo')
    parser.add_argument('--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.bam_pseudo) > 0 and len(args.bam_ref) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = 'mPing_boundary'
    if not args.bam_ref:
        args.bam_ref = '../input/RILs_ALL_bam'
    if not args.bam_pseudo:
        args.bam_pseudo = '../input/RILs_ALL_unmapped_mping_bam'

    #we use mping gff from pseudogenome as key to project to everything 
    #mping2ID_0  = defaultdict(lambda : str()) #ID_0 is the mping id from original call in HEG4, Chr1.1132977
    #id_mapping(args.gff_pseudo, mping2ID_0)
    
    #bin map and snp genotype
    #binmap = convert_MAP(args.bin_map)
    #snpmap = convert_MAP_SNP(args.bin_map)
 
    names = acc2name('rice_name_acc.list')

    #go through RILs, for each ril determine the status of each mPing
    bamcheck_file_pseudo = '%s.bamcheck_pseudo.txt' %(args.project)
    bamcheck_file_ref    = '%s.bamcheck_ref.txt' %(args.project)
    mping_status      = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str())))
    mping_status_ref  = defaultdict(lambda : defaultdict(lambda : str()))
    #mping_bin_gt  = defaultdict(lambda : defaultdict(lambda : str()))
    #mping_snp_gt  = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str())))
    #rils = [1, 2, 3, 4]
    bams = glob.glob('%s/*.NM2.bam' %(args.bam_pseudo))
    rils = get_rils(bams)
    for ril in sorted(rils):
        bam_ref    = '%s/%s_Ref.NM2.bam' %(args.bam_ref, ril)
        bam_pseudo = '%s/%s_Pseudo.NM2.bam' %(args.bam_pseudo, ril)
        #print 'ril: %s, %s' %(ril, bam_pseudo)
        ping_pseudo    = 'Chr1:1001-6341'
        l_flag, r_flag = bamcheck_ref(bam_pseudo, ping_pseudo, bamcheck_file_pseudo, ril)
        ping_ref       = 'Chr1:1305-1307'
        ref_flag       = bamcheck(bam_ref, ping_ref, bamcheck_file_ref, ril) 
        print '{}\t{}\t{}\t{}\t{}\t{}'.format(ril, names[ril], decode(ref_flag), ref_flag, decode(l_flag), decode(r_flag)) 
        #print '{}\t{}\t{}\t{}'.format(ril, decode(ref_flag), decode(ref_flag), decode(ref_flag))

if __name__ == '__main__':
    main()

