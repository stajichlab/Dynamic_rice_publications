#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import random

def usage():
    test="name"
    message='''
python Random_sample_validation.py --input 25

--input: number of random loci to select

    '''
    print message


def runjob(script, lines, cpu, queue):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines %s --interval 120 --task %s --mem 15G --time 10:00:00 --queue %s --convert no %s' %(lines, cpu, queue, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr10:21716391-21716819	mping	Ref	1	RIL220
def readtable(infile, data):
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = unit[1:]
    return data


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

    data = defaultdict(lambda: list())
    nonref = 'output_MSU7_nonreference_514_MP3_mPing_GT_Ping_code.mping_excision.list'
    ref    = 'output_MSU7_reference_57_MP3_mPing_GT_Ping_code.mping_excision.list' 
    readtable(nonref, data)
    readtable(ref, data)

    #random loci
    loci_num_total     = 0
    loci_num_excised = 0
    loci_list_excised = []
    loci_list_excised_nonref = []
    loci_list_excised_ref = []
    loci_list_random    = [] 
    for mping in sorted(data.keys()):
        loci_num_total += 1
        if int(data[mping][2]) > 0:
            loci_num_excised += 1
            loci_list_excised.append(mping)
            if data[mping][1] == 'Nonref':
                loci_list_excised_nonref.append(mping)
            else:
                loci_list_excised_ref.append(mping)
            
    #loci_list_random = random.sample(loci_list_excisized, int(args.input))
    loci_list_random.extend(random.sample(loci_list_excised_nonref, int(args.input)*2/3))
    loci_list_random.extend(random.sample(loci_list_excised_ref, int(args.input)*1/3))
    print 'Total loci: %s' %(loci_num_total) 
    print 'Excised loci: %s' %(loci_num_excised)    


    #random RILs
    print 'Loci\tTest_RIL\tTE\tType\t#RILs\tRILs'
    for mping in sorted(loci_list_random):
        rils = re.split(',', data[mping][3])
        ril  = random.sample(rils, 1)[0]
        print '%s\t%s\t%s' %(mping, ril, '\t'.join(data[mping]))

if __name__ == '__main__':
    main()

