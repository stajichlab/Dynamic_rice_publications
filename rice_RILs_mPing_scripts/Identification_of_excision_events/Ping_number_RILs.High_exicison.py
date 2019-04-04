#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir
import glob

def usage():
    test="name"
    message='''
python Ping_number_RILs.High_exicison.py --csv ../../Compare_RILs_SV/reads_map/bin/mPing_boundary_mPing

Read from csv file of bourdary check files. Add Excision code and Ping code to file.

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Pings   Ping_Code       RIL
#NA      NA      RIL158
#NA      NA      RIL242
#0       NA      RIL39
def read_ping(infile):
    data = defaultdict(lambda: list())
    r = re.compile(r'RIL(\d+)')
    with open (infile, 'r') as filehd:       
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Ping'): 
                unit = re.split(r'\t',line)
                m = r.search(unit[2])
                ril = m.groups(0)[0] if m else 'NA'
                #print ril, unit[0] 
                data[int(ril)] = [unit[0], unit[1]]
    return data

#Chr8:24622994-24622996,Genotype_Bin,Pseudo_mPing_status_up,Pseudo_mPing_status_down,Ref_mPing_status
#RIL1,HEG4,unknown,unknown,covered
def read_csv(infile, ping_code, outdir):
    prefix = os.path.split(infile)[1]
    ofile = open('%s/%s' %(outdir, prefix), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r',',line)
                if unit[0].startswith(r'Chr'):
                    unit.extend(['Excision_Code', 'Ping_Number','Ping_Code'])
                elif unit[0].startswith(r'RIL'):
                    ril = re.sub(r'RIL', r'', unit[0])
                    #print ril, ping_code[int(ril)]
                    ping_number = ping_code[int(ril)][0]
                    ping_cd   = ping_code[int(ril)][1]
                    #print ping_cd
                    flag  = 'Unknown' #final excision code, Unknown: sequence coverage not tell anything; Check: diff between pseudo and ref; Insertion: mPing; Excision: both comfirmed
                    flag1 = 0 #pseudo indicate excision if 1, insertion if 2
                    flag2 = 0 #ref indicate excision if 1, insertion if 2
                     
                    if unit[1] == 'HEG4':
                        #excision
                        if unit[2] == 'clipped' and unit[3] == 'clipped' and unit[4] == 'covered':
                            flag1 = 1
                            flag2 = 1
                        if (unit[2] == 'clipped' or unit[3] == 'clipped') and unit[4] == 'covered':
                            flag1 = 1
                            flag2 = 1
                        if unit[4] == 'covered':
                            flag2 = 1
                        #insertion
                        if unit[2] == 'covered' and unit[3] == 'covered' and unit[4] == 'clipped':
                            flag1 = 2
                            flag2 = 2
                        if (unit[2] == 'covered' or unit[3] == 'covered') and unit[4] == 'clipped':
                            flag1 = 2
                            flag2 = 2
                        if unit[4] == 'clipped': 
                            flag2 = 2
                        #status
                        if flag1 == 1 and flag2 == 1:
                            flag = 'Excision'
                        elif flag1 == 2 and flag2 == 2: 
                            flag = 'Insertion'
                        elif flag1 == 1 and flag2 == 2:
                            flag = 'Check'
                        elif flag1 == 2 and flag2 == 1:
                            flag = 'Check'
                        elif flag1 == 1 and flag2 == 0:
                            flag = 'Excision'
                        elif flag1 == 0 and flag2 == 1:
                            flag = 'Excision'
                        elif flag1 == 2 and flag2 == 0:
                            flag = 'Insertion'
                        elif flag1 == 0 and flag2 == 2:
                            flag = 'Insertion'
                        else:
                            flag = 'Unknown'

                        #conflic
                        if unit[2] == 'clipped' and unit[3] == 'covered':
                            flag = 'Check'
                        elif unit[2] == 'covered' and unit[3] == 'clipped':
                            flag = 'Check'
                        elif unit[2] == unit[4] and not unit[4] == 'unknown':
                            flag = 'Check'
                        elif unit[3] == unit[4] and not unit[4] == 'unknown':
                            flag = 'Check'

                    elif unit[1] == 'NB':
                        flag = 'NB_GT'
                    #print ril, unit[1], flag, flag1, flag2
                    unit.extend([flag, ping_number, ping_cd])
                else:
                    pass
                print >> ofile, ','.join(map(str, unit))
    ofile.close()


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--csv')
    parser.add_argument('-p', '--ping_code')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.csv) > 0
    except:
        usage()
        sys.exit(2)
    
    #outdir
    if not args.output:
        args.output = 'mPing_boundary_mPing_GT_Ping_code'
    createdir(args.output)

    #read ping code
    ping_code = read_ping(args.ping_code)
    
    #read and write csv
    csvs = glob.glob('%s/*.csv' %(args.csv))
    for csv in csvs:
        read_csv(csv, ping_code, args.output)

if __name__ == '__main__':
    main()

