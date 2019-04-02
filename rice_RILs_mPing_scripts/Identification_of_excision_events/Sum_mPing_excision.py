#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python Sum_excision_distance.py --dir mPing_boundary_mPing --distance ../input/mPing_dist.100kb.list.sorted --blacklist blacklist.ril

--dir: directory of results of mPing_Boundary_Coverage.py
--distance: distance between pairs of mPing and their strand
Chr3.29404858   Chr3.29404901   43      -       +

Summary the presence/absence of mPing of these pairs within 100kb from results of mPing_Boundary_Coverage.py.
The results will be table of format (+ is covered and - is clipped):
Pairs		RILs_count	++	+-	+-	--
mping1:mping2 	N		N	N	N	N
We will manually correct results of mPing_Boundary_Coverage.py and update this results agasin to see if consistent

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr1:6806763-6806761,Genotype_Bin,Pseudo_mPing_status_up,Pseudo_mPing_status_down,Ref_mPing_status
#RIL1,NB,clipped,unknown,unknown
#RIL2,HEG4,covered,covered,unknown
def readcsv(infile, mping_gt):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'RIL'): 
                unit = re.split(r',',line)
                if unit[1] == 'HEG4' and mping_gt == 'HEG4':
                    #both known and consistent
                    if unit[2] == unit[3] and unit[3] != 'unknown':
                        data[unit[0]] = unit[2]
                    #both unknonw
                    elif unit[2] == unit[3] and unit[3] == 'unknown':
                        if unit[4] != 'unknown':
                            ##use reference instead of pseudoref to determine the status
                            if unit[4] == 'covered':
                                data[unit[0]] = 'clipped'
                            elif unit[4] == 'clipped':
                                data[unit[0]] = 'covered'
                        else:
                            data[unit[0]] = 'unknown'
                    #both known and not euqal
                    elif unit[2] != 'unknown' and unit[3] != 'unknown':
                        if unit[4] != 'unknown':
                            ##use reference instead of pseudoref to determine the status
                            if unit[4] == 'covered':
                                data[unit[0]] = 'clipped'
                            elif unit[4] == 'clipped':
                                data[unit[0]] = 'covered'
                        else:
                            data[unit[0]] = 'unknown' 
                    elif unit[2] != 'unknown':
                        if unit[4] != 'unknown':
                            ##use reference instead of pseudoref to determine the status
                            if unit[4] == 'covered' and unit[2] == 'clipped':
                                data[unit[0]] = 'clipped'
                            elif unit[4] == 'clipped' and unit[2] == 'covered':
                                data[unit[0]] = 'covered'
                            else:
                                data[unit[0]] = 'unknown'
                        else:
                            data[unit[0]] = unit[2]
                    elif unit[3] != 'unknown':
                        if unit[4] != 'unknown':
                            ##use reference instead of pseudoref to determine the status
                            if unit[4] == 'covered' and unit[3] == 'clipped':
                                data[unit[0]] = 'clipped'
                            elif unit[4] == 'clipped' and unit[3] == 'covered':
                                data[unit[0]] = 'covered'
                            else:
                                data[unit[0]] = 'unknown'
                        else:
                            data[unit[0]] = unit[3] 
                    else:
                        if unit[4] != 'unknown':
                            ##use reference instead of pseudoref to determine the status
                            if unit[4] == 'covered':
                                data[unit[0]] = 'clipped'
                            elif unit[4] == 'clipped':
                                data[unit[0]] = 'covered'
                        else:
                            data[unit[0]] = 'unknown'
    return data

#RIL3,HEG4,unknown,clipped,covered,Excision,6,ABCEFG*
def readcsv_excision(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'RIL'): 
                unit = re.split(r',',line)
                #already take care of genotype in identifying excision (Ping_number_RILs.High_exicison_Ref.py)
                if unit[1] == 'NB' or unit[1] == 'HEG4' or unit[1] == 'NA':
                    data[unit[0]] = unit[5]
    return data

#read blacklist of problem RILs
def read_blacklist(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                ril  = re.sub(r'RIL', r'RIL', unit[0])
                data[ril] = 1
    return data

#Chr1:10903901-10903903	HEG4
#Chr1:11226135-11226137	HEG4
def read_mpinglist(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                unit[0] = re.sub(r':', r'_', unit[0])
                unit[0] = re.sub(r'-', r'_', unit[0])
                if unit[1] == 'HEG4':
                    data[unit[0]] = 'HEG4'
                else:
                    data[unit[0]] = 'NB'
    return data


def summary(directory, mpings, prefix):
    #ofile = open ('%s.excision_distance.list' %(prefix), 'w')
    #ofile1 = open ('%s.distance_accumulation_excision.list' %(prefix), 'w')
    ofile2 = open ('%s.mping_excision.list' %(prefix), 'w')
    #ofile3 = open ('%s.distance_accumulation_excision.kstest' %(prefix), 'w')0
    #print >> ofile1, 'Distance\tExcision\tAccumulated_Excision\tFraction_Accumulated_Excision\ttAccumulated_Control\tFraction_Accumulated_Control\tFraction_Control'
    print >> ofile2, 'mPing\t#Excision'
    
    #distance_int = 10000
    #distance_excision = defaultdict(lambda : int())
    #distance_control  = defaultdict(lambda : int())
    #mping_num_control = 0
    #excision_num_control = 0
    for mping in sorted(mpings.keys()):
        #print '%s, %s, %s, %s, %s, %s' %(pair, mpings[pair][0], mpings[pair][1], mpings[pair][2], mpings[pair][3], mpings[pair][4])
        #print mping
        #mping_status = readcsv('%s/%s.matrix.csv' %(directory, mping))
        if not os.path.isfile('%s/%s.matrix.csv' %(directory, mping)):
            print 'matrix.csv not found for %s' %(mping)
            continue
        mping_status = readcsv_excision('%s/%s.matrix.csv' %(directory, mping))
        excision = 0
        rils = [] 
        for ril in mping_status.keys():
            #if blacklist.has_key(ril):
            #    continue
            #if mping_status[ril] == 'clipped':
            if mping_status[ril] == 'Excision':
                excision += 1
                ril = re.sub(r'RIL', r'', ril)
                rils.append('RIL%s' %(ril))
        #print mping
        mping_temp = re.split('_', mping)
        mping_temp1= '%s:%s-%s' %(mping_temp[0], mping_temp[1], mping_temp[2])
        print >> ofile2, '%s\t%s\t%s\t%s\t%s' %(mping_temp1, mpings[mping][0], mpings[mping][1], excision, ','.join(rils))
        #if excision < 2000:
            #print >> ofile, '%s\t%s\t%s' %(mping, mpings[mping], excision)
            #index = int(float(mpings[mping])/distance_int)
            #distance_excision[index] += excision
            #distance_control[index] += 1
        #    mping_num_control += 1
    #excision_accum = 0
    #excision_control_accum = 0
    #for i in range(0, max(distance_excision.keys())):
    #    excision_n = distance_excision[i] if distance_excision.has_key(i) else 0
    #    excision_accum = excision_accum if excision_n == 0 else excision_accum + excision_n
    #    #control number is the mPing * average excison per mPing 
    #    excision_control = distance_control[i]*(float(sum(distance_excision.values()))/mping_num_control) #mping number * (total exicison/total mping)
    #    excision_control_accum += excision_control
    #    print >> ofile1, '%s\t%s\t%s\t%s\t%s\t%s\t%s' %((i+1)*distance_int, excision_n, excision_accum, float(excision_accum)/sum(distance_excision.values()), excision_control, excision_control_accum, float(excision_control_accum)/sum(distance_excision.values()))
    #ofile.close()
    #ofile1.close()
    ofile2.close()
    #ofile3.close() 

#Chr3.29404858   Chr3.29404901   43      -       +
def readtable(infile):
    data = defaultdict(lambda : int())
    pair = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping1 = re.split(r'\.', unit[5])
                mping2 = re.split(r'\.', unit[6])
                mping1_idx = '%s_%s_%s' %(mping1[0], mping1[1], mping1[2])
                mping2_idx = '%s_%s_%s' %(mping2[0], mping2[1], mping2[2])
                #mping1_idx = '%s_%s_%s' %(mping1[0], mping1[1], str(int(mping1[1]) + 2))
                #mping2_idx = '%s_%s_%s' %(mping2[0], mping2[1], str(int(mping2[1]) + 2))
                #mping1_idx = '%s_%s_%s' %(mping1[0], str(int(mping1[1]) + 2), mping1[1])
                #mping2_idx = '%s_%s_%s' %(mping2[0], str(int(mping2[1]) + 2), mping2[1])
                #if data.has_key(mping1_idx):
                #    print 'present in more than 1 pairs: %s' %(mping1_idx)
                #if data.has_key(mping1_idx):
                #    print 'present in more than 1 pairs: %s' %(mping2_idx)
                if data.has_key(mping1_idx):
                    data[mping1_idx] = int(unit[2]) if int(unit[2]) < data[mping1_idx] else data[mping1_idx]
                else:
                    data[mping1_idx] = int(unit[2])
                if data.has_key(mping2_idx):
                    data[mping2_idx] = int(unit[2]) if int(unit[2]) < data[mping2_idx] else data[mping2_idx]
                else:
                    data[mping2_idx] = int(unit[2])
                #data[mping1_idx] = 1
                #data[mping2_idx] = 1
                #index       = '%s:%s' %(unit[0], unit[1])
                #pair[index] = [mping1_idx, mping2_idx, unit[2], unit[3], unit[4]]
    #print 'linked mPing: %s' %(str(len(data.keys())))
    return data

def readgff(infile):
    data = defaultdict(lambda : str())
    #pt   = re.compile(r'Shared')
    pn   = re.compile(r'TE=(\w+?);') 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                name  = 'NA'
                unit  = re.split(r'\t',line)
                locus = '%s_%s_%s' %(unit[0], unit[3], unit[4])
                m = pn.search(unit[8])
                if m:
                    name = m.groups(0)[0]
                data[locus] = [name, 'Nonref']
                #if pt.search(unit[8]):
                #    data[locus] = [name, 'Shared']
                #else:
                #    data[locus] = [name, 'Ref']
                #print line
                #print '%s\t%s' %(locus, data[locus])
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-dir', '--dir')
    parser.add_argument('-g', '--gff')
    #parser.add_argument('-b', '--blacklist')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.dir) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = args.dir

    #read list of mpings and their genotypes (HEG4, Nipponbare or shared)
    #mpingt = read_mpinglist('HEG4_NB_ALL.mPing.list')
    #mpings = readtable(args.distance)
    #blacklist_ril = defaultdict(lambda : int())
    #if args.blacklist:
    #    blacklist_ril = read_blacklist(args.blacklist) 
    mpings = readgff(args.gff)
    summary(args.dir, mpings, args.project)

if __name__ == '__main__':
    main()

