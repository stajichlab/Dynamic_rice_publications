#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import subprocess

def usage():
    test="name"
    message='''
python footprint_events.py --input mping.excision.draw.example

Read list of mPing and excised RILs. Identify footprint for each RIL and summary excision events according to different footprints.
--input:
mPing   NumberOfRILs    RILs
Chr1:36267659-36267661  9       101,112,95,119,66,89,127,143,151
Chr1:36270511-36270514  5       122,125,63,71,98

    '''
    print message

def get_fasta_seq(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = str(record.seq)
    return fastaid

def read_blacklist(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                ril  = re.sub(r'RIL', r'', unit[0])
                data[ril] = 1
    return data


##default flanking of mping
def get_flank_seq(mping, flank_len):
    reg = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match = reg.search(mping)
    chro  = match.groups(0)[0]
    start = int(match.groups(0)[1]) - flank_len
    end   = int(match.groups(0)[2]) + flank_len
    refseq= get_fasta_seq('/rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa')
    flank = refseq[chro][start:end]
    return flank

##flanking of mping with footprint
def get_footprint_str(mping, flank, fp_dict, flank_len):
    reg = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match = reg.search(mping)
    chro  = match.groups(0)[0]
    start = int(match.groups(0)[1])
    end   = int(match.groups(0)[2])
    flank_str = list(flank)
    flank_ins = ''
    for fp_start in sorted(fp_dict.keys(), key=int):
        if fp_dict[fp_start][0] == 'D':
            for i in range(fp_start-start-(flank_len + 3), fp_start-start- (flank_len + 3) +int(fp_dict[fp_start][1])):
                flank_str[i] = '-'
        if fp_dict[fp_start][0] == 'I':
            flank_ins = ' '*(fp_start-start+ (flank_len - 3)) + '|' + fp_dict[fp_start][1]
    return [flank_ins, ''.join(flank_str)]

#mPing   NumberOfRILs    RILs
#Chr1:36267659-36267661  9       101,112,95,119,66,89,127,143,151
#Chr1:36270511-36270514  5       122,125,63,71,98
def readtable(infile, output, blacklist):
    flank_len = 50 #short length may introduce error when the deletion size are large
    ofile1 = open('%s.footprint.sequence.txt' %(output), 'w')
    ofile2 = open('%s.footprint.list.txt' %(output), 'w')
    print >> ofile2, 'mPing\tNumberOfEvents\tNumberOfEvents_Footprint\tNumberOfRILs\tNumberOfRILs_Footprint\tFootprint:RILs'
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'):
                data = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : list()))) 
                unit = re.split(r'\t',line)
                mping = unit[0]
                if unit[1] == '0' or unit[1] == 'NA':
                    #skip these without excision
                    print >> ofile2, '%s\t0\t0\t0\t0\t0' %(mping)
                    continue
                flank = get_flank_seq(mping, flank_len)
                rils  = re.split(r',', unit[2])
                print >> ofile1, '>%s' %(mping)
                print >> ofile1, '%s\t%s' %('{0:15}'.format('Nipponbare'), flank)
                for ril in rils:
                    ril = re.sub(r'RIL', r'', ril)
                    bam = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/GN%s.bam' %(ril)
                    if blacklist.has_key(ril):
                        continue
                    #data.update(subbam(mping, ril, bam, output))
                    subbam(mping, ril, bam, output, data)
                    fp_dict = data[mping][ril]
                    flank_ril = flank
                    fp_str  = get_footprint_str(mping, flank_ril, fp_dict, flank_len)
                    if fp_str[0] is not '':
                        print >> ofile1, '%s\t%s' %('{0:15}'.format(' '), fp_str[0])
                    else:
                        print >> ofile1, ' '
                    print >> ofile1, '%s\t%s' %('{0:15}'.format('RIL%s' %(ril)), fp_str[1])
                #sum footprint and excision events for mPing
                events = defaultdict(lambda : list())
                events_str = []
                rils_footprint = 0
                events_footprint = 0
                #print mping
                for ril in sorted(data[mping].keys(), key=int):
                    fp_dict = data[mping][ril]
                    ril_fp_type = []
                    fp_signature = 'Perfect'
                    if len(fp_dict.keys()) > 0:
                        #have footprints
                        rils_footprint += 1
                        for start in sorted(fp_dict.keys(), key=int):
                            fp_signature = '_'.join([str(start), fp_dict[start][0], fp_dict[start][1]])
                            ril_fp_type.append(fp_signature)
                            #print fp_signature
                    else:
                        ril_fp_type.append(fp_signature)
                    #print '%s\t%s\t%s' %(mping, ril, fp_signature)
                    events['-'.join(ril_fp_type)].append('RIL%s' %(ril))
                ##merge all footprint for one mPing, fp:ril,ril;fp:ril
                for t in sorted(events.keys()):
                    if not t == 'Perfect':
                        events_footprint += 1
                    events_str.append('%s:%s' %(t, ','.join(events[t])))
                if len(events.keys()) == 0:
                    print >> ofile2, '%s\t%s\t%s\t%s\t%s\t0' %(mping, len(events.keys()), events_footprint, len(rils), rils_footprint)
                else: 
                    print >> ofile2, '%s\t%s\t%s\t%s\t%s\t%s' %(mping, len(events.keys()), events_footprint, len(rils), rils_footprint, ';'.join(events_str))
    ofile1.close()
    ofile2.close()

#for one mping in ril, identify footprint, return by dict: mping->ril->footprint_startonref->[footprint_type, footprint_base]
def subbam(mping, ril, bam, output, data):
    #data = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : list())))
    reg = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match = reg.search(mping)
    start = int(match.groups(0)[1])
    end   = int(match.groups(0)[2])
    chro    = match.groups(0)[0]
    region  = '%s:%s-%s' %(chro, start-100000, end+100000)
    outdir  = './%s/%s' %(output, mping)
    if not os.path.exists(output):
        os.mkdir(output)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    test_bam = '%s/%s_%s.bam' %(outdir, ril, mping)
    if not os.path.isfile(test_bam):
        os.system('samtools view -hb %s %s > %s' %(bam, region, test_bam))
        os.system('samtools index %s' %(test_bam))

    cmd = 'samtools view %s %s' %(bam, mping)
    out = subprocess.check_output(cmd, shell=True)

    ##parse align around mPing
    #print 'mPing: %s' %(mping)
    #print '%s' %(cmd)
    total = 0
    footprint = []
    pattern = re.compile(r'([0-9]+)([A-Z]+)')
    lines = re.split(r'\n', out)
    for line in lines:
        unit = re.split(r'\t', line)
        if len(unit) < 2:
            continue
        #print 'Read: %s' %(unit[0])
        total += 1
        ref_start = int(unit[3])
        for (base, match) in re.findall(pattern, unit[5]):
            if match == 'M':
                ref_start += int(base)
            elif match == "I":
                ins_s   = ref_start - int(unit[3])
                ins_e   = ins_s + int(base)
                ins_seq = unit[9][ins_s:ins_e]
                footprint.append([ref_start, 'I', ins_seq])
            elif match == 'D':
                footprint.append([ref_start, 'D', base])
                ref_start += int(base)
    ##remove low frequency footprint
    temp_dict = defaultdict(lambda : int())
    fp_dict   = defaultdict(lambda : list())
    for fp in footprint:
        index = '_'.join(map(str, fp))
        temp_dict[index] += 1
    for i in range(len(footprint)):
        #print i
        index = '_'.join(map(str, footprint[i]))
        if temp_dict[index] > total*0.3:
            fp_dict[footprint[i][0]] = [footprint[i][1], footprint[i][2]]

    ##output footprint table
    for start in sorted(fp_dict.keys(), key=int): 
        #print '%s\t%s\t%s' %(start, fp_dict[start][0], fp_dict[start][1])
        data[mping][ril][start] = [fp_dict[start][0], fp_dict[start][1]] 

    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-b', '--blacklist')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = 'test_example_draw'
    
    blacklist_ril = defaultdict(lambda : int())
    if args.blacklist:
        blacklist_ril = read_blacklist(args.blacklist)
    
    readtable(args.input, args.output, blacklist_ril)

if __name__ == '__main__':
    main()

