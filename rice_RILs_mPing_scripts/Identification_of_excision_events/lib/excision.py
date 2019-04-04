#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import subprocess

def usage():
    test="name"
    message='''
python Excision.py --input ../input/10092013.mpings.gff --mode Non_ref

--input: input gff file of all insertions in RILs. input only non_ref inseriton if Non_ref mode. input reference if Reference mode. Shared for Shared mode
--mode:  detect excision for Non_ref, Reference or shared insertions. Non_ref or Reference or Shared

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
Chr1    RelocaTE        mPing   10903901        10903903
'''
def readmPing(gff, f):
    data = defaultdict(lambda: int)
    if int(f) == 0:
        data.update(readfile(gff, 0))
    elif int(f) == 1:
        data.update(readfile(gff, 1))
    elif int(f) == 2:
        data.update(readfile(gff, 2)) 
    return data

def readfile(infile, flag):
    data = defaultdict(lambda: int)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = '%s:%s-%s' %(unit[0], unit[3], unit[4])
                data[mping] = flag
    return data

'''
Chr1    RIL10_0 RelocaTE        1132975 1132977 .       .       .       ID=mPing_1;Strain=RIL10_0;TSD=TAA;
'''
def readtable(infile):
    data = defaultdict(lambda: defaultdict(lambda: int()))
    p = re.compile(r'Strain=RIL(\d+)\_\d+;')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = '%s:%s-%s' %(unit[0], unit[3], unit[4])
                m = p.search(unit[8])
                strain = m.groups(0)[0] if m else 'NA'
                #print strain, mping, line
                data[strain][mping] = 1
    return data


'''
Chr1    RIL10_0 RelocaTE        1132975 1132977 .       .       .       ID=mPing_1;Strain=RIL10_0;TSD=TAA;
'''
def mping_frequency(infile):
    data = defaultdict(lambda: defaultdict(lambda: int))
    rils = defaultdict(int)
    inf  = defaultdict()
    p = re.compile(r'Strain=RIL(\d+)\_\d+;')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = '%s:%s-%s' %(unit[0], unit[3], unit[4])
                m = p.search(unit[8])
                strain = m.groups(0)[0] if m else 'NA'
                data[mping][strain] = 1
                rils[strain] =1
                inf[mping] = [unit[0], unit[3], unit[4]]
    total = len(rils.keys())
    #print total
    mping_frq = defaultdict(lambda: float)
    for m in data.keys():
        count = len(data[m].keys())
        frq = float(count)/total
        #print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(inf[m][0], inf[m][1], inf[m][2], m, '+', str(count), str(frq))
        mping_frq[m] = frq
    return mping_frq


#Convert BIN MAP
'''
""      "GN1"   "GN10"  "GN100" "GN101" "GN102" "GN103" "GN104" "GN105" "GN106" "GN107" "GN108
"0100222046"    1       1       0       0       0       0       1       1       1       1
"0100500860"    1       1       0       0       0       0       1       1       1       1
'''

def convert_MAP(infile):
    rils = []
    data = defaultdict(lambda : defaultdict(lambda: defaultdict(lambda: str)))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            line = line.replace('"', '')
            if line[0:1].isdigit():
                unit = re.split(r'\t',line)
                #print '%s\t%s\t%s' %(chrs, int(unit[0][2:]), str(chr_start[chrs]))
                chrs = 'Chr%s' %(str(int(unit[0][0:2])))
                pos = int(unit[0][2:])
                for i in range(1,len(unit)):
                    #print i, rils[i], chrs, pos
                    ril = rils[i]
                    data[ril][chrs][pos] = unit[i]
            else:
                unit = re.split(r'\t',line)
                unit[0] = 'Position'
                #print unit[0], unit[1], unit[2]
                rils.extend(unit)
   
    #for t in sorted(data['GN204']['Chr10'].keys(), key=int):
    #    print t
    return data

#Convert SNP map
#""      "GN1"   "GN10"  "GN100" "GN101" "GN102" "GN103" "GN104" "GN105" "GN106" 
#"0100021547A"   NA      NA      NA      0       0       0       NA      NA      
#"0100031071A"   NA      1       0       0       0       0       1       1       
#"0100031478C"   1       1       0       0       0       0       NA      1       

def convert_MAP_SNP(infile):
    rils = []
    data = defaultdict(lambda : defaultdict(lambda: defaultdict(lambda: str)))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            line = line.replace('"', '')
            if line[0:1].isdigit():
                unit = re.split(r'\t',line)
                #print '%s\t%s\t%s' %(chrs, int(unit[0][2:]), str(chr_start[chrs]))
                chrs = 'Chr%s' %(str(int(unit[0][0:2])))
                unit[0] = unit[0][:-1] #remove reference base
                pos = int(unit[0][2:])
                for i in range(1,len(unit)):
                    #print i, rils[i], chrs, pos
                    ril = rils[i]
                    data[ril][chrs][pos] = unit[i]
            else:
                unit = re.split(r'\t',line)
                unit[0] = 'Position'
                #print unit[0], unit[1], unit[2]
                rils.extend(unit)
   
    #for t in sorted(data['GN204']['Chr10'].keys(), key=int):
    #    print t
    return data



def binarySearch(data, val):
    highIndex = len(data)-1
    lowIndex = 0
    while highIndex > lowIndex:
            index = (highIndex + lowIndex) / 2
            sub = int(data[index])
            #print highIndex, index, lowIndex, sub, val
            if data[lowIndex] == val:
                    return [lowIndex, lowIndex]
            elif sub == val:
                    return [index, index]
            elif data[highIndex] == val:
                    return [highIndex, highIndex]
            elif sub > val:
                    if highIndex == index:
                            return sorted([highIndex, lowIndex])
                    highIndex = index
            else:
                    if lowIndex == index:
                            return sorted([highIndex, lowIndex])
                    lowIndex = index
    return sorted([highIndex, lowIndex])
 

def findbin(start, binmap, ril, chrs):
    
    array = []
    array.extend(sorted(binmap[ril][chrs].keys(), key=int))
    #mping after last bin on chromosome, return 0 mean genotype unknown
    if int(start) > int(array[-1]):
        return 0
    index = binarySearch(array, int(start))
    #bin overlap with mping is a transition bin on recombination block, which mean genotype may not be precise, return 0
    #if transition[ril][chrs].has_key(array[index[1]]):
    #    return 0
    #return the bin that overlap with mping
    #if int(start) == 6566243 and ril == 'GN204':
    #    for t in sorted(binmap[ril][chrs].keys(), key=int):
    #        print t
    #print start, array[index[0]], array[index[1]]
    return array[index[1]]

#[1,0,0,0,0]
#get genotype of five snp, return genotype and frequency
def snp_type(snp_index, snpmap, ril, chrs):
    data = defaultdict(lambda : int())
    for i in snp_index:
        s = snpmap[ril][chrs][i]
        if not s == 'NA':
            data[s] += 1
            print 'S: %s; Data[s]: %s' %(s, data[s])
    genotypes = sorted(data, key=data.get, reverse=True)
    if len(genotypes) >= 1:
        print 'GT: %s' %('\t'.join(map(str, genotypes)))
        print 'GT %s Rate %s' %(genotypes[0], float(data[genotypes[0]])/sum(data.values()))
        return [genotypes[0], float(data[genotypes[0]])/sum(data.values())]
    else:
        return [3, 0]

#def findbin_SNP(start, binmap, snpmap, ril, chrs):
#    print 'start: %s' %(start)
#    a = [0, 'NA']
#    return a

def findbin_SNP(start, binmap, snpmap, ril, chrs):
    array = []
    array.extend(sorted(binmap[ril][chrs].keys(), key=int))
    array_snp = []
    array_snp.extend(sorted(snpmap[ril][chrs].keys(), key=int))
    #mping after last bin on chromosome, return 0 mean genotype unknown
    if int(start) > int(array[-1]):
        return [0, 'NA']
    index = binarySearch(array, int(start))
    index_snp = binarySearch(array_snp, int(start))
    
    print 'Find Bin %s %s' %(index[0], index[1])
    print 'Find SNP %s %s' %(index_snp[0], index_snp[1])

    #check five flanking SNP if consistent
    pos_gt = []
    backward_snp_idx = array_snp[(index_snp[0]-4):(index_snp[0]+1)]
    forward_snp_idx  = array_snp[index_snp[1]:(index_snp[1]+5)]
    gt_1 = snp_type(backward_snp_idx, snpmap, ril, chrs)
    gt_2 = snp_type(forward_snp_idx, snpmap, ril, chrs)
    
    print 'backward SNP %s' %(backward_snp_idx)
    print 'forwad   SNP %s' %(forward_snp_idx)
    if len(backward_snp_idx) >= 1:
        print '5SNP %s %s' %(backward_snp_idx[-1], snpmap[ril][chrs][backward_snp_idx[-1]])
    else:
        print '5SNP NA NA'
    if len(forward_snp_idx) >= 1:
        print '3SNP %s %s' %(forward_snp_idx[0], snpmap[ril][chrs][forward_snp_idx[0]])
    else:
        print '3SNP NA NA'

    #block genotype unknown, return
    print 'RIL: %s' %(ril)
    print 'Chrs %s Start: %s' %(chrs, start)
    print 'GT upstream: %s' %(gt_1)
    print 'GT downstream: %s' %(gt_2)
    print 'Index:%s %s' %(index[0], index[1])
    print 'Bin %s %s' %(array[index[0]], array[index[1]])
    print 'Bin Genotype: %s %s' %(binmap[ril][chrs][array[index[0]]], binmap[ril][chrs][array[index[1]]])
    if str(binmap[ril][chrs][array[index[1]]]) == 'NA':
        pos_gt.extend([array[index[1]], 'NA'])
    else:
        #block genotype known, use snp to comfirm
        if gt_1[0] == gt_2[0] and str(gt_1[0]) == str(binmap[ril][chrs][array[index[1]]]):
            #top SNP genotype consistent with bin genotype, then check the flanking SNP to remove any possibility of recombination at mPing loci
            if str(snpmap[ril][chrs][backward_snp_idx[-1]]) == str(binmap[ril][chrs][array[index[1]]]) and str(snpmap[ril][chrs][forward_snp_idx[0]]) == str(binmap[ril][chrs][array[index[1]]]):
                print 'Consistent: SNP genotype %s; Bin genotype %s' %(gt_1[0], binmap[ril][chrs][array[index[1]]])
                #snp consistent with genotype block
                pos_gt.extend([array[index[1]], gt_1[0]])
            elif str(snpmap[ril][chrs][backward_snp_idx[-1]]) == 'NA':
                # 5' flank SNP is NA, we use the second last SNP to compare with bin genotype
                if str(snpmap[ril][chrs][backward_snp_idx[-2]]) == str(binmap[ril][chrs][array[index[1]]]) and str(snpmap[ril][chrs][forward_snp_idx[0]]) == str(binmap[ril][chrs][array[index[1]]]):
                    pos_gt.extend([array[index[1]], gt_1[0]])
                else:                
                    print 'Inconsistent: flanking SNP genotype %s %s; Bin genotype %s' %(str(snpmap[ril][chrs][backward_snp_idx[-2]]), str(snpmap[ril][chrs][forward_snp_idx[0]]), binmap[ril][chrs][array[index[1]]])
                    pos_gt.extend([array[index[1]], 'NA'])
            elif str(snpmap[ril][chrs][forward_snp_idx[0]]) == 'NA':
                # 3' flank SNP is NA, we use the second last SNP to compare with bin genotype
                if str(snpmap[ril][chrs][backward_snp_idx[-1]]) == str(binmap[ril][chrs][array[index[1]]]) and str(snpmap[ril][chrs][forward_snp_idx[1]]) == str(binmap[ril][chrs][array[index[1]]]):
                    pos_gt.extend([array[index[1]], gt_1[0]])
                else:
                    print 'Inconsistent: flanking SNP genotype %s %s; Bin genotype %s' %(str(snpmap[ril][chrs][backward_snp_idx[-1]]), str(snpmap[ril][chrs][forward_snp_idx[1]]), binmap[ril][chrs][array[index[1]]])
                    pos_gt.extend([array[index[1]], 'NA'])
            else:
                print 'Inconsistent: flanking SNP genotype %s %s; Bin genotype %s' %(str(snpmap[ril][chrs][backward_snp_idx[-1]]), str(snpmap[ril][chrs][forward_snp_idx[0]]), binmap[ril][chrs][array[index[1]]])
                pos_gt.extend([array[index[1]], 'NA'])
        else:
            #snp inconsistent with genotype block
            print 'Inconsistent: SNP genotype %s; Bin genotype %s' %(gt_1[0], binmap[ril][chrs][array[index[1]]])
            pos_gt.extend([array[index[1]], 'NA'])
    return pos_gt



def genotyping(ril, mping, binmap):
    ril = 'GN%s' %(ril)
    p = re.compile(r'(\w+):(\d+)\-(\d+)')
    m = p.search(mping)
    chrs = ''
    start = 0
    end   = 0
    if m:
        chrs  = m.groups(0)[0]
        start = m.groups(0)[1]
        end   = m.groups(0)[2]
    #print ril, chrs, start, len(binmap[ril][chrs].keys())
    pos = findbin(start, binmap, ril, chrs)
    #print 'find bin', ril, chrs, start, pos
    #pos is 0 when genotype of bin that overlap with mping is unknown
    #we should use raw bin not filled or uniq bin here and check if the genotype is NA.
    #genotype = binmap[ril][chrs][pos] if pos > 0 else '3'
    if pos > 0 and binmap[ril][chrs][pos] != 'NA':
        genotype = binmap[ril][chrs][pos]
    else:
        genotype = '3'
    return genotype

def genotyping_SNP(ril, mping, binmap, snpmap):
    ril = 'GN%s' %(ril)
    p = re.compile(r'(\w+):(\d+)\-(\d+)')
    m = p.search(mping)
    chrs = ''
    start = 0
    end   = 0
    if m:
        chrs  = m.groups(0)[0]
        start = m.groups(0)[1]
        end   = m.groups(0)[2]
    #print ril, chrs, start, len(binmap[ril][chrs].keys())
    pos, genotype = findbin_SNP(start, binmap, snpmap, ril, chrs)
    #print 'find bin', ril, chrs, start, pos
    #pos is 0 when genotype of bin that overlap with mping is unknown
    #we should use raw bin not filled or uniq bin here and check if the genotype is NA.
    #genotype = binmap[ril][chrs][pos] if pos > 0 else '3'
    #print 'pos: %s, genotype: %s' %(pos, genotype)
    if pos > 0 and genotype != 'NA':
        return genotype
    else:
        genotype = '3'


def validmap(binmap):
    last = 0
    data = defaultdict(lambda : defaultdict(lambda: defaultdict(lambda: int))) 
    for ril in sorted(binmap.keys()):
        for chrs in sorted(binmap[ril].keys()):
            for pos in sorted(binmap[ril][chrs].keys(), key=int):
                if binmap[ril][chrs][pos] != binmap[ril][chrs][last]:
                    data[ril][chrs][last] = 1
                    data[ril][chrs][pos] = 1
                    last = pos
                    #print ril, chrs, pos
    return data

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
    l_flag  = bamcheck_simple(bam, l_mping, bamck_file, 5)
    r_flag  = bamcheck_simple(bam, r_mping, bamck_file, 3)
    region  = '%s:%s-%s' %(chro, l_start-500, r_end+500)  
    test_bam = './test_bam/%s_%s.bam' %(ril, mping)
    os.system('samtools view -hb %s %s > %s' %(bam, region, test_bam))
    os.system('samtools index %s' %(test_bam))
    #flag is for excision of non_ref insertion. But we can use here to check breakpoint reads
    #flag: 0 excision with footprint, 1 support clip reads of non_ref insertion , 2 coverage too low, 3 other case, 4 excision with perfect breakpoint
    print 'l_flag: %s; r_flag: %s' %(l_flag, r_flag)
    if l_flag == 4 and r_flag == 4:
        ##read cover both boundary of ref insertion, suggesting insertion exists
        return 0
    elif l_flag == 4 or r_flag == 4:
        ##read cover one boundary of ref insertion, suggesting insertion exists
        return 0
    elif l_flag == 0 and r_flag == 0:
        ##read cover both boundary of ref insertion, suggesting insertion exists
        return 0
    elif l_flag == 0 or r_flag == 0:
        ##read cover one boundary of ref insertion, suggesting insertion exists
        return 0
    elif l_flag == 1 or r_flag == 1:
        ##cliped reads at boundary, suggesting excision of reference insertions
        return 1
    elif l_flag == 2 and r_flag == 2:
        ##coverage too low
        return 2
    else:
        return 3
   
def bamcheck(bam, mping, bamck_file, ril):
    reg = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match = reg.search(mping)
    start = int(match.groups(0)[1])
    end   = int(match.groups(0)[2])

    chro    = match.groups(0)[0]
    region  = '%s:%s-%s' %(chro, start-500, end+500)
    mping_5 = '%s:%s-%s' %(chro, start-5, end+5) 
    #test_bam = './test_bam/%s_%s.bam' %(ril, mping)
    #os.system('samtools view -hb %s %s > %s' %(bam, region, test_bam))
    #os.system('samtools index %s' %(test_bam))

    cmd = '/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools view %s %s' %(bam, mping_5)
    ofile = open(bamck_file, 'a') 
    print >> ofile, bam, mping
    print >> ofile, cmd
    out = subprocess.check_output(cmd, shell=True)
    lines = re.split(r'\n', out)
    total = 0
    total_count = 0
    covered = 0
    clipped = 0
    footprint = 0
    flag = 0 # flag indicate there are indel and softclips in the reads, probably solfclip only, supports for mping insertion
    if len(lines) < 2:
        return 2
    pattern = re.compile(r'([0-9]+)([A-Z]+)')
    for line in lines:
        print >> ofile, line
        total += 1
        unit = re.split(r'\t', line)
        if len(unit) < 2:
            continue
        if int(unit[4]) < 29:
            continue
        matches = []
        indel = 0
        soft  = 0
        for (base, match) in re.findall(pattern, unit[5]):
            if match == 'I' or match == 'D':
                indel += 1
                continue
            elif match == 'S' and int(base) > 10:
                soft += 1
                matches.append([base, match])
            else:
            #print >> ofile, base, match
                matches.append([base, match])
            #print >> ofile, matches[0]
        # merge neighbor matches that have some types after removing I and D
        j = 0
        for i in range(1, len(matches)):
            i = i - j
            if matches[i][1] == matches[i-1][1]:
                matches[i-1][0] = str(int(matches[i][0]) + int(matches[i-1][0]))
                del matches[i]
                j = j + 1
        
        print >> ofile, matches
        #flank should not be too long so we won't miss good reads that covers the junction. There is no TE problems here.
        flank = 10
        if indel > 0 and soft > 0:
            flag += 1
        elif indel > 0:
            footprint += 1
        if len(matches) == 1 and matches[0][1] == 'M': # perfect match
            if int(unit[3]) < start - flank and int(unit[3]) + int(matches[0][0]) > end + flank: # read cover mping insertion site
                covered += 1
        elif len(matches) == 2: # may have soft clip
            if matches[0][1] == 'S' and int(matches[0][0]) > 10 and matches[1][1] == 'M' and int(matches[1][0]) > 10: # clip at start
                if int(unit[3]) > start - flank and int(unit[3]) < end + flank: # read at mping insertion site, so the clip should be due to mping insertion
                    #covered += 1
                    clipped += 1
                elif int(unit[3]) <= start - flank and int(unit[3]) + int(matches[1][0]) >= end + flank: # read cover mping insertion site
                    covered += 1
            elif matches[1][1] == 'S' and int(matches[1][0]) > 10 and matches[0][1] == 'M' and int(matches[0][0]) > 10: # clip at end
                if int(unit[3]) + int(matches[0][0]) > start - flank and int(unit[3]) + int(matches[0][0]) < end + flank: # read at mping insertion site, so the clip should be due to mping insertion
                    #covered += 1
                    clipped += 1
                elif int(unit[3]) <= start - flank and int(unit[3]) + int(matches[0][0]) >= end + flank: # read cover mping insertion site
                    covered += 1
        elif len(matches) == 3 and matches[0][1] == 'S' and matches[1][1] == 'M' and matches[2][1] == 'S': # may have soft clip, but the other end of reads have clip too
            if int(unit[3]) > start - flank and int(unit[3]) < end + flank and int(matches[0][0]) > 10: # read at mping insertion site, so the clip should be on the left if due to mping
                #covered += 1
                clipped += 1
            elif int(unit[3]) + int(matches[1][0]) > start - flank and int(unit[3]) + int(matches[1][0]) < end + flank and int(matches[2][0]) > 10: # read start before mping insertion site, but clipped at mping
                #covered += 1
                clipped += 1
            elif int(unit[3]) <= start - flank and int(unit[3]) + int(matches[1][0]) >= end + flank: # read cover the mping insertion site
                covered += 1
        print >> ofile, covered, clipped
 
    #Use only covered and clipped reads for each junctions. These are useful reads. Other reads are not covering the junction or TEs. 
    total_count = covered + clipped
    total       = total_count
    if int(total) == 0:
        return 2 # coverage too low
    print >> ofile, covered, clipped
    print >> ofile, total, footprint, flag, float(float(footprint)/total)
    rate = float(float(clipped)/covered) if covered > 0 else 0
    if total <= 2:
        print >> ofile, 2
        return 2 # coverage too low
    #Determine covered reads first, which support excisions (our target). 
    #Even there are clipped reads (het with both covered and clipped reads) we call these excisions.
    elif float(float(covered)/total) >= 0.3:
        print >> ofile, 0
        if footprint >= 3:
            return 0 # support mping excision, many read cover the breapoint suggests precise excision, footprint
        else:
            return 4 # support mping excision, many read cover the breapoint suggests precise excision, perfect
    elif float(float(clipped)/total) >= 0.3:
        print >> ofile, 1
        return 1 # support mping insertion with clipped reads
    elif total > 2 and clipped < 1 and covered < 1 and flag < 1 and footprint < 1:
        print >> ofile, 2 # no supporting reads for anything
        return 2 # 
    #elif float(float(footprint)/total) > 0.3 and flag == 0:
    #    print >> ofile, 0
    #    return 0 # support mping excision, indel is possibly footprint of excision
    #elif float(float(covered)/total) >= 0.3:
	#print >> ofile, 0
        #if footprint >= 3:
        #    return 0 # support mping excision, many read cover the breapoint suggests precise excision, footprint
        #else:
        #    return 4 # support mping excision, many read cover the breapoint suggests precise excision, perfect
    else:
        print >> ofile, 3
        return 3 # other case

 
def bamcheck_simple(bam, mping, bamck_file, orientation):
    reg = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match = reg.search(mping)
    chro = match.groups(0)[0]
    start = int(match.groups(0)[1])
    end   = int(match.groups(0)[2])
    mping_5 = '%s:%s-%s' %(chro, start-5, end+5)
    cmd = '/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools view %s %s' %(bam, mping_5)
    ofile = open(bamck_file, 'a') 
    print >> ofile, bam, mping
    print >> ofile, cmd
    out = subprocess.check_output(cmd, shell=True)
    lines = re.split(r'\n', out)
    total = 0
    total_count = 0
    covered = 0
    clipped = 0
    footprint = 0
    flag = 0 # flag indicate there are indel and softclips in the reads, probably solfclip only, supports for mping insertion
    if len(lines) < 2:
        return 2
    pattern = re.compile(r'([0-9]+)([A-Z]+)')
    for line in lines:
        #print >> ofile, line
        #total += 1
        unit = re.split(r'\t', line)
        if len(unit) < 2:
            continue
        #skip if map quality is low, test for excision of reference insertion
        if int(unit[4]) < 29:
            continue
        print >> ofile, line
        total += 1
        matches = []
        indel = 0
        soft  = 0
        for (base, match) in re.findall(pattern, unit[5]):
            if match == 'I' or match == 'D':
                indel += 1
                continue
            elif match == 'S' and int(base) > 10:
                soft += 1
                matches.append([base, match])
            else:
            #print >> ofile, base, match
                matches.append([base, match])
            #print >> ofile, matches[0]
        # merge neighbor matches that have some types after removing I and D
        j = 0
        for i in range(1, len(matches)):
            i = i - j
            if matches[i][1] == matches[i-1][1]:
                matches[i-1][0] = str(int(matches[i][0]) + int(matches[i-1][0]))
                del matches[i]
                j = j + 1
          
        print >> ofile, matches
        #This flank need to be longer. As TE inserted in the pseudogenome we do not want to get incorrectly aligned reads from TE or flanking sequences.
        flank = 20
        if indel > 0 and soft > 0:
            flag += 1
        elif indel > 0:
            footprint += 1
        if len(matches) == 1 and matches[0][1] == 'M': # perfect match
            if int(unit[3]) < start - flank and int(unit[3]) + int(matches[0][0]) > end + flank: # read cover mping insertion site
                covered += 1
        elif len(matches) == 2: # may have soft clip
            if matches[0][1] == 'S' and int(matches[0][0]) > 10 and matches[1][1] == 'M' and int(matches[1][0]) > 10: # clip at start
                if int(unit[3]) > start - flank and int(unit[3]) < end + flank: # read at mping insertion site, so the clip should be due to mping insertion
                    #covered += 1
                    #clipped += 1
                    if int(orientation) == 5: 
                        #this is the 5' junction of mPing, these reads are clipped reads of mPing elements
                        continue
                    elif int(orientation) == 3:
                        #this is the 3' junction of mPing, these reads are clipped reads of flanking sequences
                        clipped += 1
                elif int(unit[3]) <= start - flank and int(unit[3]) + int(matches[1][0]) >= end + flank: # read cover mping insertion site
                    covered += 1
            elif matches[1][1] == 'S' and int(matches[1][0]) > 10 and matches[0][1] == 'M' and int(matches[0][0]) > 10: # clip at end
                if int(unit[3]) + int(matches[0][0]) > start - flank and int(unit[3]) + int(matches[0][0]) < end + flank: # read at mping insertion site, so the clip should be due to mping insertion
                    #covered += 1
                    #clipped += 1
                    if int(orientation) == 5:
                        #this is the 5' junction of mPing, these reads are clipped reads of flanking sequences
                        clipped += 1
                    elif int(orientation) == 3:
                        #this is the 3' junction of mPing, these reads are clipped reads of mPing elements
                        continue
                elif int(unit[3]) <= start - flank and int(unit[3]) + int(matches[0][0]) >= end + flank: # read cover mping insertion site
                    covered += 1
        elif len(matches) == 3 and matches[0][1] == 'S' and matches[1][1] == 'M' and matches[2][1] == 'S': # may have soft clip, but the other end of reads have clip too
            if int(unit[3]) > start - flank and int(unit[3]) < end + flank and int(matches[0][0]) > 10: # read at mping insertion site, so the clip should be on the left if due to mping
                #covered += 1
                #clipped += 1
                if int(orientation) == 5:
                    #this is the 5' junction of mPing, these reads are clipped reads of mPing elements
                    continue
                elif int(orientation) == 3:
                    #this is the 3' junction of mPing, these reads are clipped reads of flanking sequences
                    clipped += 1 
            elif int(unit[3]) + int(matches[1][0]) > start - flank and int(unit[3]) + int(matches[1][0]) < end + flank and int(matches[2][0]) > 10: # read start before mping insertion site, but clipped at mping
                #covered += 1
                #clipped += 1
                if int(orientation) == 5:
                    #this is the 5' junction of mPing, these reads are clipped reads of flanking
                    clipped += 1
                elif int(orientation) == 3:
                    #this is the 3' junction of mPing, these reads are clipped reads of mPing elements
                    continue 
            elif int(unit[3]) <= start - flank and int(unit[3]) + int(matches[1][0]) >= end + flank: # read cover the mping insertion site
                covered += 1
        print >> ofile, covered, clipped

    #Use only covered and clipped reads for each junctions. These are useful reads. Other reads are not covering the junction or TEs.
    total_count = covered + clipped
    total       = total_count 
    if int(total) == 0:
        return 2 # coverage too low
    print >> ofile, covered, clipped
    print >> ofile, total, covered, clipped, flag, float(float(covered)/total), float(float(clipped)/total)
    rate = float(float(clipped)/covered) if covered > 0 else 0
    if total <= 2:
        print >> ofile, 2
        return 2 # coverage too low
    #here we determine clipped reads first, which support excision.
    #even there are covered reads supporting insertion we still call these excisions.
    elif float(float(clipped)/total) >= 0.3:
        print >> ofile, 1
        return 1 # support mping insertion with clipped reads
    elif total > 2 and clipped < 1 and covered < 1 and flag < 1 and footprint < 1:
        print >> ofile, 2 # no supporting reads for anything
        return 2 # 
    #elif float(float(footprint)/total) > 0.3 and flag == 0:
    #    print >> ofile, 0
    #    return 0 # support mping excision, indel is possibly footprint of excision
    elif float(float(covered)/total) >= 0.3:
        print >> ofile, 0
        if footprint >= 3:
            return 0 # support mping excision, many read cover the breapoint suggests precise excision, footprint
        else:
            return 4 # support mping excision, many read cover the breapoint suggests precise excision, perfect
    else:
        print >> ofile, 3
        return 3 # other case



def excision(mPing_ancestor, mPing_rils, mPing_frq, num_file, bamck_file):
    mPing_excision = defaultdict(lambda : defaultdict(lambda : int()))
    binmap = convert_MAP('MPR.geno.bin')
    #transition = validmap(binmap)
    print 'Find excision: %s' %(len(mPing_ancestor.keys()))
    ofile_num = open(num_file, 'w')
    print >> ofile_num, 'mPing\tPresence\tAbsence\tUnknown\tFootprint'
    for mping in sorted(mPing_ancestor.keys()):
        print 'mPing: %s' %(mping)
        num_present   = 0
        num_unsure    = 0
        num_absent    = 0
        num_excision  = 0
        ril_absent    = []
        ril_unsure    = []
        ril_footprint = []
        num_footprint = 0 # reference insertion, footprint can not be determined, set to 0
        if not mPing_frq.has_key(mping) or mPing_frq[mping] < 0.05:
            continue
        else:
            print 'frequency: %s' %(mPing_frq[mping])
        #print 'in rils'
        for ril in sorted(mPing_rils.keys()):
            #print 'TT',ril,len(binmap[ril].keys())
            if not binmap.has_key('GN%s' %(ril)):
                continue
            ##genotype: 0 ref, 1 non_ref, 3 unknown
            genotype = genotyping(ril, mping, binmap)
            bam = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/GN%s.bam' %(ril)
            flag = 2 # 0 mean footprint, 4 mean perfect
            #print 'Check:', mping, ril, genotype
            #if mping == 'Chr1:1715117-1715119' and int(ril) == 10:
                #print 'EX', mping, ril, genotype, mPing_rils[ril][mping], mPing_ancestor[mping]
            if (int(genotype) == 0 and int(mPing_ancestor[mping]) == 0): # ril has genotype of reference and mping has genotype of reference (not include shared). Check excision of referecne only insertion
                #print 'S1', mPing_rils[ril].has_key(mping)
                if not mPing_rils[ril].has_key(mping):
                    if os.path.isfile(bam):
                        ##bamcheck need to changed to check reference insertion, bamcheck_ref
                        flag = bamcheck_ref(bam, mping, bamck_file, ril)
                        if flag == 1: ## bam check showed no mping insertion in this ril, excision
                            mPing_excision[mping][ril] = 1
                            num_absent += 1
                            ril_absent.append(ril)
                        elif flag == 2:
                            num_unsure += 1
                            ril_unsure.append(ril)
                        elif flag == 3:
                            num_unsure += 1
                            ril_unsure.append(ril)
                        elif flag == 0:
                            num_present += 1
                    #print mPing_excision[mping]
                else:
                    num_present += 1
            elif int(mPing_ancestor[mping]) == 2: #mping is shared, we do not care about genotype. becasue it is there anyway.
                if not mPing_rils[ril].has_key(mping):
                    if os.path.isfile(bam):
                        ##bamcheck need to changed to check reference insertion, bamcheck_ref
                        print 'checking in bam'
                        print 'Check:', mping, ril, genotype
                        flag = bamcheck_ref(bam, mping, bamck_file, ril)
                        if flag == 1: ## bam check showed no mping insertion in this ril, excision
                            mPing_excision[mping][ril] = 1
                            num_absent += 1
                            ril_absent.append(ril)
                            print 'excision'
                        elif flag == 2:
                            num_unsure += 1
                            ril_unsure.append(ril)
                            print 'coverage too low'
                        elif flag == 3:
                            num_unsure += 1
                            ril_unsure.append(ril)
                            print 'other case'
                        elif flag == 0:
                            num_present += 1
                            print 'have reference insertion'
                else:
                    num_present += 1
            elif (int(genotype) == 1 and int(mPing_ancestor[mping]) == 1): # ril has genotype of nonref and mping has genotype of nonref. Check excision for non_reference insertion
                #print 'S2', mPing_rils[ril].has_key(mping)
                if not mPing_rils[ril].has_key(mping):
                    if os.path.isfile(bam):
                        flag = bamcheck(bam, mping, bamck_file, ril)
                        if flag == 0: ## bam check showed no mping insertion in this ril
                            mPing_excision[mping][ril] = 1
                            num_absent += 1
                            ril_absent.append(ril)
                            num_footprint += 1
                            ril_footprint.append(ril)
                        elif flag == 4:
                            mPing_excision[mping][ril] = 2
                            num_absent += 1
                            ril_absent.append(ril)
                        elif flag == 2:
                            ##coverage too low
                            num_unsure += 1
                            ril_unsure.append(ril)
                        elif flag == 3:
                            ##other case
                            num_unsure += 1
                            ril_unsure.append(ril)
                        elif flag == 1:
                            ##have insertion
                            num_present += 1
                    #print mPing_excision[mping]
        print >> ofile_num, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(mping, num_present, num_absent, num_unsure, num_footprint, ','.join(ril_absent), ','.join(ril_unsure), ','.join(ril_footprint))
    return mPing_excision

def preplot(frq, ancestor, excision, output):
    ofile = open(output, 'w')
    for m in excision.keys():
        print '>', m, len(excision[m].keys())
        print >> ofile, m, len(excision[m].keys()), frq[m]
        for i in sorted(excision[m].keys()):
            print i, excision[m][i]
    ofile.close()

'''
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-m', '--mode')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    
    if not args.mode:
        args.mode = 'Non_Ref'

    if args.mode == 'Non_Ref':
        mPing_ancestor = readmPing('HEG4.ALL.mping.non-ref.gff', 1)
        mPing_rils     = readtable(args.input)
        mPing_frq      = mping_frequency(args.input)
        mPing_excision = excision(mPing_ancestor, mPing_rils, mPing_frq, 'mping.excision.non_ref.number', 'mping.excision.non_ref.bamcheck')
        preplot(mPing_frq, mPing_ancestor, mPing_excision, 'mping.excision.non_ref.table')
    elif args.mode == 'Ref' or args.mode == 'Reference':
        mPing_ancestor = readmPing('HEG4.mping.ref_only.gff', 0)
        mPing_rils     = readtable(args.input)
        mPing_frq      = mping_frequency(args.input)
        mPing_excision = excision(mPing_ancestor, mPing_rils, mPing_frq, 'mping.excision.ref_only.number', 'mping.excision.ref_only.bamcheck')
        preplot(mPing_frq, mPing_ancestor, mPing_excision, 'mping.excision.ref_only.table')
    elif args.mode == 'Shared':
        mPing_ancestor = readmPing('HEG4.mping.shared.gff', 2)
        mPing_rils     = readtable(args.input)
        mPing_frq      = mping_frequency(args.input)
        mPing_excision = excision(mPing_ancestor, mPing_rils, mPing_frq, 'mping.excision.shared.number', 'mping.excision.shared.bamcheck')
        preplot(mPing_frq, mPing_ancestor, mPing_excision, 'mping.excision.shared.table')
'''

#if __name__ == '__main__':
    #main()

