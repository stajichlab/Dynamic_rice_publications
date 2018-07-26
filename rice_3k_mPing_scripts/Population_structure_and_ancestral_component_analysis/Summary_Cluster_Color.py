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
python Summary_Cluster.py --strain_list core_v0.7.pruneddata.structure.nosex --meanQ_list testoutput_simple_4.4.meanQ

The order of strains in two input files should be same

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines 2 --interval 120 --task 1 --mem 15G --time 100:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#B001    6       0       6       379     378     0       B001    blue    Heibiao|B001|Temp       Heibiao China   Temperate jap
def readtable(infile):
    data = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int())))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Taxa'): 
                unit = re.split(r'\t',line)
                if unit[0] is not '':
                    #G-type Ping
                    if int(unit[1]) > 1 and int(unit[2]) > 1 and int(unit[3]) == 0:
                        data['ping']['G'][unit[12]] += 1
                    #A-type Ping
                    elif int(unit[1]) > 1 and int(unit[2]) == 0 and int(unit[3]) > 1:
                        data['ping']['A'][unit[12]] += 1
                    #GA-type Ping
                    elif int(unit[1]) > 1 and int(unit[2]) > 1 and int(unit[3]) > 1:
                        data['ping']['GA'][unit[12]] += 1
                    #G-tu[e mPing
                    if int(unit[4]) > 1 and int(unit[5]) > 1 and int(unit[6]) == 0:
                        data['mping']['G'][unit[12]] += 1
                    #A-type mPing
                    elif int(unit[4]) > 1 and int(unit[5]) == 0 and int(unit[6]) > 1:
                        data['mping']['A'][unit[12]] += 1
                    #GA-type mPing
                    elif int(unit[4]) > 1 and int(unit[5]) > 1 and int(unit[6]) > 1:
                        data['mping']['GA'][unit[12]] += 1
    for e in ['ping', 'mping']:
        for pop in sorted(data[e]['G'].keys()):
            print '%s G-type %s strains: %s' %(e, pop, data[e]['G'][pop])
        for pop in sorted(data[e]['A'].keys()):
            print '%s A-type %s strains: %s' %(e, pop, data[e]['A'][pop])
        for pop in sorted(data[e]['GA'].keys()):
            print '%s GA-type %s strains: %s' %(e, pop, data[e]['GA'][pop])
    return data

#Taxa	Color	Label	Name	Origin	Group	mPing	Ping	Pong
#B001	blue	Heibiao|B001|Temp	Heibiao	China	Temperate jap	71	1	8
def read_anno(infile):
    data = defaultdict(lambda : str())
    pop  = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Taxa'): 
                unit = re.split(r'\t',line)
                data[unit[0]] = unit[5]
                pop[unit[5]]  = 1
    return data, pop

def read_table(infile):
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Taxa'): 
                unit = re.split(r'\s+', line)
                data.append(unit)
    return data

def sum_anno(anno):
    pop_sum = defaultdict(lambda : int())
    for s in anno.keys():
        pop_sum[anno[s]] += 1
    return pop_sum

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--strain_list')
    parser.add_argument('--meanQ_list')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.strain_list) > 0 and len(args.meanQ_list) > 0 
    except:
        usage()
        sys.exit(2)

    meta = defaultdict(lambda : defaultdict(lambda : int()))
    tk_anno, pop = read_anno('rice_line_ALL_3000.anno.list')
    strains      = read_table(args.strain_list)
    cluster      = read_table(args.meanQ_list)
    for i in range(len(strains)):
        #print i
        s = strains[i][0]
        c = cluster[i]
        #print c
        c_max = max(map(float, c))
        if float(c_max) > 0.8:
            c_index = c.index('%0.6f' %(c_max))
            meta[c_index][tk_anno[s]] += 1
    #indica, aus/boro, aro, temperate, tropical
    color = {'Indica':'green', 'Aus/boro':'chocolate', 'Basmati/sadri':'darkorchid', 'Temperate jap':'blue', 'Tropical jap':'cornflowerblue'}
    #print '%s\t%s' %('Cluster', '\t'.join(sorted(pop.keys()))) 
    ofile = open('%s.color.txt' %(args.meanQ_list), 'w')
    for c_index in sorted(meta.keys()):
        pop_distr = []
        n_max     = 0
        p_max     = ''
        for p in sorted(pop.keys()):
            #pop_distr.append(meta[c_index][p])
            if meta[c_index][p] > n_max:
                n_max = meta[c_index][p]
                p_max = p
        print >> ofile, '%s\t%s\t%s' %(str(c_index), p_max, color[p_max])
    ofile.close()

if __name__ == '__main__':
    main()

