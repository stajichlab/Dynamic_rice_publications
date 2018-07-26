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
python Assign_Individual_Pop.py --Qtable core_v0.7.pruneddata3.8.Q.strains_order.list --colortable core_v0.7.pruneddata3.8.Q.color.txt
Assign individual to predefined subgroups using ADMIXTURE results.
1. Add the proportion of all cluster from each predefined subgroup
2. Assign an individual to a subgroup if the proportion is larger than 70%
3. If an individual is not assign to any subgroup and assign it to admixture if it has two and more subgroup larger than 20%

#core_v0.7.pruneddata3.8.Q.strains_order.list has proportion of strain origin from each columen (cluster)
IRIS_313-10628	0.000010 0.000010 0.000010 0.999930 0.000010 0.000010 0.000010 0.000010
IRIS_313-10640	0.000010 0.000010 0.000010 0.999930 0.000010 0.000010 0.000010 0.000010

#core_v0.7.pruneddata3.8.Q.color.txt has a predefined population for each columen (clustser)
0	Indica	green
1	Temperate jap	blue
2	Tropical jap	cornflowerblue
3	Indica	green
4	Indica	green
5	Basmati/sadri	darkorchid
6	Aus/boro	chocolate
7	Indica	green


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

#IRIS_313-10628	0.000010 0.000010 0.000010 0.999930 0.000010 0.000010 0.000010 0.000010
#IRIS_313-10640	0.000010 0.000010 0.000010 0.999930 0.000010 0.000010 0.000010 0.000010
def readQtable(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t| ',line)
                data[unit[0]] = unit[1:]
                #print '%s\t%s' %(unit[0], data[unit[0]])
    return data

#0	Indica	green
#1	Temperate jap	blue
#2	Tropical jap	cornflowerblue
#3	Indica	green
#4	Indica	green
#5	Basmati/sadri	darkorchid
#6	Aus/boro	chocolate
#7	Indica	green
def readcolortable(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[int(unit[0])] = unit[1]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--Qtable')
    parser.add_argument('--colortable')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.Qtable) > 0 and len(args.colortable) > 0
    except:
        usage()
        sys.exit(2)
    
    pop   = readQtable(args.Qtable)
    color = readcolortable(args.colortable) 
    
    subpop = list(set(color.values()))
    print 'Strain\tSubgroup\t%s' %('\t'.join(sorted(subpop)))
    for strain in sorted(pop.keys()):
        subpop_meta = defaultdict(lambda : float())
        for i in range(5):
            #print '%s\t%s' %(i, pop[strain][i])
            subpop_meta[color[i]] += float(pop[strain][i])
            #print '%s\t%s' %(color[i], subpop_meta[color[i]])
        line  = []
        max_rate = 0.0
        max_pop  = ''
        mix_pop  = 0
        for p in sorted(subpop):
            line.append(str(subpop_meta[p]))
            if subpop_meta[p] > max_rate:
                max_rate = subpop_meta[p]
                max_pop  = p
            if subpop_meta[p] > 0.3:
                mix_pop += 1
        #top subgroup contribute to more than 80% of the genome
        if max_rate >= 0.8:
            print '%s\t%s\t%s' %(strain, max_pop, '\t'.join(line))
        #have two or more subgroups admixture
        elif mix_pop >= 2:
            print '%s\t%s\t%s' %(strain, 'admixture', '\t'.join(line))
        #other admixture
        else:
            print '%s\t%s\t%s' %(strain, 'admixture', '\t'.join(line))
        #print '%s\t%s' %(strain, '\t'.join(line))
        

if __name__ == '__main__':
    main()

