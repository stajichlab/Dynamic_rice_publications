#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir
import networkx as nx

def usage():
    test="name"
    message='''
python Blast_classifier.py --input rice3k.9.mPing.clean.Nclean.mPing_var.blast.table

Use blast results to classfy sequence and generate fasta of representive sequence

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines 2 --interval 120 --task 1 --mem 15G --time 100:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def sub_fasta_seq(fastafile, seq_list):
    fasta_s = defaultdict(lambda : str())
    for record in SeqIO.parse(fastafile,"fasta"):
        if seq_list.has_key(record.id):
            fasta_s[record.id] = record.seq
    return fasta_s

#1	B001	China	Temperate japonica	ERS470219
def read_download_table(infile):
    data = defaultdict(lambda : str())
    pop  = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                #strain = re.sub(r'IRIS', r'IRIS_', unit[1])
                name = unit[3]
                if unit[3] == 'Aromatic (basmati/sandri type':
                    name = 'Basmati/sadri'
                if unit[3] == 'Aus':
                    name = 'Aus/boro'
                if unit[3] == 'japonica':
                    name = 'Japonica'
                data[unit[4]] = name
                pop[name]  = 1
    return data, pop

#summary population distribution of given list of strains
def population_summary_list(threek_anno, strain_list):
    data = defaultdict(lambda : int())
    for strain in strain_list:
        if threek_anno.has_key(strain):
            data[threek_anno[strain]] += 1
    return data
#
#ERS467823.repeat_Chr2_25755510_25755512
def summary_strain_and_loci(seq_ids):
    strain_list = defaultdict(lambda : int())
    locus_list  = defaultdict(lambda : int())
    for seq_id in seq_ids:   
        unit = re.split(r'\.', seq_id)
        strain_list[unit[0]] += 1
        locus_list[unit[1]]  += 1
    return strain_list, locus_list 

#Query_id        Query_length    Query_start     Query_end       Subject_id      Subject_length  Subject_start   Subject_end     Identity        Positive        Gap     Align_length    Score
#1_rice3k.fa.9   450     1       450     mPingD  450     1       450     1       --      0       450     892     0.0     Query:mPing_TSD-UNKSbjct:ERS468296
def blast_classifier(infile):
    #data = defaultdict(str)
    threek_anno, threek_pop = read_download_table('rice_line_3000.download.list')
    blast_g = nx.Graph()
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Query_id'): 
                unit = re.split(r'\t',line)
                #not self hit
                if not unit[0] == unit[4]:
                    #query and target are 100% aligned without gap and mismatch
                    if int(unit[11]) == int(unit[1]) and int(unit[11]) == int(unit[5]) and int(unit[10]) == 0 and float(unit[8]) == 1:
                        blast_g.add_edge(unit[0], unit[4])
    #write representive id
    blast_g_sub = list(nx.connected_component_subgraphs(blast_g))
    rank = 0
    reprensentive = defaultdict(lambda : str())
    ofile_class = open('%s.class' %(infile), 'w')
    print >> ofile_class, 'Node\tNumberOfNode\t%s\tNumberOfStrains\tNumberOfLoci\tStrainDetails\tLocusDetails\tRepresentiveLocus\tListOfNode' %('\t'.join(sorted(threek_pop.keys())))
    ofile_representive = open('%s.representive.fa' %(infile), 'w')
    for subgraph in blast_g_sub:
        rank += 1
        node = 'node%s' %(rank)
        subgraph_strain, subgraph_loci = summary_strain_and_loci(subgraph.nodes())
        subgraph_strain_distr = population_summary_list(threek_anno, subgraph_strain)
        ####
        strain_string = []
        for strain in sorted(subgraph_strain.keys()):
            strain_string.append('%s:%s' %(strain, str(subgraph_strain[strain])))
        locus_string  = []
        for locus in sorted(subgraph_loci.keys()):
            locus_string.append('%s:%s' %(locus, str(subgraph_loci[locus])))
        pop_string    = []
        for pop in sorted(threek_pop.keys()):
            count = subgraph_strain_distr[pop] if subgraph_strain_distr.has_key(pop) else 0
            pop_string.append(str(count)) 
        print >> ofile_class, 'node%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(str(rank), str(subgraph.number_of_nodes()), '\t'.join(pop_string), len(subgraph_strain.keys()), len(subgraph_loci.keys()), ','.join(strain_string), ','.join(locus_string),  sorted(subgraph.nodes())[0], ','.join(subgraph.nodes()))
        #output only graph have more than 10 nodes into fasta 
        if subgraph.number_of_nodes() >= 10:
            rep_id = sorted(subgraph.nodes())[0]
            reprensentive[rep_id] = node
    #write representive fasta
    fastafile = re.sub(r'.blast.table', r'.fa', infile)
    rep_seq = sub_fasta_seq(fastafile, reprensentive)
    for seq_id in rep_seq.keys():
        rep_record = SeqRecord(rep_seq[seq_id], id=seq_id, name=seq_id, description=reprensentive[seq_id])
        SeqIO.write(rep_record, ofile_representive, 'fasta') 
    ofile_class.close()
    ofile_representive.close()


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

    #3k_anno = read_download_table('rice_line_3000.download.list')
    blast_classifier(args.input)

if __name__ == '__main__':
    main()

