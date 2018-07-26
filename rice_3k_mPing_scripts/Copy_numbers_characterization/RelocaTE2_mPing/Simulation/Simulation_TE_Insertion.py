#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random

def usage():
    test="name"
    message='''
python Simulation_TE_Insertion.py --repeat mping.fa --genome MSU_r7.fa --te mPing  --chr Chr4 --number 200

Given the TE and genome sequence, simulate random insertion site in the genome or individual chromosome

    '''
    print message

##
def write_fasta(ref, chrs, filename):
    ofile = open(filename, "w")
    for chrn in sorted(ref.keys()):
        if chrs == 'ALL':
            seq = Seq(ref[chrn])
            newrecord = SeqRecord(seq, id=chrn,description="")
            SeqIO.write(newrecord, ofile, "fasta")
        elif chrs == chrn:
            seq = Seq(ref[chrn])
            newrecord = SeqRecord(seq, id=chrn,description="")
            SeqIO.write(newrecord, ofile, "fasta")
    ofile.close()

##repeat fasta
def fasta_te(fastafile, te):
    fastaid = defaultdict(str)
    fastatsd = defaultdict(list)
    ids = []
    for record in SeqIO.parse(fastafile, "fasta"):
        if te == str(record.id):
            fastaid[record.id] = str(record.seq)
            ids.append(str(record.id))
    #get tsd
    s = re.compile(r'>(.\S+)\s+TSD\=(\w+)')
    disc = SeqIO.index(fastafile, "fasta")
    for d in ids:
        lines = re.split('\n', disc.get_raw(d).decode())
        head  = lines[0]
        m = s.search(head)
        repid = ''
        tsd   = ''
        if m:
            repid = m.groups(0)[0]
            tsd   = m.groups(0)[1]
            fastatsd[repid] = [fastaid[repid], tsd]
            #print repid, tsd
    return fastatsd

##reference fasta
def fasta_ref(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = str(record.seq)
    return fastaid

##random choose pos from chromosome
def random_pos(chrn, ref, tsd):
    seq_len = len(ref[chrn])
    rand_pos= random.randint(1,seq_len)
    strand  = 1 if random.randint(1,10) > 5 else 0
    #tsd
    tsdstart = rand_pos - len(tsd)
    chrseq   = ref[chrn]
    tsdseq   = chrseq[tsdstart:rand_pos]
    #tsdseq   = chrseq[tsdstart:rand_pos] if strand == 1 else str(Seq(chrseq[tsdstart:rand_pos]).reverse_complement())
    return [chrn, rand_pos, rand_pos+1, strand, tsdseq]

##except for the for insertion, all the following insertion need to add inserted sequence to their position
def update_pos(data, repseq, tsd):
    ##tsd only add one copy of itself
    add_len = len(repseq) + 1*len(tsd)
    add_num = defaultdict(int)
    for i in range(0, len(data)):
        chrn  = data[i][0]
        start = data[i][1]
        end   = data[i][2]
        #print chrn, start, end
        if add_num.has_key(chrn):
            data[i][1] = start + add_num[chrn]*add_len
            data[i][2] = end   + add_num[chrn]*add_len
            add_num[chrn] += 1
        else:
            add_num[chrn] += 1
        #print data[i][0], data[i][1], data[i][2]
    return data

##insertion element into genome
def insert_element(pos, repseq, tsd, ref):
    chrn = pos[0]
    chrseq = ref[chrn]
    half1 = chrseq[:pos[1]]
    half2 = chrseq[pos[1]:]
    ##we choose sequence at target site as tsd, not use tsd provided
    tsdstart = pos[1] - len(tsd)
    tsdseq   = chrseq[tsdstart:pos[1]]
    newseq   = ''
    if pos[3] == 1:
        newseq = half1 + repseq + tsdseq + half2
    else:
        repseq_seq = Seq(repseq)
        repseq_rec = repseq_seq.reverse_complement()
        newseq = half1 + str(repseq_rec) + tsdseq + half2
    ref[chrn] = newseq

##write simulated insertion into gff format
def writegff(data, te, gff_out):
    ofile = open(gff_out, 'w') 
    for pos in data:
        chrn   = pos[0]
        start  = pos[1]
        end    = pos[2]
        strand = '+' if pos[3] == 1 else '-'
        tsd    = pos[4]
        te_id  = '%s.%s.%s' %(te, chrn, str(start))
        print >> ofile, '%s\tMSU7\t%s\t%s\t%s\t.\t%s\t.\tID=%s;TSD=%s;' %(chrn, te, str(start), str(end), strand, te_id, tsd)

##main function of simulation
def simulate(ref, repseq, tsd, te, chrs, number, prefix):
    data = []
    #print repseq
    for n in range(int(number)):
        #print 'Rank %s' %(n)
        if chrs == 'ALL':
            c = random.randint(1,12)
            chrn = 'Chr%s' %(c)
            ins  = random_pos(chrn, ref, tsd)
            data.append(ins)
            #print ins[0], ins[1], ins[2]
        else:
            ins  = random_pos(chrs, ref, tsd)
            data.append(ins)
 
    data = sorted(data, key = lambda x: (x[0], x[1]))
    gff_out = '%s.gff' %(prefix)
    writegff(data, te, gff_out)
    datac = data
    datac = update_pos(datac, repseq, tsd)
    for pos in datac:
        #print pos[0], pos[1], pos[2]
        insert_element(pos, repseq, tsd, ref)
    fasta_out = '%s.fasta' %(prefix)
    write_fasta(ref, chrs, fasta_out)    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--repeat')
    parser.add_argument('-g', '--genome')
    parser.add_argument('-c', '--chr')
    parser.add_argument('-t', '--te')
    parser.add_argument('-n', '--number')
    parser.add_argument('-p', '--prefix')
    parser.add_argument('-o', '--outdir')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.repeat) > 0 or len(args.genome) > 0
    except:
        usage()
        sys.exit(2)

    if not args.chr:
        args.chr = "ALL"

    if not args.number:
        args.number = 200

    if not args.te:
        args.te = 'mPing'
   
    if not args.outdir:
        args.outdir = 'MSU7.%s.%s' %(args.chr, args.te)
    
    if not args.prefix:
        args.prefix = '%s.1' %(args.outdir)
    
    ref = fasta_ref(args.genome)
    element  = fasta_te(args.repeat, args.te)
    #print args.te,element[args.te][0], element[args.te][1]
    simulate(ref, element[args.te][0], element[args.te][1], args.te, args.chr, args.number, args.prefix)
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    cmd = 'mv %s.* %s' %(args.prefix, args.outdir)
    os.system(cmd)

if __name__ == '__main__':
    main()

