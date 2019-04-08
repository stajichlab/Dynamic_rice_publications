#!/opt/Python/2.7.3/bin/python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import re
import argparse
import sys

def usage():
    message='Python ChrName.py --input INPUT --output OUTPUT'
    print message

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.output) > 1
    except:
        usage()
        sys.exit(2)
    
    scaf    = defaultdict(lambda : str())
    scaflen = defaultdict(lambda : int())
    ofile = open(args.output,"w")
    ofile1 = open('%s.id.table' %(args.output), "w")
    for record in SeqIO.parse(args.input,"fasta"):
        if str(record.id) in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9']:
            SeqIO.write(record, ofile, "fasta")
        else:
            scaf[str(record.id)]    = record.seq
            scaflen[str(record.id)] = len(str(record.seq))
    rank = 9
    for scafid in sorted(scaflen, key=scaflen.get, reverse=True):
        if len(scaf[scafid]) < 1000:
            continue
        rank += 1
        newrecord = SeqRecord(scaf[scafid], id='scaffold%s' %(rank), description="")
        print >> ofile1, '%s\t%s\t%s' %(scafid, 'scaffold%s' %(rank), len(scaf[scafid]))
        SeqIO.write(newrecord, ofile, "fasta")
    ofile.close()
    ofile1.close()
