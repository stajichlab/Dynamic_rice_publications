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
python Extract_mPing_varients.py --input rice3k_mPing_target

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines 2 --interval 120 --task 1 --mem 15G --time 100:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

def convert_name(fastafile):
    #264_rice3k.fa.6 Query:mPing_TSD-UNK Sbjct:ERS468130.repeat_Chr5_3221898_3221900_len476 Length:430 Location:(59 - 477) Direction:minus
    s = re.compile(r'Sbjct:(.*)\_len\d+')
    #s = re.compile(r'Sbjct:(\w+)\.repeat_(.*)_len*')
    newfile = re.sub(r'.Nclean.fa', r'.Nclean.name.fa', fastafile)
    ofile = open(newfile, 'w')
    for record in SeqIO.parse(fastafile,"fasta"):
        unit = re.split(r' ', record.description)
        m = s.search(str(unit[2]))
        print unit[2]
        if m:
            print 'match'
            record.description = record.id
            record.id = m.groups(0)[0]
            SeqIO.write(record, ofile, 'fasta')
    ofile.close()
    #for i in old.keys():
    #    if not new.has_key(i):
    #        print '%s missed' %(i)

def summary_mPing_assembly(fasta_files, output):
    fastaid = defaultdict(str)
    s = re.compile(r'Sbjct:(\w+)\.repeat_(.*)_len*')
    qc_seq_num = [0, 0]
    summary_strain = defaultdict(lambda : int())
    summary_locus  = defaultdict(lambda : int())
    print output
    ofile = open('%s.mPing_assembly.summary' %(output), 'w')
    ofile_strain = open('%s.mPing_assembly.strain' %(output), 'w')
    ofile_locus  = open('%s.mPing_assembly.locus' %(output), 'w')
    for fastafile in fasta_files:
        for record in SeqIO.parse(fastafile, "fasta"):
            #print record.id
            #print record.description
            #264_rice3k.fa.6 Query:mPing_TSD-UNK Sbjct:ERS468130.repeat_Chr5_3221898_3221900_len476 Length:430 Location:(59 - 477) Direction:minus
            unit = re.split(r' ', record.description)
            m = s.search(unit[2])
            qc_seq_num[0] += 1
            if m:
                strain = m.groups(0)[0]
                locus  = m.groups(0)[1]
                qc_seq_num[1] += 1
                summary_strain[strain] = 1
                summary_locus[locus]   = 1
    print >> ofile, 'Assembled mPing from %s strains' %(len(summary_strain.keys()))
    print >> ofile, 'Assembled mPing from %s loci' %(len(summary_locus.keys()))
    print >> ofile_strain, '\n'.join(sorted(summary_strain.keys()))
    print >> ofile_locus, '\n'.join(sorted(summary_locus.keys()))
    ofile_strain.close()
    ofile_locus.close()
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

#rice3k_mPing_target/Target_Run_mPing_rice3k.fa.10_2017_03_11_205722/mping_split1/mping_split1.flank_filter-2.0_under
def clean_mPing_element(indir, output, cpu):
    output = os.path.abspath(output)
    final_clean_mPing_flank_fa = []
    final_clean_mPing_N_fa = []
    fa_files = glob.glob('%s/*/mping_split1/mping_split1.flank_filter-*_under' %(indir))
    s = re.compile(r'_(rice3k.fa.\d+)_')
    for fa in sorted(fa_files):
        m = s.search(fa)
        if m:
            prefix = m.groups(0)[0]
            prefix = re.sub(r'.fa', r'', prefix)
            fa = os.path.abspath(fa)
            if not os.path.exists('%s.mPing.clean.Nclean.fa' %(prefix)):
                os.system('ln -s %s %s.mPing.fa' %(fa, prefix))
                os.system('python clean_target_candidate.py --input %s.mPing.fa' %(prefix))
                os.system('python clean_target_candidate_N.py --input %s.mPing.clean.fa' %(prefix))
                convert_name('%s.mPing.clean.Nclean.fa' %(prefix))
            final_clean_mPing_flank_fa.append('%s.mPing.clean.fa' %(prefix))
            final_clean_mPing_N_fa.append('%s.mPing.clean.Nclean.fa' %(prefix))
    #return final_clean_mPing_flank_fa, final_clean_mPing_N_fa

    if not os.path.exists('%s.mPing.clean.Nclean.name.blast.table' %(output)):
        os.system('cat *.mPing.clean.Nclean.name.fa > %s.mPing.clean.Nclean.name.fa' %(output))
        os.system('/opt/linux/centos/7.x/x86_64/pkgs/ncbi-blast/2.2.26/bin/formatdb -i %s.mPing.clean.Nclean.name.fa -p F' %(output))
        os.system('/opt/linux/centos/7.x/x86_64/pkgs/ncbi-blast/2.2.26/bin/blastall -p blastn -a %s -i %s.mPing.clean.Nclean.name.fa -d %s.mPing.clean.Nclean.name.fa -o %s.mPing.clean.Nclean.name.blast' %(cpu, output, output, output))
        os.system('perl ~/BigData/software/bin/blast_parser.pl %s.mPing.clean.Nclean.name.blast > %s.mPing.clean.Nclean.name.blast.table' %(output, output))
    if not os.path.exists('%s.mPing.clean.Nclean.name.blast.table.representive.mPing_var.fa' %(output)):
        os.system('python Blast_classifier.py --input %s.mPing.clean.Nclean.name.blast.table' %(output))
        os.system('cat %s.mPing.clean.Nclean.name.blast.table.representive.fa mPing_variants.fa > %s.mPing.clean.Nclean.name.blast.table.representive.mPing_var.fa' %(output, output))
    if not os.path.exists('%s.mPing.clean.Nclean.name.blast.table.representive.mPing_var.vcf' %(output)):
        os.system('/opt/linux/centos/7.x/x86_64/pkgs/muscle/3.8.425/bin/muscle -in %s.mPing.clean.Nclean.name.blast.table.representive.mPing_var.fa -out %s.mPing.clean.Nclean.name.blast.table.representive.mPing_var.msa' %(output, output))
        os.system('/opt/linux/centos/7.x/x86_64/pkgs/java/jdk1.8.0_45/bin/java -jar /rhome/cjinfeng/BigData/software/jvarkit/dist/msa2vcf.jar -a %s.mPing.clean.Nclean.name.blast.table.representive.mPing_var.msa > %s.mPing.clean.Nclean.name.blast.table.representive.mPing_var.all.vcf' %(output, output))
        os.system('/opt/linux/centos/7.x/x86_64/pkgs/java/jdk1.8.0_45/bin/java -jar /rhome/cjinfeng/BigData/software/jvarkit/dist/msa2vcf.jar %s.mPing.clean.Nclean.name.blast.table.representive.mPing_var.msa > %s.mPing.clean.Nclean.name.blast.table.representive.mPing_var.vcf' %(output, output))
        os.system('python vcf2haplotype_plot.py %s.mPing.clean.Nclean.name.blast.table.representive.mPing_var.vcf' %(output))
        os.system('/rhome/cjinfeng/BigData/software/fasttree/FastTree -noml -nome -nt %s.mPing.clean.Nclean.name.blast.table.representive.mPing_var.msa > %s.mPing.clean.Nclean.name.blast.table.representive.mPing_var.tree' %(output, output))
        #ofile = open('%s.mPing.clean.Nclean.mPing_var.sh' %(output), 'w')
        #print >> ofile, '\n'.join(cmd)
        #ofile.close() 
    return final_clean_mPing_flank_fa, final_clean_mPing_N_fa

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('--cpu')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.cpu:
        args.cpu = 1

    if not args.output:
        args.output = args.input[:-1] if args.input[-1] == '/' else args.input
    final_clean_mPing_flank_fa, final_clean_mPing_N_fa  = clean_mPing_element(args.input, args.output, args.cpu)
    #final_clean_mPing_fa = ['rice3k.6.mPing.clean.Nclean.fa']
    summary_mPing_assembly( final_clean_mPing_flank_fa, '%s.clean_flank' %(args.output))
    summary_mPing_assembly( final_clean_mPing_N_fa, '%s.clean_N' %(args.output))


if __name__ == '__main__':
    main()

