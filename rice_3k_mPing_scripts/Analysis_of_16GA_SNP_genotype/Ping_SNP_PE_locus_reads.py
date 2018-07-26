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
import subprocess
import pysam

def usage():
    test="name"
    message='''
python Ping_SNP.py --input fq_RelocaTE2

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines %s --interval 120 --task 1 --mem 15G --time 10:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

def pop_sum(strains, strains3k_data, acc2name, tp):
    data_sum = defaultdict(lambda : int())
    for s in strains.keys():
        name = acc2name[s]
        pop  = strains3k_data[name]
        data_sum[pop] += 1
        print '%s\t%s\t%s\t%s\t%s' %(tp, s, name, pop, data_sum[pop])
    return data_sum

#Taxa	Color	Label	Name	Origin	Group	mPing	Ping	Pong
#B001	blue	Heibiao|B001|Temp	Heibiao	China	Temperate jap	71	1	8
def sum_3k_pop(infile):
    data = defaultdict(lambda : str())
    pop  = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Taxa'): 
                unit = re.split(r'\t',line)
                data[unit[0]] = unit[5]
                pop[unit[5]]  += 1
    return data, pop

#1	B001	China	Temperate japonica	ERS470219
def read_download_table(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                strain = re.sub(r'IRIS', r'IRIS_', unit[1])
                data[unit[4]] = strain
               
    return data

#Chr9    ERS469546       RelocaTE_i      15710721        15710723        .       +       .       ID=repeat_Chr9_15710721_15710723
def read_intact_gff(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                rep  = 'repeat_%s_%s_%s' %(unit[0], unit[3], unit[4])
                acc  = unit[1]
                print "reading: %s" %(rep)
                data[rep] = 1
                data[acc] = 1
    return data

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr3	ERS469922	RelocaTE_i	33213691	33213693
def readgff(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[1]] = 1
    return data

def count_nucleotides(dna, nucleotide):
    return dna.lower().count(nucleotide.lower())

#input bam file of SNP splitted
#return dict with read cover the SNP
def snp_reads_from_bam(bam, snp):
    fsam = pysam.AlignmentFile(bam, 'rb')
    data = defaultdict(lambda : int())
    for record in fsam.fetch(until_eof = True):
        if not record.is_unmapped:
            name   = record.query_name
            start  = int(record.reference_start) + 1
            end    = int(record.reference_end) + 1
            if start < 16 and end > 16:
                data[name] = 1
    return data

#overlap SNP reads and Ping
#input ping_reads, G_reads and A_reads
#return dict with SNP->reads->1 
def overlap_reads(ping_reads, G_reads, A_reads):
    ping_reads_class = defaultdict(lambda : defaultdict(lambda : int()))
    for ping_read in sorted(ping_reads.keys()):
        if G_reads.has_key(ping_read):
            ping_reads_class['G'][ping_read] = 1
        elif A_reads.has_key(ping_read):
            ping_reads_class['A'][ping_read] = 1
        else:
            ping_reads_class['UN'][ping_read] = 1
    return ping_reads_class 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-g', '--gff')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.gff:
        args.gff = 'Rice3k_3000_RelocaTEi_Ping.CombinedGFF.ALL.high_conf.gff'
    args.gff = os.path.abspath(args.gff)    

    home_dir = os.path.split(os.path.realpath(__file__))[0]
    ping = os.path.abspath('ping.fa')
    relocate2_dirs = glob.glob("%s/*_RelocaTEi" %(args.input))
    ofile = open('%s.run_ping_SNP.sh' %(args.input), 'w')
    number_of_unit     = 0
    number_of_job_per  = 0
    ping_strain = readgff(args.gff)
    for strain_dir in sorted(relocate2_dirs):
        #print >> ofile, strain_dir
        strain = re.sub(r'_RelocaTEi', r'', os.path.split(strain_dir)[1])
        #strain = re.sub(r'_RelocaTE2', r'', os.path.split(strain_dir)[1])
        if not ping_strain.has_key(strain):
            continue
        strain_dir = os.path.abspath(strain_dir)
        strain     = os.path.abspath(strain)
        cmd = []
        #../fq_RelocaTE2_Ping/nivara_IRGC105327_RelocaTE2/repeat/results/Chr1.repeat.reads.list
        assembly_dir = '%s_mPing_assembly' %(strain)
        prefix       = '%s/%s' %(assembly_dir, os.path.split(strain)[1])
        cmd.append('mkdir %s' %(assembly_dir))
        cmd.append('cat %s/repeat/results/*.repeat.reads.list > %s.repeat.reads.list' %(strain_dir, prefix))
        cmd.append('cat %s/repeat/te_containing_fq/*_1.te_repeat.ContainingReads.fq > %s_1.te_repeat.ContainingReads.fq' %(strain_dir, prefix))
        cmd.append('cat %s/repeat/te_containing_fq/*_2.te_repeat.ContainingReads.fq > %s_2.te_repeat.ContainingReads.fq' %(strain_dir, prefix))
        cmd.append('python %s/Change_Fastq_Header.py --fastq %s_1.te_repeat.ContainingReads.fq --header read1' %(home_dir, prefix))
        cmd.append('python %s/Change_Fastq_Header.py --fastq %s_2.te_repeat.ContainingReads.fq --header read2' %(home_dir, prefix))
        #cmd.append('cat %s/repeat/flanking_seq/*_1.te_repeat.flankingReads.fq.matched > %s_1.te_repeat.flankingReads.fq.matched' %(strain_dir, prefix))
        #cmd.append('cat %s/repeat/flanking_seq/*_2.te_repeat.flankingReads.fq.matched > %s_2.te_repeat.flankingReads.fq.matched' %(strain_dir, prefix))
        #cmd.append('cat %s/repeat/flanking_seq/*_1.te_repeat.flankingReads.unPaired.fq > %s_2.te_repeat.flankingReads.unPaired.fq' %(strain_dir, prefix))
        cmd.append('cat %s_1.te_repeat.ContainingReads.name.fq %s_2.te_repeat.ContainingReads.name.fq > %s.te_repeat.ContainingReads.fq' %(prefix, prefix, prefix))
        cmd.append('python %s/Ping_locus_reads_list.py --input %s.repeat.reads.list --gff %s' %(home_dir, prefix, args.gff))
        cmd.append('python %s/Get_List_Fastq_runner.py --input %s --ping %s' %(home_dir, prefix, ping))       
        #cmd.append('rm %s.*.list %s.*.fq %s_*.fq' %(prefix, prefix, prefix))
        #
        number_of_job_per = len(cmd)
        number_of_unit   += 1
        ##assembly
        #cmd.append('cat %s_1.te_repeat.ContainingReads.fa %s_2.te_repeat.ContainingReads.fa > %s.te_reads_all.fa' %(strain, strain, strain))
        #cmd.append('python ~/BigData/software/bin/fasta2fastq.py %s.te_reads_all.fa %s.te_reads_all.fq' %(strain, strain))
        #cmd.append('python ~/BigData/software/bin/Split_Fastq2PE.py --input %s.te_reads_all.fq' %(strain))
        #cmd.append('/rhome/cjinfeng/BigData/software/Velvet/velvet/velveth %s.assembly 31 -shortPaired -fastq -separate %s.te_reads_all_1.fq %s.te_reads_all_2.fq' %(strain, strain, strain))
        #cmd.append('/rhome/cjinfeng/BigData/software/Velvet/velvet/velvetg %s.assembly -ins_length 400 -exp_cov 50 -min_contig_lgth 200 -scaffolding yes' %(strain))
        ##
        #cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/bwa/0.7.12/bin/bwa mem -p %s %s.te_reads_smart.fq > %s.sam' %(ping, strain, strain))
        #cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools view -bS -o %s.raw.bam %s.sam' %(strain, strain)) 
        #cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools sort -m 1000000000 %s.raw.bam %s' %(strain, strain))
        #cmd.append("/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools view -h %s.bam | perl -lane 'if($F[11] =~ /^NM:i:(\d+)$/){print if $1<=2}else{print}'| /opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools view -bS - -o %s.NM2.bam" %(strain, strain))
        #cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools index %s.bam' %(strain))
        #cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools index %s.NM2.bam' %(strain))
        #cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/samtools/0.1.19/bin/samtools mpileup %s.NM2.bam > %s.NM2.mpileup' %(strain, strain))
        #cmd.append('rm %s.raw.bam %s.sam %s_1.te_repeat.5_3_prime.fa %s_2.te_repeat.5_3_prime.fa %s_1.te_repeat.ContainingReads.fq %s_2.te_repeat.ContainingReads.fq %s_1.te_repeat.ContainingReads.fa %s_2.te_repeat.ContainingReads.fa %s_1.te_repeat.ContainingReads_middle.fa %s_2.te_repeat.ContainingReads_middle.fa' %(strain, strain, strain, strain, strain, strain, strain, strain, strain, strain))
        print >> ofile, '\n'.join(cmd)
        #cat fq_RelocaTE2/rufipogon_W1715_RelocaTE2/repeat/te_only_read_portions_fa/*.five_prime.fa > rufipogon_W1715_1.te_repeat.five_prime.fa
    ofile.close()

    jobs_per_run = number_of_job_per
    if number_of_unit > 30:
        jobs_per_run = number_of_job_per*(int(number_of_unit/30))
    print 'Total number of strains: %s' %(number_of_unit)
    print 'Number of job per strain: %s' %(number_of_job_per)
    print 'Number of job per run by qsub: %s' %(jobs_per_run)
    #runjob('%s.run_ping_SNP.sh' %(args.input), jobs_per_run)

    intact_ping = read_intact_gff("Rice3k_3000_RelocaTEi_Ping.CombinedGFF.merge.intact.gff")

    strain3k_data, strain3k_sum    = sum_3k_pop('rice_line_ALL_3000.anno.list')
    acc2name         = read_download_table('rice_line_3000.download.list')

    ofile = open('%s.Ping_Locus_16th_SNP.summary' %(args.input), 'w')
    print >> ofile, 'Strain\tPing_loci\t#G\t#T\t#A\t#C'
    #ERS467761_mPing_assembly/ERS467761.repeat_Chr3_14307409_14307411.NM2.mpileup
    mpileup_files = glob.glob('*_mPing_assembly/*.NM2.mpileup')
    G_strains = defaultdict(lambda : int())
    A_strains = defaultdict(lambda : int())
    for f in sorted(mpileup_files):
        cmd = "cat %s | awk '{if($2==16){print $5}}'" %(f) 
        sequence = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read().rstrip()
        g = count_nucleotides(sequence, 'G')
        t = count_nucleotides(sequence, 'T')
        a = count_nucleotides(sequence, 'A')
        c = count_nucleotides(sequence, 'C')
        strain = re.split(r'\.', os.path.split(f)[1])[0]
        locus  = re.sub(r'.NM2.mpileup', r'', re.split(r'\.', os.path.split(f)[1])[1])
        print locus
        if intact_ping.has_key(locus) and intact_ping.has_key(strain):
           
            #print f, os.path.split(f)[1], strain, locus
            print >> ofile, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(strain, acc2name[strain], strain3k_data[acc2name[strain]], locus, g, t, a, c)
            if g >= 2:
                G_strains[strain] = 1
            if a >= 2:
                A_strains[strain] = 1
    ofile.close()

    #summary
    intact_strains   = readgff('Rice3k_3000_RelocaTEi_Ping.CombinedGFF.merge.intact.gff')
    strain3k_data, strain3k_sum    = sum_3k_pop('rice_line_ALL_3000.anno.list')
    acc2name         = read_download_table('rice_line_3000.download.list')
    ping_strains_sum = pop_sum(intact_strains, strain3k_data, acc2name, '3000_strains')
    G_strains_sum    = pop_sum(G_strains, strain3k_data, acc2name, 'G_strains')
    A_strains_sum    = pop_sum(A_strains, strain3k_data, acc2name, 'A_strains') 
  
    ofile = open('%s.Ping_Locus_16th_SNP.pop.summary' %(args.input), 'w')
    print >> ofile, 'Population\t3000\tPing\tPing_G\tPing_A'
    for p in sorted(strain3k_sum.keys()):
        ping_strains_temp = ping_strains_sum[p] if ping_strains_sum.has_key(p) else 'NA'
        G_strains_temp = G_strains_sum[p] if G_strains_sum.has_key(p) else '0'
        A_strains_temp = A_strains_sum[p] if A_strains_sum.has_key(p) else '0'
        print >> ofile, '%s\t%s\t%s\t%s\t%s' %(p, strain3k_sum[p], ping_strains_temp, G_strains_temp, A_strains_temp)
    ofile.close()

'''
    #summary 16th SNP read count
    ofile = open('run_ping_SNP.16th_SNP.summary', 'w')
    mpileup_files = glob.glob('*.mpileup')
    for f in sorted(mpileup_files):
        cmd = "cat %s | awk '{if($2==16){print $5}}'" %(f) 
        sequence = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read().rstrip()
        g = count_nucleotides(sequence, 'G')
        t = count_nucleotides(sequence, 'T')
        a = count_nucleotides(sequence, 'A')
        c = count_nucleotides(sequence, 'C')
        strain = re.sub(r'.mpileup', r'', f)
        print >> ofile, '%s\t%s\t%s\t%s\t%s' %(strain, g, t, a, c)
    ofile.close()

    #split 16th SNP reads
    #os.system('python splitSNP.py ERS470370.NM2.bam ERS470370.NM2.SNPreads ping:16:A:G')

    #summary 16th SNP paired-end reads: ping or mPing?
    ofile_log = open('%s.Ping_SNP_reads.log' %(args.input), 'w')
    ofile_sum = open('%s.Ping_SNP_reads.sum.txt' %(args.input), 'w')
    print >> ofile_sum, 'Taxa\tTotal_Ping_reads\tPing_reads_16thSNP_G\ttPing_reads_16thSNP_A\ttTotal_mPing_reads\tmPing_reads_16thSNP_G\ttmPing_reads_16thSNP_A'
    bam_files = glob.glob('*.NM2.bam')
    for bam in sorted(bam_files):
        print >> ofile_log, 'bam file: %s' %(bam)
        strain = re.sub(r'.NM2.bam', r'', bam)
        #split 16th SNP reads
        prefix = re.sub(r'.bam', r'.SNPreads', bam)
        os.system('python splitSNP.py %s %s ping:16:A:G' %(bam, prefix))
        #
        fsam = pysam.AlignmentFile(bam, 'rb')
        data = defaultdict(lambda : defaultdict(lambda : list()))
        bam_paired_reads_class = '%s.paired_reads_class.txt' %(bam)
        ofile_class = open(bam_paired_reads_class, 'w')
        for record in fsam.fetch(until_eof = True):
            if not record.is_unmapped:
                name   = record.query_name
                start  = int(record.reference_start) + 1
                #end    = int(start) + int(length) - 1 #should not allowed for indel or softclip
                end    = int(record.reference_end) + 1
                read1  = 1 if record.is_read1 else 2
                pair_mapped = 1 if (record.is_paired and not record.mate_is_unmapped) else 0
                #print 'alignment: %s, %s, %s' %(name, start, end)
                if pair_mapped:
                    data[name][read1] = [str(start), str(end)]
        #process each paired-end that mapped on ping element
        summary = defaultdict(lambda : int())
        ping_reads = defaultdict(lambda : int())
        mping_reads = defaultdict(lambda : int())
        for read_pair in sorted(data.keys()):
            #some reads with their pairs filtered out because of too much mismatch, remove these here
            if len(data[read_pair][1]) == 2 and len(data[read_pair][2]) == 2:
                summary['paired'] += 1
                #print '%s\t%s\t%s' %(read_pair, '\t'.join(data[read_pair][1]), '\t'.join(data[read_pair][2]))
                #read1 covers 16th SNP
                if int(data[read_pair][1][0]) < 16 and int(data[read_pair][1][1]) > 16:
                    #read2 entend into Ping
                    if int(data[read_pair][2][1]) >= 300 and int(data[read_pair][2][1]) <= 1000:
                        summary['Ping_reads'] += 1
                        ping_reads[read_pair] = 1
                        print >> ofile_class, '%s\t%s\t%s\tPing_reads' %(read_pair, '\t'.join(data[read_pair][1]), '\t'.join(data[read_pair][2]))                        
                    elif int(data[read_pair][2][1]) < 300:
                        summary['Short_insert'] += 1
                        print >> ofile_class, '%s\t%s\t%s\tShort_insert' %(read_pair, '\t'.join(data[read_pair][1]), '\t'.join(data[read_pair][2]))
                    elif int(data[read_pair][2][1]) > 1000:
                        summary['mPing_reads'] += 1
                        mping_reads[read_pair] = 1
                        print >> ofile_class, '%s\t%s\t%s\tmPing_reads' %(read_pair, '\t'.join(data[read_pair][1]), '\t'.join(data[read_pair][2]))
                    else:
                        summary['error1'] += 1
                        print >> ofile_class, '%s\t%s\t%s\terror1' %(read_pair, '\t'.join(data[read_pair][1]), '\t'.join(data[read_pair][2]))
                #read2 covers 16th SNP
                elif int(data[read_pair][2][0]) < 16 and int(data[read_pair][2][1]) > 16:
                    #read1 extend into Ping
                    if int(data[read_pair][1][1]) >= 300 and int(data[read_pair][1][1]) <= 1000:
                        summary['Ping_reads'] += 1
                        ping_reads[read_pair] = 1
                        print >> ofile_class, '%s\t%s\t%s\tPing_reads' %(read_pair, '\t'.join(data[read_pair][1]), '\t'.join(data[read_pair][2]))
                    elif int(data[read_pair][1][1]) < 300:
                        summary['Short_insert'] += 1
                        print >> ofile_class, '%s\t%s\t%s\tShort_insert' %(read_pair, '\t'.join(data[read_pair][1]), '\t'.join(data[read_pair][2]))
                    elif int(data[read_pair][1][1]) > 1000:
                        summary['mPing_reads'] += 1
                        mping_reads[read_pair] = 1
                        print >> ofile_class, '%s\t%s\t%s\tmPing_reads' %(read_pair, '\t'.join(data[read_pair][1]), '\t'.join(data[read_pair][2]))
                    else:
                        summary['error2'] += 1
                        print >> ofile_class, '%s\t%s\t%s\terror2' %(read_pair, '\t'.join(data[read_pair][1]), '\t'.join(data[read_pair][2]))
                else:
                    summary['Not_16thSNP_reads'] += 1
                    print >> ofile_class, '%s\t%s\t%s\tNot_16thSNP_reads' %(read_pair, '\t'.join(data[read_pair][1]), '\t'.join(data[read_pair][2]))
        print >> ofile_log, 'Total mapped paired reads: %s' %(summary['paired'])
        print >> ofile_log, 'Ping reads: %s' %(summary['Ping_reads'])
        print >> ofile_log, 'mPing reads: %s' %(summary['mPing_reads'])
        print >> ofile_log, 'Short insert: %s' %(summary['Short_insert'])
        print >> ofile_log, 'error1: %s' %(summary['error1']) 
        print >> ofile_log, 'error2: %s' %(summary['error2'])
        print >> ofile_log, 'Not_16thSNP_reads: %s' %(summary['Not_16thSNP_reads'])
        #overlap SNPreads and Pingreads
        G_reads = snp_reads_from_bam('%s.alt.G.bam' %(prefix), 'G')   
        A_reads = snp_reads_from_bam('%s.ref.A.bam' %(prefix), 'A')
        ping_reads_class = overlap_reads(ping_reads, G_reads, A_reads) 
        mping_reads_class = overlap_reads(mping_reads, G_reads, A_reads)
        print >> ofile_log, 'Ping G-type 16th SNP reads: %s, %s' %(str(len(ping_reads_class['G'].keys())), ','.join(ping_reads_class['G'].keys()))
        print >> ofile_log, 'Ping A-type 16th SNP reads: %s, %s' %(str(len(ping_reads_class['A'].keys())), ','.join(ping_reads_class['A'].keys()))
        print >> ofile_log, 'Ping UN-type 16th SNP reads: %s, %s' %(str(len(ping_reads_class['UN'].keys())), ','.join(ping_reads_class['UN'].keys()))
        print >> ofile_log, 'mPing G-type 16th SNP reads: %s, %s' %(str(len(mping_reads_class['G'].keys())), ','.join(mping_reads_class['G'].keys()))
        print >> ofile_log, 'mPing A-type 16th SNP reads: %s, %s' %(str(len(mping_reads_class['A'].keys())), ','.join(mping_reads_class['A'].keys()))
        print >> ofile_log, 'mPing UN-type 16th SNP reads: %s, %s' %(str(len(mping_reads_class['UN'].keys())), ','.join(mping_reads_class['UN'].keys()))
        #summary
        print >> ofile_sum, '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(strain, summary['Ping_reads'], str(len(ping_reads_class['G'].keys())), str(len(ping_reads_class['A'].keys())), summary['mPing_reads'], str(len(mping_reads_class['G'].keys())), str(len(mping_reads_class['A'].keys()))) 
    ofile_log.close()
    ofile_sum.close()
'''
if __name__ == '__main__':
    main()

