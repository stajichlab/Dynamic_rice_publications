#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob

def usage():
    test="name"
    message='''
RelocaTEi: improved version of RelocaTE for calling transposable element insertions

bam mode:
python relocaTE.py --bam MSU7.Chr4.ALL.rep1_reads_2X_100_500.bam --genome_fasta MSU7.Chr4.fa --te_fasta mping.fa --reference_ins MSU_r7.fa.RepeatMasker.Chr4.out --outdir RelocaTE_output_multiTE_bam

fastq mode:
python relocaTE.py --fq_dirMSU7.Chr4.ALL.rep1_reads_5X_100_500 --genome_fasta MSU7.Chr4.fa --te_fasta mping.fa --reference_ins MSU7.Chr4.fa.RepeatMasker.out --outdir RelocaTE_output_mPing_gz

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def createdir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)

def writefile(outfile, lines):
    ofile = open(outfile, 'w')
    print >> ofile, lines
    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam')
    parser.add_argument('-t', '--te_fasta')
    parser.add_argument('-d', '--fq_dir')
    parser.add_argument('-g', '--genome_fasta')
    parser.add_argument('-r', '--reference_ins')
    parser.add_argument('-f', '--fastmap')
    parser.add_argument('-o', '--outdir')
    parser.add_argument('-s', '--size')
    parser.add_argument('-1', '--mate_1_id', default='_1')
    parser.add_argument('-2', '--mate_2_id', default='_2')
    parser.add_argument('-u', '--unpaired_id', default='.unParied')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    
    try:
        os.path.isfile(args.bam) and os.path.isfile(args.te_fasta) and os.path.isfile(args.genome_fasta)
    except:
        try:
            os.path.exists(args.fq_dir) and os.path.isfile(args.te_fasta) and os.path.isfile(args.genome_fasta)
        except:
            usage()
            sys.exit(2)

    #absolute path for all files, directory and scripts
    RelocaTE_bin = os.path.split(os.path.abspath(__file__))[0]
    reference    = os.path.abspath(args.genome_fasta)
    te_fasta     = os.path.abspath(args.te_fasta)
    mode         = 'bam'
    bam          = ''
    fastq_dir    = ''
    
    try: 
        os.path.isfile(args.bam)
        bam  = os.path.abspath(args.bam)
    except:
        try:
            if os.path.abspath(args.fq_dir):
                fastq_dir = os.path.abspath(args.fq_dir)
                mode      = 'fastq'
        except:
            usage()
            exit(2)

    print fastq_dir
    print mode

    #Prepare directory and script
    if args.outdir is None:
        args.outdir = '%s/RelocaTE_output' %(os.getcwd())
        createdir(args.outdir)
    else:
        args.outdir = os.path.abspath(args.outdir)
        createdir(args.outdir)

    if args.size is None:
        args.size = 500

    samtools = ''
    bedtools = ''
    bwa      = ''
    
    try:
        subprocess.check_output('which samtools', shell=True)
        samtools = subprocess.check_output('which samtools', shell=True)
        samtools = re.sub(r'\n', '', samtools)
    except:
        samtools = '/opt/samtools-0.1.16/samtools'

    try:
        subprocess.check_output('which bedtools', shell=True)
        bedtools = subprocess.check_output('which bedtools', shell=True)
        bedtools = re.sub(r'\n', '', bedtools)
    except:
        bedtools = '/opt/bedtools/2.17.0-25-g7b42b3b/bin//bedtools'

    try:
        subprocess.check_output('which bwa', shell=True)
        bwa = subprocess.check_output('which bwa', shell=True)
        bwa = re.sub(r'\n', '', bwa)
    except:
        bwa = '/opt/tyler/bin/bwa'
 
    try:
        subprocess.check_output('which blat', shell=True)
        blat = subprocess.check_output('which blat', shell=True)
        blat = re.sub(r'\n', '', blat)
    except:
        blat = '/usr/local/bin/blat'   
    
    try:
        subprocess.check_output('which seqtk', shell=True)
        seqtk = subprocess.check_output('which seqtk', shell=True)
        seqtk = re.sub(r'\n', '', seqtk)
    except:
        seqtk = '/rhome/cjinfeng/software/tools/seqtk-master//seqtk'

    #overwrite tools
    blat = '/opt/linux/centos/7.x/x86_64/pkgs/blat/35/bin/blat'
    #bwa = '/opt/bwa/0.7.9/bin/bwa'
    bwa  = '/rhome/cjinfeng/BigData/00.RD/RelocaTE2/tools/bwa-0.6.2/bwa'
    bowtie2  = '/opt/bowtie2/2.2.3/bowtie2'
    bedtools = '/opt/bedtools/2.17.0-25-g7b42b3b/bin//bedtools'
    samtools = '/opt/samtools/0.1.19/bin/samtools'
    seqtk = '/rhome/cjinfeng/BigData/software/seqtk-master/seqtk'
    fastq_split = '%s/fastq_split.pl' %(RelocaTE_bin)

    #MSU_r7.fa.bwt
    if not os.path.isfile('%s.bwt' %(reference)):
        print 'Reference need to be indexed by bwa: %s' %(reference)
        exit()
 
    run_std = '%s/run.std' %(args.outdir)

    writefile('%s/regex.txt' %(args.outdir), '%s\t%s\t%s\tUNK' %(args.mate_1_id, args.mate_2_id, args.unpaired_id))
    createdir('%s/shellscripts' %(args.outdir))
    createdir('%s/repeat' %(args.outdir))
    createdir('%s/repeat/blat_output' %(args.outdir))
    createdir('%s/repeat/flanking_seq' %(args.outdir))
    createdir('%s/repeat/te_containing_fq' %(args.outdir))
    createdir('%s/repeat/te_only_read_portions_fa' %(args.outdir))

    shells = []
    #step0 existing TE blat
    reference_ins_flag = 'NONE'
    if args.reference_ins is None or args.reference_ins == '0':
        step0_file = '%s/shellscripts/step_0_do_not_call_reference_insertions' %(args.outdir)
        writefile(step0_file, '')
    elif os.path.isfile(args.reference_ins):
        step0_file = '%s/shellscripts/step_0_te_annotation_provided' %(args.outdir)
        reference_ins_flag = args.reference_ins
        writefile(step0_file, '')
    elif args.reference_ins == '1':
        createdir('%s/shellscripts/step_0' %(args.outdir))
        step0_file = '%s/shellscripts/step_0/step_0.existingTE_blat.sh' %(args.outdir)
        shells.append('sh %s' %(step0_file))
        existingTE_blat = '%s %s %s %s/existingTE.blatout 1> %s/existingTE.blat.stdout' %(blat, reference, te_fasta, args.outdir, args.outdir)
        reference_ins_flag = '%s/existingTE.blatout' %(args.outdir)
        writefile(step0_file, existingTE_blat)

    #step1 format reference genome

    #step2 fastq to fasta
    fastas = defaultdict(lambda : str)
    if mode == 'fastq':
        fastqs = glob.glob('%s/*.f*q*' %(fastq_dir))
        step2_flag = 0
        step2_count= 0
        for fq in fastqs:
            fa    = ''
            if os.path.splitext(fq)[-1] == '.gz':
                fa    = '%s.fa' %(os.path.splitext(os.path.splitext(fq)[0])[0])
            else:
                fa    = '%s.fa' %(os.path.splitext(fq)[0])
            fastas[fa] = fq
            if not os.path.isfile(fa):
                createdir('%s/shellscripts/step_2' %(args.outdir))
                fq2fa = '%s seq -A %s > %s' %(seqtk, fq, fa)
                step2_file = '%s/shellscripts/step_2/%s.fq2fq.sh' %(args.outdir, step2_count)
                shells.append('sh %s' %(step2_file))
                writefile(step2_file, fq2fa)
                step2_flag == 1
                step2_count += 1
        if step2_flag == 0:
            step2_file = '%s/shellscripts/step_2_not_needed_fq_already_converted_2_fa' %(args.outdir)
            writefile(step2_file, '')        
    elif mode == 'bam':
        #print 'Add module of obtaining reads from bam then prepare as fa files'
        cmd_step2 = []
        fastq_dir = '%s/repeat/fastq' %(args.outdir)
        createdir(fastq_dir)
        subbam = '%s/%s.subset.bam' %(fastq_dir, os.path.splitext(os.path.split(bam)[1])[0])
        fq1 = '%s/%s_1.fq' %(fastq_dir, os.path.splitext(os.path.split(bam)[1])[0])
        fq2 = '%s/%s_2.fq' %(fastq_dir, os.path.splitext(os.path.split(bam)[1])[0])
        fa1 = '%s.fa' %(os.path.splitext(fq1)[0])
        fa2 = '%s.fa' %(os.path.splitext(fq2)[0])
        fastas[fa1] = fq1
        fastas[fa2] = fq2
        if not os.path.isfile(subbam):
            #cmd_step2.append('%s view -h %s | awk \'$5<60\' | samtools view -Shb - | samtools sort -m 500000000 -n - %s 2> %s' %(samtools, bam, os.path.splitext(subbam)[0], run_std))
            cmd_step2.append('%s sort -m 1000000000 -n %s %s 2> %s' %(samtools, bam, os.path.splitext(subbam)[0], run_std))
        if not os.path.isfile(fq1) and not os.path.isfile(fq2):
            cmd_step2.append('%s bamtofastq -i %s -fq %s -fq2 %s 2> %s' %(bedtools, subbam, fq1, fq2, run_std))
        cmd_step2.append('%s seq -A %s > %s' %(seqtk, fq1, fa1))
        cmd_step2.append('%s seq -A %s > %s' %(seqtk, fq2, fa2))
        step2_flag = 0
        if not os.path.isfile(fa1) and not os.path.isfile(fa2):
            createdir('%s/shellscripts/step_2' %(args.outdir))
            step2_file = '%s/shellscripts/step_2/0.bam2fa.sh' %(args.outdir)
            shells.append('sh %s' %(step2_file))
            writefile(step2_file, '\n'.join(cmd_step2))
        else:
            step2_file = '%s/shellscripts/step_2_not_needed_fq_already_converted_2_fa' %(args.outdir)
            writefile(step2_file, '')

    #step3 blat fasta to repeat
    step3_count = 0
    for fa in sorted(fastas.keys()):
        createdir('%s/shellscripts/step_3' %(args.outdir))
        fq      = fastas[fa]
        #fq      = '%s.fq' %(os.path.splitext(fa)[0]) if os.path.isfile('%s.fq' %(os.path.splitext(fa)[0])) else '%s.fastq' %(os.path.splitext(fa)[0])
        fa_prefix = os.path.split(os.path.splitext(fa)[0])[1]
        blatout = '%s/repeat/blat_output/%s.te_repeat.blatout' %(args.outdir, fa_prefix)
        blatstd = '%s/repeat/blat_output/blat.out' %(args.outdir)
        blatcmd = '%s -minScore=10 -tileSize=7 %s %s %s 1>>%s 2>>%s' %(blat ,te_fasta, fa, blatout, blatstd, blatstd)
        #if args.fastmap is not None:
        #    blat = 'blat -minScore=10 -tileSize=7 -fastMap %s %s %s 1>> %s' %(te_fasta, fa, blatout, blatstd)
        flank= '%s/repeat/flanking_seq/%s.te_repeat.flankingReads.fq' %(args.outdir, fa_prefix)
        trim = 'python %s/relocaTE_trim.py %s %s 10 1 > %s' %(RelocaTE_bin, blatout, fq, flank)
        step3_file = '%s/shellscripts/step_3/%s.te_repeat.blat.sh' %(args.outdir, step3_count)
        if not os.path.isfile(blatout) or os.path.getsize(blatout) == 0:
            shells.append('sh %s' %(step3_file))
            step3_cmds = '%s\n%s' %(blatcmd, trim)
            writefile(step3_file, step3_cmds)
            step3_count += 1
        elif not os.path.isfile(flank) or os.path.getsize(flank) == 0:
            shells.append('sh %s' %(step3_file))
            step3_cmds = '%s' %(trim)
            writefile(step3_file, step3_cmds)
            step3_count += 1

    #step4 align TE trimed reads to genome
    ref = os.path.split(os.path.splitext(reference)[0])[1]
    createdir('%s/shellscripts/step_4' %(args.outdir))
    step4_file= '%s/shellscripts/step_4/step_4.%s.repeat.align.sh' %(args.outdir, ref)
    shells.append('sh %s' %(step4_file))
    step4_cmd = 'python %s/relocaTE_align.py %s %s/repeat %s %s %s/regex.txt repeat not.given 0' %(RelocaTE_bin, RelocaTE_bin, args.outdir, reference, fastq_dir, args.outdir)
    writefile(step4_file, step4_cmd)
    
    #step5 find insertions
    ids = fasta_id(reference)
    createdir('%s/shellscripts/step_5' %(args.outdir))
    step5_count = 0
    for chrs in ids:
        step5_cmd = 'python %s/relocaTE_insertionFinder.py %s/repeat/bwa_aln/%s.repeat.bwa.sorted.bam %s %s repeat %s/regex.txt not.give 100 %s 0 0 %s' %(RelocaTE_bin, args.outdir, ref, chrs, reference, args.outdir, reference_ins_flag, args.size)
        step5_file= '%s/shellscripts/step_5/%s.repeat.findSites.sh' %(args.outdir, step5_count)
        shells.append('sh %s' %(step5_file))
        writefile(step5_file, step5_cmd)
        step5_count +=1
    
    #step6 find transposons on reference: reference only or shared
    createdir('%s/shellscripts/step_6' %(args.outdir))
    step6_count = 0
    if mode == 'fastq' or mode == 'bam':
        for chrs in ids:
            step6_cmd = 'python %s/relocaTE_absenceFinder.py %s/repeat/bwa_aln/%s.repeat.bwa.sorted.bam %s %s repeat %s/regex.txt not.give 100 %s 0 0 %s' %(RelocaTE_bin, args.outdir, ref, chrs, reference, args.outdir, reference_ins_flag, args.size)
            step6_file= '%s/shellscripts/step_6/%s.repeat.absence.sh' %(args.outdir, step6_count)
            shells.append('sh %s' %(step6_file))
            writefile(step6_file, step6_cmd)
            step6_count +=1
    #elif mode == 'bam':
    #    pass
 
    #step7 characterize homozygous, heterozygous and somatic insertion
    createdir('%s/shellscripts/step_7' %(args.outdir))
    step7_cmd = []
    step7_cmd.append('cat %s/repeat/results/*.all_nonref_insert.gff > %s/repeat/results/ALL.all_nonref_insert.gff' %(args.outdir, args.outdir))
    step7_cmd.append('cat %s/repeat/results/*.all_nonref_insert.txt | grep "^TE" -v > %s/repeat/results/ALL.all_nonref_insert.txt' %(args.outdir, args.outdir))
    step7_cmd.append('cat %s/repeat/results/*.all_ref_insert.txt > %s/repeat/results/ALL.all_ref_insert.txt' %(args.outdir, args.outdir))
    step7_cmd.append('cat %s/repeat/results/*.all_ref_insert.gff > %s/repeat/results/ALL.all_ref_insert.gff' %(args.outdir, args.outdir))
    if os.path.exists(bam):
        step7_cmd.append('perl %s/characterizer.pl -s %s/repeat/results/ALL.all_nonref_insert.txt -b %s -g %s --samtools %s' %(RelocaTE_bin, args.outdir, bam, reference, samtools))
    step7_file= '%s/shellscripts/step_7/0.repeat.characterize.sh' %(args.outdir)
    shells.append('sh %s' %(step7_file)) 
    writefile(step7_file, '\n'.join(step7_cmd))


    #step8 clean temp files
    shells_clean = []
    clean_cmd = []
    clean_cmd.append('rm %s/*.fa' %(fastq_dir))
    #clean_cmd.append('rm -R %s/repeat/blat_output' %(args.outdir))
    #clean_cmd.append('rm -R %s/repeat/flanking_seq' %(args.outdir))
    #clean_cmd.append('rm -R %s/repeat/te_containing_fq' %(args.outdir))
    #clean_cmd.append('rm -R %s/repeat/te_only_read_portions_fa' %(args.outdir))
    #clean_cmd.append('rm -R %s/repeat/fastq_split' %(args.outdir))
    #clean_cmd.append('rm %s/repeat/bwa_aln/*.mates.bam* %s/repeat/bwa_aln/*.unPaired.bam* %s/repeat/bwa_aln/*.bwa.bam*' %(args.outdir, args.outdir, args.outdir))
    clean_file = '%s/clean_intermediate_files.sh' %(args.outdir)
    writefile(clean_file, '\n'.join(clean_cmd))
    shells.append('sh %s' %(clean_file))
    shells_clean.append('sh %s' %(clean_file))


    #write script
    writefile('%s/run_these_jobs.sh' %(args.outdir), '\n'.join(shells))

if __name__ == '__main__':
    main()

