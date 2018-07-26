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
import urllib2

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message


def runjob(script, lines):
    #cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 100 --lines %s --interval 120 --resource nodes=1:ppn=1,walltime=200:00:00,mem=10G --convert no %s' %(lines, script)
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 100 --lines %s --interval 120 --task 1 --mem 15G --time 10:00:00 --convert no --queue intel %s' %(lines, script)
    #print cmd 
    os.system(cmd)


def read_cov(infile):
    name     = os.path.split(infile)[1]
    coverage = 0
    data     = []
    count   = 0
    covered = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                flag = 0 
                unit = re.split(r'\t',line)
                if unit[0] == 'chr06' and int(unit[1]) >= 23520641 and int(unit[1]) <= 23527981:
                    if int(unit[2]) > 0:
                        data.append('10')
                    else:
                        data.append('0')
                    count += 1
                    if int(unit[2]) > 0:
                        covered += 1
    coverage = covered
    return coverage, data

#samtools view http://s3.amazonaws.com/3kricegenome/Nipponbare/IRIS_313-10271.realigned.bam
def ping_coverage(acc):
    #cmd = '/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.3/bin/samtools view http://s3.amazonaws.com/3kricegenome/Nipponbare/%s.realigned.bam > test.sam' %(acc)
    bam = 'http://s3.amazonaws.com/3kricegenome/Nipponbare/%s.realigned.bam' %(acc)
    file_ok    = 0
    try: 
        urllib2.urlopen(bam)
        file_ok = 1
        print '%s: okay' %(bam)
    except urllib2.HTTPError, e:
        print '%s: %s' %(bam, e.code)
        file_ok = 0
    except urllib2.URLError, e:
        print '%s: %s' %(bam, e.code)
        file_ok = 0
   
    return file_ok 


#Taxa    Color   Label   Name    Origin  Group   mPing   Ping    Pong
#B001    blue    Heibiao|B001|Temp       Heibiao China   Temperate jap   71      1       8
def test_bam(infile):
    file_ok = 0
    file_wrong = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not unit[0] == 'Taxa':
                    file_status = ping_coverage(unit[0])
                    if file_status == 0:
                        file_wrong +=1
                    elif file_status == 1:
                        file_ok    += 1
    print 'Okay: %s' %(file_ok)
    print 'Wrong: %s' %(file_wrong)  


def stowaway_check(infile):
    ofile = open('stowaway_check.sh', 'w')
    count = 0
    samtools = '/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.3/bin/samtools'
    bamtools = '/opt/linux/centos/7.x/x86_64/pkgs/bamtools/2.4.0/bin/bamtools'
    bedtools = '/rhome/cjinfeng/BigData/software/bedtools2-2.19.0/bin/bedtools'
    bedping  = '/rhome/cjinfeng/Rice/Rice_population_sequence/Rice_3000/analysis/Ping_Pong_deletion_derivative/pacbio_NB/fastq/WGS/Rice3k_CAAS/genome.ping.bed'
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Taxa'):
                unit = re.split(r'\t',line)
                acc = unit[0]
                #bam = 'http://s3.amazonaws.com/3kricegenome/Nipponbare/%s.realigned.bam' %(acc)
                bam = '/rhome/cjinfeng/Rice/Rice_population_sequence/Rice_3000/analysis/Ping_Pong_deletion_derivative/pacbio_NB/fastq/WGS/Rice3k_CAAS/ping_coverage_3k/%s.realigned.bam' %(acc)
                subbam = 'stowaway_check/%s.special_ping_locus.bam' %(acc) 
                #coverage  = 'ping_coverage_3k/%s.ping.coverage.txt' %(acc)
                #coveragez = 'ping_coverage_3k/%s.ping.coverage.txt.gz' %(acc)
                #coverage1 = 'ping_coverage_3k/%s.ping.coverage.clean.txt' %(acc)
                #if not os.path.exists(coverage1) or int(os.path.getsize(coverage1)) == 0:
                if not os.path.exists(subbam) or int(os.path.getsize(subbam)) == 0:
                    cmd   = '''%s view -b %s chr01:2540500-2740500 > %s''' %(samtools, bam, os.path.abspath(subbam))
                    index = '''%s index %s''' %(samtools, os.path.abspath(subbam))
                    print >> ofile, cmd
                    print >> ofile, index
                    #print >> ofile, sed
                    count += 1
    ofile.close()
    if count > 0:
        runjob('stowaway_check.sh', 300)
   
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
    
 
    #test_bam(args.input)
    stowaway_check(args.input)


if __name__ == '__main__':
    main()

