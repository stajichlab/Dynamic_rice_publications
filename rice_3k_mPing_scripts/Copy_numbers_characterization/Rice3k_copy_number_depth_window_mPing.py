#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from scipy import stats

def usage():
    test="name"
    message='''
python Rice3k_copy_number_depth.py --input fq_RelocaTE2_Ping_NM2 --output fq_RelocaTE2_Ping_NM2.Ping_copy.txt

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 100 --lines %s --interval 120 --resource nodes=1:ppn=1,walltime=200:00:00,mem=10G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

#ping    2       N       11      GGGGGGGgg^]G^]g IIIIIIIIIII
def ping_avg_mpileup(infile, genome_depth, name):
    data = defaultdict(lambda : list())
    #data = []
    mping = 0
    win   = 50
    step  = 40
    nwin  = 430/float(step)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if int(unit[1]) >=1 and int(unit[1]) <= 430:
                    #key1 marks starting point of each window
                    key1 = int(int(unit[1])/step)
                    #key2 marks overlapping point between windows and can be used to mark end point of each window
                    key2 = key1 if int(int(unit[1])%step) >= win - step else key1 - 1 
                    if key1 == key2: 
                        data[key1].append([int(unit[1]), int(unit[3])])
                    else:
                        data[key1].append([int(unit[1]), int(unit[3])])
                        data[key2].append([int(unit[1]), int(unit[3])])
                    #print 'Key1: %s, Key2: %s, %s' %(key1, key2, line)
                else:
                    mping += 1
    data_region = []
    for region in sorted(data.keys(), key=int):
        win_depth = []
        for pos in data[region]:
            #print '%s\t%s\t%s' %(region, str(pos[0]), str(pos[1]))
            win_depth.append(pos[1])
        if len(win_depth) >= 30:
            win_mean   = np.mean(win_depth)
            win_median = np.median(win_depth)
            win_std    = np.std(win_depth)
            #print '%s\t%s\t%s\t%s' %(region, win_mean, win_median, win_std)
            data_region.append(win_mean)
    print '%s\t%s\t%s' %(name, genome_depth[name], len(data_region))
    #if len(data_region) < 0.7*(nwin):
    #    return [name, 'NA']
    #target copy number mean and std
    strain_genome_depth = float(genome_depth[name])
    strain_target_depth = 0.0
    strain_mean  = 0.0
    strain_std   = 0.0
    ttest1sample = ['nan', 'nan']
    if len(data_region) > 1:
        strain_target_depth = np.mean(data_region)
        strain_mean  = np.mean([x/strain_genome_depth for x in data_region])
        strain_std   = np.std([x/strain_genome_depth for x in data_region])
        #target copy number one sample t-test
        data_region_array = np.array(data_region)
        ttest1sample      = stats.ttest_1samp(data_region_array, float(strain_genome_depth))
        #print '%s\t%s\t%s\t%s\t%s' %(name, str(strain_target_depth), str(strain_genome_depth), str(strain_mean), str(strain_std))
        output = [name, str(strain_target_depth), str(strain_genome_depth), str(strain_mean), str(strain_std), str(ttest1sample[0]), str(ttest1sample[1])]
    else:
        output = [name, str(strain_target_depth), str(strain_genome_depth), str(strain_mean), str(strain_std), str(ttest1sample[0]), str(ttest1sample[1])]
    #avg, med, sd = 0, 0, 0
    #if len(data) > 0: 
    #    avg = np.sum(data)/430.00
    #    #med = np.median(data)
    #    #sd  = np.std(data) 
    return output

#Sample  #Read   Depth   Mapped_Depth    Mapped_rate
#nivara_IRGC105327       66904761        8.816764895     8.36437303618   0.948689585783
def read_depth(infile):  
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Sample'): 
                unit = re.split(r'\t',line)
                strain = re.sub(r'.realigned', r'', unit[0])
                data[strain] = unit[3]
    return data

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
 
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-d', '--depth')
    parser.add_argument('-l', '--list')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
        len(args.depth) > 0
    except:
        usage()
        sys.exit(2)
    
    if not args.output:
        args.output = '%s.window_copy_number.txt' %(args.input.strip("/"))
 
    if not args.list:
        args.list = 'rice_line_3000.download.list'

    #test_bam(args.input)
    genome_depth = read_depth(args.depth)
    acc2name     = read_download_table(args.list)

    #####summarize coverage into matrix
    #ping    2       N       11      GGGGGGGgg^]G^]g IIIIIIIIIII
    covs = glob.glob('%s/*.NM2.mpileup' %(args.input))
    ofile = open(args.output, 'w')
    print >> ofile, "Taxa\tTarget_Depth\tGenome_Depth\tTarget_Copy_Number_Mean\tTarget_Copy_Number_STD\tT-statistic\tP-value"
    for ping in sorted(covs):
        name  = os.path.split(ping)[1]
        name  = re.sub(r'_RelocaTE2\.NM2.mpileup', r'', name)
        name  = re.sub(r'\.NM2.mpileup', r'', name)
        name  = acc2name[name]
        ping_avg = ping_avg_mpileup(ping, genome_depth, name)
        print >> ofile, '%s' %('\t'.join(ping_avg)) 
    ofile.close()

if __name__ == '__main__':
    main()

