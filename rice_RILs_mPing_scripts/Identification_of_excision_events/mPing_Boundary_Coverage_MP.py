#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import subprocess
import multiprocess as mp
sys.path.append('%s/lib' %(os.getcwd()))
from excision import bamcheck, bamcheck_simple, convert_MAP, genotyping, convert_MAP_SNP, genotyping_SNP
import glob

def usage():
    test="name"
    message='''
python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam --gff_ref ../input/HEG4.ALL.mping.non-ref.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing.gff

Check read coverage at mPing boundary from bam files of read2reference and read2mping_flanking. We will create matrix
of RILs and mPing, which gives information of each end of mPing.

    '''
    print message

#Chr1    PseudoGenome    Transposable_element    1132975 1133405 -       .       .       ID=Chr1_1132975_1133405;Original_ID=Chr1.1132977.spanners;TE=mping;TSD=TAA;
#Chr1    PseudoGenome    Transposable_element    2642232 2647573 +       .       .       ID=Chr1_2642232_2647573;Original_ID=Chr1.2640500;TE=ping;TSD=TAA;
def id_mapping(infile, mping2ID_0):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                start = int(unit[3]) 
                end   = int(unit[4])
                chro  = unit[0]
                strand= unit[6]
                temp  = defaultdict(str)
                attrs = re.split(r';', unit[8])
                for attr in attrs:
                    if not attr == '':
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value
                temp['Original_ID'] = re.sub(r'.spanners', r'', temp['Original_ID'])
                temp_ids   = re.split(r'_', temp['ID'])
                temp['ID'] = '%s:%s-%s' %(temp_ids[0], temp_ids[1], temp_ids[2])
                temp_oids  = re.split(r'\.', temp['Original_ID'])
                #print temp_oids
                temp['Original_ID'] = '%s:%s-%s' %(temp_oids[0], str(int(temp_oids[1])-2), temp_oids[1])
                mping2ID_0[temp['ID']]  = temp['Original_ID']
                #print '%s\t%s' %(temp['ID'], mping2ID_0[temp['ID']])
    return data

def bamcheck_ref(bam, mping, bamck_file, ril):
    reg     = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match   = reg.search(mping)
    chro    = match.groups(0)[0]
    l_start = int(match.groups(0)[1]) - 2
    l_end   = int(match.groups(0)[1]) + 2
    l_mping = '%s:%s-%s' %(chro, l_start, l_end)
    r_start = int(match.groups(0)[2]) - 2
    r_end   = int(match.groups(0)[2]) + 2
    r_mping = '%s:%s-%s' %(chro, r_start, r_end)
    l_flag  = bamcheck_simple(bam, l_mping, bamck_file, 5)
    r_flag  = bamcheck_simple(bam, r_mping, bamck_file, 3)
    return (l_flag, r_flag)

def decode(flag):
    if int(flag) == 0 or int(flag) == 4:
        return 'covered'
    elif int(flag) == 1:
        return 'clipped'
    else:
        return 'unknown'

def decode_gt(gt):
    if gt == 'NA':
        return 'NA'
    elif str(gt) == '0':
        return 'NB'
    elif str(gt) == '1':
        return 'HEG4'
    else:
        return 'NA'

def sort_mping_chr(mpings):
    data = defaultdict(lambda : defaultdict(lambda : str()))
    for mping in mpings.values():
        reg     = re.compile(r'Chr(\d+):(\d+)\-(\d+)')
        match   = reg.search(mping)
        chro    = match.groups(0)[0]
        start   = match.groups(0)[1]
        end     = match.groups(0)[2]
        data[chro][start] = mping
    sorted_mping = []
    for c in sorted(data.keys(), key=int):
        for s in sorted(data[c].keys(), key=int):
            sorted_mping.append(data[c][s])
    return sorted_mping

def get_rils(bams):
    rils = []
    for bam in bams:
        ril = os.path.split(bam)[1]
        ril = re.sub(r'.bam', r'', ril)
        ril = re.sub(r'RIL', r'', ril)
        rils.append(ril)
    return rils

##run function with parameters using multiprocess of #cpu
def multiprocess_pool(parameters, cpu):
    pool = mp.Pool(int(cpu))
    imap_it = pool.map(mping_genotyper_mp_helper, tuple(parameters))
    collect_list = []
    for x in imap_it:
        #print 'status: %s' %(x)
        collect_list.append(x)
    #clean up
    #https://timothyawiseman.wordpress.com/2012/12/21/a-really-simple-multiprocessing-python-example/
    pool.close()
    pool.join()
    return collect_list

def mping_genotyper_mp_helper(args):
    return mping_genotyper_mp_runner(*args)

def mping_genotyper_mp_runner(ril, bam_ref, bam_pseudo, mping2ID_0, binmap, snpmap, bamcheck_file_ref, bamcheck_file_pseudo):
    genotype = 3
    l_flag = 3
    r_flag = 3
    ref_flag = 3
    results  = []
    bamcheck_file_pseudo = '%s.%s.txt' %(bamcheck_file_pseudo, ril)
    bamcheck_file_ref    = '%s.%s.txt' %(bamcheck_file_ref, ril) 
    if os.path.isfile(bam_ref):
        if os.path.isfile(bam_pseudo):
            for mping in sorted(mping2ID_0.keys()):
                #genotype: 0 ref, 1 non_ref, 3 unknown
                genotype = genotyping_SNP(ril, mping2ID_0[mping], binmap, snpmap)
                l_flag, r_flag = bamcheck_ref(bam_pseudo, mping, bamcheck_file_pseudo, ril)
                ref_flag       = bamcheck(bam_ref, mping2ID_0[mping], bamcheck_file_ref, ril)
                results.append([mping, genotype, l_flag, r_flag, ref_flag]) 
                print '%s\t%s\t%s\t%s' %(ril, mping, l_flag, r_flag)
        else:
            print 'pseudo bam file not found for rils: RIL%s' %(ril)   
    else:
        print 'reference bam file not found for rils: RIL%s' %(ril)
    return [ril, results]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_ref')
    parser.add_argument('--bam_pseudo')
    parser.add_argument('--gff_ref')
    parser.add_argument('--gff_pseudo')
    parser.add_argument('--bin_map')
    parser.add_argument('--snp_map')
    parser.add_argument('--cpu') 
    parser.add_argument('--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.bam_pseudo) > 0 and len(args.bam_ref) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = 'mPing_boundary'
    if not args.gff_ref:
        args.gff_ref = '../input/HEG4.ALL.mping.non-ref.gff'
    if not args.gff_pseudo:
        args.gff_pseudo = '../input/MSU_r7.Pseudo_mPing.gff'
    if not args.bam_ref:
        args.bam_ref = '../input/RILs_ALL_bam'
    if not args.bam_pseudo:
        args.bam_pseudo = '../input/RILs_ALL_unmapped_mping_bam'
    if not args.bin_map:
        args.bin_map = 'MPR.geno.bin'
    if not args.snp_map:
        args.snp_map = 'MPR.geno.data'
    if not args.cpu:
        args.cpu = 16

    cpu = args.cpu
    #we use mping gff from pseudogenome as key to project to everything 
    mping2ID_0  = defaultdict(lambda : str()) #ID_0 is the mping id from original call in HEG4, Chr1.1132977
    id_mapping(args.gff_pseudo, mping2ID_0)
    
    #bin map and snp genotype
    binmap = convert_MAP(args.bin_map)
    snpmap = convert_MAP_SNP(args.snp_map)

    #go through RILs, for each ril determine the status of each mPing
    bamcheck_file_pseudo = '%s.bamcheck_pseudo' %(args.project)
    bamcheck_file_ref    = '%s.bamcheck_ref' %(args.project)
    mping_status      = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str())))
    mping_status_ref  = defaultdict(lambda : defaultdict(lambda : str()))
    mping_bin_gt  = defaultdict(lambda : defaultdict(lambda : str()))
    mping_snp_gt  = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str())))
    #rils = [1, 2, 3, 4]
    bams = glob.glob('%s/*.bam' %(args.bam_pseudo))
    rils = get_rils(bams)
    parameters = []
    for ril in sorted(rils):
        bam_ref    = '%s/GN%s.bam' %(args.bam_ref, ril)
        bam_pseudo = '%s/RIL%s.bam' %(args.bam_pseudo, ril)
        print 'ril: %s, %s' %(ril, bam_pseudo)
        parameters.append([ril, bam_ref, bam_pseudo, mping2ID_0, binmap, snpmap, bamcheck_file_ref, bamcheck_file_pseudo])

    #run multiprocesses
    collect_results = multiprocess_pool(parameters, cpu)
    #process results
    for ril, results in collect_results:
        for mping, genotype, l_flag, r_flag, ref_flag in results:
            print '%s\t%s\t%s\t%s\t%s\t%s' %(ril, mping, genotype, l_flag, r_flag, ref_flag)
            mping_status[ril][mping2ID_0[mping]]['up']   = decode(l_flag)
            mping_status[ril][mping2ID_0[mping]]['down'] = decode(r_flag)
            mping_status_ref[ril][mping2ID_0[mping]]     = decode(ref_flag)
            mping_bin_gt[ril][mping2ID_0[mping]] = decode_gt(genotype)

    #output matrix into file
    matrix_file = '%s.mping_status.matrix.txt' %(args.project)
    ofile = open(matrix_file, 'w')
    mping_ranked= sort_mping_chr(mping2ID_0)
    #mping names
    #print >> ofile, 'mPing,%s' %(','.join(mping_ranked))
    mping_lines = ['mPing']
    for m in mping_ranked:
        mping_lines.append('%s,Genotype,Pseudo_up,Pseudo_down,Ref' %(m))
    print >> ofile, ','.join(mping_lines)

    #mping genotype
    #mping status, matirx
    for ril in sorted(mping_status.keys(), key=int):
        inf_line = ['RIL%s' %(ril)]
        for mping in mping_ranked:
            inf_line.append(mping_bin_gt[ril][mping])
            inf_line.append(mping_status[ril][mping]['up'])
            inf_line.append(mping_status[ril][mping]['down'])
            inf_line.append(mping_status_ref[ril][mping])
            #status = '%s:%s' %(mping_status[ril][mping]['up'], mping_status[ril][mping]['down'])
            #inf_line.append(status)
        print >> ofile, ','.join(inf_line)
    ofile.close()

    #output matrix for individual mPing file
    outdir = '%s_mPing' %(os.path.abspath(args.project))
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    sum_lines_matrix = defaultdict(lambda : list())
    for mping in mping_ranked:
        sum_gt           = [0, 0, 0]
        sum_ref_coverage = defaultdict(lambda : defaultdict(lambda : int()))
        sum_ref_covered  = defaultdict(lambda : defaultdict(lambda : int()))
        sum_ref_clipped  = defaultdict(lambda : defaultdict(lambda : int()))
        sum_nonref_coverage = defaultdict(lambda : defaultdict(lambda : int()))
        sum_nonref_covered  = defaultdict(lambda : defaultdict(lambda : int()))
        sum_nonref_clipped  = defaultdict(lambda : defaultdict(lambda : int()))
        mping_name = re.sub(r':', r'_', mping)
        mping_name = re.sub(r'-', r'_', mping_name)
        ofile = open('%s/%s.matrix.csv' %(outdir, mping_name), 'w')
        #print >> ofile, 'RILs\tDepth(X)\tGenotype_Bin\tDistance_SNP5\tGenotype_SNP5\tDistance_SNP3\tGenotype_SNP3\tmPing_status'
        print >> ofile, '%s,Genotype_Bin,Pseudo_mPing_status_up,Pseudo_mPing_status_down,Ref_mPing_status' %(mping)
        for ril in sorted(mping_status.keys(), key=int):
            #status in pseudo and ref genome have different results
            #in pseudogenome, cover mean junction was covered by reads, which indicates insertion
            #in refgenome, cover mean junction was covered by reads, which indicates no insertion or excision
            #status     = '%s:%s' %(mping_status[ril][mping]['up'], mping_status[ril][mping]['down'])
            status_up  = mping_status[ril][mping]['up']
            status_down= mping_status[ril][mping]['down']
            status_ref = mping_status_ref[ril][mping]
            print >> ofile, 'RIL%s,%s,%s,%s,%s' %(ril, mping_bin_gt[ril][mping], status_up, status_down, status_ref)
            if mping_bin_gt[ril][mping] == 'NB':
                sum_gt[0] += 1
            elif mping_bin_gt[ril][mping] == 'HEG4':
                sum_gt[1] += 1
            else:
                sum_gt[2] += 1
            if 1:
                gt = mping_bin_gt[ril][mping]
                #summary on mping status in ref mapping bam, 0 mean only one status for each mPing
                if status_ref == 'covered':
                    sum_ref_coverage[gt][0] += 1
                    sum_ref_covered[gt][0]  += 1
                elif status_ref == 'clipped':
                    sum_ref_coverage[gt][0] += 1
                    sum_ref_clipped[gt][0]  += 1
                #summary on mping status in pseudo mapping bam, 0 mean upstream status and 1 mean downstream status
                if mping_status[ril][mping]['up'] == 'covered':
                    sum_nonref_coverage[gt][0] += 1
                    sum_nonref_covered[gt][0]  += 1
                elif mping_status[ril][mping]['up'] == 'clipped':
                    sum_nonref_coverage[gt][0] += 1
                    sum_nonref_clipped[gt][0]  += 1
                if mping_status[ril][mping]['down'] == 'covered':
                    sum_nonref_coverage[gt][1] += 1
                    sum_nonref_covered[gt][1]  += 1
                elif mping_status[ril][mping]['down'] == 'clipped':
                    sum_nonref_coverage[gt][1] += 1
                    sum_nonref_clipped[gt][1]  += 1
        #total line, perfectage of coverage
        total_gt = (sum_gt[0] + sum_gt[1])/float(sum(sum_gt))
        total_up   = sum_nonref_coverage['NB'][0]/float(sum(sum_gt)) + sum_nonref_coverage['HEG4'][0]/float(sum(sum_gt)) + sum_nonref_coverage['NA'][0]/float(sum(sum_gt))
        total_down = sum_nonref_coverage['NB'][1]/float(sum(sum_gt)) + sum_nonref_coverage['HEG4'][1]/float(sum(sum_gt)) + sum_nonref_coverage['NA'][1]/float(sum(sum_gt))
        total_ref  = sum_ref_coverage['NB'][0]/float(sum(sum_gt)) + sum_ref_coverage['HEG4'][0]/float(sum(sum_gt)) + sum_ref_coverage['NA'][0]/float(sum(sum_gt))
        print >> ofile, 'Total,%s,%s,%s,%s' %(total_gt, total_up, total_down, total_ref)
        sum_lines_matrix['total'].append('%s,%s,%s,%s' %(total_gt, total_up, total_down, total_ref)) 
        #NB line, covered/clipped/unknown
        nb_gt   = sum_gt[0]
        nb_up   = '%s:%s:%s' %(sum_nonref_covered['NB'][0], sum_nonref_clipped['NB'][0], sum_gt[0]-sum_nonref_covered['NB'][0]-sum_nonref_clipped['NB'][0])
        nb_down = '%s:%s:%s' %(sum_nonref_covered['NB'][1], sum_nonref_clipped['NB'][1], sum_gt[0]-sum_nonref_covered['NB'][1]-sum_nonref_clipped['NB'][1])
        nb_ref  = '%s:%s:%s' %(sum_ref_covered['NB'][0], sum_ref_clipped['NB'][0], sum_gt[0]-sum_ref_covered['NB'][0]-sum_ref_clipped['NB'][0])
        print >> ofile, 'NB,%s,%s,%s,%s' %(nb_gt, nb_up, nb_down, nb_ref)
        sum_lines_matrix['nb'].append('%s,%s,%s,%s' %(nb_gt, nb_up, nb_down, nb_ref))
        #HEG4 line, covered/clipped/unknown
        heg4_gt = sum_gt[1]
        heg4_up   = '%s:%s:%s' %(sum_nonref_covered['HEG4'][0], sum_nonref_clipped['HEG4'][0], sum_gt[1]-sum_nonref_covered['HEG4'][0]-sum_nonref_clipped['HEG4'][0])
        heg4_down = '%s:%s:%s' %(sum_nonref_covered['HEG4'][1], sum_nonref_clipped['HEG4'][1], sum_gt[1]-sum_nonref_covered['HEG4'][1]-sum_nonref_clipped['HEG4'][1])
        heg4_ref  = '%s:%s:%s' %(sum_ref_covered['HEG4'][0], sum_ref_clipped['HEG4'][0], sum_gt[1]-sum_ref_covered['HEG4'][0]-sum_ref_clipped['HEG4'][0])
        print >> ofile, 'HEG4,%s,%s,%s,%s' %(heg4_gt, heg4_up, heg4_down, heg4_ref)
        sum_lines_matrix['heg4'].append('%s,%s,%s,%s' %(heg4_gt, heg4_up, heg4_down, heg4_ref))
        #NA line, covered/clipped/unknown
        na_gt = sum_gt[2]
        na_up   = '%s:%s:%s' %(sum_nonref_covered['NA'][0], sum_nonref_clipped['NA'][0], sum_gt[2]-sum_nonref_covered['NA'][0]-sum_nonref_clipped['NA'][0])
        na_down = '%s:%s:%s' %(sum_nonref_covered['NA'][1], sum_nonref_clipped['NA'][1], sum_gt[2]-sum_nonref_covered['NA'][1]-sum_nonref_clipped['NA'][1])
        na_ref  = '%s:%s:%s' %(sum_ref_covered['NA'][0], sum_ref_clipped['NA'][0], sum_gt[2]-sum_ref_covered['NA'][0]-sum_ref_clipped['NA'][0])
        print >> ofile, 'NA,%s,%s,%s,%s' %(na_gt, na_up, na_down, na_ref)
        sum_lines_matrix['na'].append('%s,%s,%s,%s' %(na_gt, na_up, na_down, na_ref))
        ofile.close()
    #add summary to big matrix
    ofile = open(matrix_file, 'a')
    print >> ofile, '%s,%s' %('Total', ','.join(sum_lines_matrix['total']))
    print >> ofile, '%s,%s' %('NB',    ','.join(sum_lines_matrix['nb']))
    print >> ofile, '%s,%s' %('HEG4',  ','.join(sum_lines_matrix['heg4']))
    print >> ofile, '%s,%s' %('NA',    ','.join(sum_lines_matrix['na']))
    ofile.close()

if __name__ == '__main__':
    main()

