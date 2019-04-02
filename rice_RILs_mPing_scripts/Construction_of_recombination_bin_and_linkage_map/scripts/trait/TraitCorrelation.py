#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
from scipy.stats.stats import pearsonr
import re
import os
import argparse

def usage():
    test="name"
    message='''
Python TraitCorrelation.py --input ../input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.ALL
    '''
    print message

def killnan(list1, list2):
    list3 = []
    list4 = []
    for i in range(len(list1)):
        if isnan(list1[i]) or isnan(list2[i]):
            continue
        else:
            #print list1[i], list2[i]
            list3.append(list1[i]) 
            list4.append(list2[i]) 
    #print len(list1), len(list2), len(list3), len(list4)
    return list3, list4 

def trait(infile):
    data = defaultdict(list)
    corr = defaultdict(lambda : defaultdict(lambda : list()))
    '''read trait table'''
    with open (infile, 'r') as filefh:
        header  = filefh.readline()
        headers = re.split('\t',header.rstrip())
        headers = headers[1:]
        #print headers[0], headers[1]
        for line in filefh:
            line = line.rstrip()
            unit = re.split('\t',line)
            for i in range(1,len(unit)):                
                value = re.sub(r'NA',r'NaN',unit[i])
                #print i, value 
                data[i].append(float(value))
    '''Pearson correlation'''
    for rank in sorted(data.keys()):
        #print rank, ','.join(data[rank])
        for rank2 in range(rank + 1, len(data.keys())+1):
            #print rank, rank2
            list1, list2 = killnan(data[rank],data[rank2])
            correlation = pearsonr(list1,list2)
            #print correlation[0], correlation[1]
            corr[rank][rank2] = [correlation[0], correlation[1]]
    '''output matirx, correlation and p-value'''
    print ",".join(headers)
    for rank in sorted(data.keys()):
        outline = []
        #outline.append(str(rank))
        outline.append(headers[rank-1])
        for rank2 in sorted(data.keys()):
            #print rank, rank2
            #r = corr[rank2][rank][0]
            #p = corr[rank2][rank][1]
            if rank2 < rank:
                r = corr[rank2][rank][0]
                outline.append(r)
            elif rank2 > rank:
                p = corr[rank][rank2][1]
                outline.append(p)
            else:
                outline.append(" ")
        print ",".join(map(str,outline))
        
            

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

    trait(args.input)

if __name__ == '__main__':
    main()

