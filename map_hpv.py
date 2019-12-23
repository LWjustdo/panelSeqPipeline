#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/11/7
@Description :HPV/HBV/EBV基因组比对，排序，去重，
'''

import os
import time
import sys
from profile import Profile
import functools
import pandas as pd
import csv

def timefly(func):
    @functools.wraps(func)
    def wrapper(*args,**kw):
        s=time.time()
        res = func(*args,**kw)
        e=time.time()
        print('{} runtime: {}'.format(func.__name__,TransTime(e-s)))
        return res
    return wrapper

def TransTime(seconds):
    h = seconds//3600
    m = seconds%3600//60
    s = seconds%60
    return '{}h {}min {:.0f}s'.format(h,m,s)

class MapSortRedupBQSR():
    def __init__(self,prefix,threads):
        self.var_path = Profile()
        self.prefix = prefix
        self.threads = threads
        self.path = 'MappingHPV'
        self.genome = '/home/longzhao/panelSel/src/hpv_hbv_ebv_genome.fa'

    def mapping(self):
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        os.system("{bwa} mem -t {th} -M -Y -R \'@RG\\tID:{T}\\tPL:ILLUMINA\\tLB:{T}\\tSM:{T}\'  "
                         "{genome} DataQC/InterFiles/{T}_R1.clean.fastq.gz DataQC/InterFiles/{T}_R2.clean.fastq.gz  1> MappingHPV/{T}.sam "
                  "2>> MappingHPV/{T}.map.log".format(T=self.prefix,th=self.threads,genome=self.genome,**self.var_path))

    def sort(self):
        os.system("{samtools} view -@ {th} -bS MappingHPV/{T}.sam -o MappingHPV/{T}.bam >>MappingHPV/{T}.map.log 2>&1".format(T=self.prefix,th=self.threads,**self.var_path))
        os.system("{samtools} sort -@ {th} MappingHPV/{T}.bam MappingHPV/{T}.sort >>MappingHPV/{T}.map.log 2>&1".format(T=self.prefix,th=self.threads,**self.var_path))

    def redup(self):
        os.system("{gatk} MarkDuplicates -I MappingHPV/{T}.sort.bam  -O MappingHPV/{T}.sort.mrkdup.bam --REMOVE_DUPLICATES true "
              "-M MappingHPV/{T}.MarkDuplicates.xls >> MappingHPV/{T}.map.log 2>&1".format(T=self.prefix,**self.var_path))
        os.system("rm MappingHPV/{T}.sam MappingHPV/{T}.bam".format(T=self.prefix))

    def blat(self):
        os.system("/data/biosoft/samtools-1.4/samtools fasta MappingHPV/{T}.sort.mrkdup.bam > MappingHPV/{T}.fasta".format(T=self.prefix))
        os.system("/data/biosoft/blat {genome} MappingHPV/{T}.fasta -out=blast8 MappingHPV/{T}.blast8".format(T=self.prefix,genome=self.genome,**self.var_path))

    @timefly
    def MapPipeline(self):
        self.mapping()
        self.sort()
        self.redup()
        self.blat()


def trans_col(filename):
    df = pd.read_csv(filename,sep='\t',header=None)
    df.columns=['Query id','Subject id','identity','alignment length','mismatches','gap openings','q. start','q. end','s. start','s. end','e-value','bit score']
    df['succeed alignment']=df['alignment length']-df['mismatches']-df['gap openings']
    #df=df.sort_values('identity',ascending=False)
    df=df[df['identity'] >=80]
    df=df[df['succeed alignment']>=70]
    #df=df.sort_values('succeed alignment',ascending=False)
    return df

if __name__  == '__main__':
    if len(sys.argv) < 3:
        print('\nusage:  python {} [prefix] [threads] \n'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    threads = sys.argv[2]
    print(time.ctime(), 'Map begin')
    MapSortRedupBQSR(prefix,threads).MapPipeline()
    filename = 'MappingHPV/%s.blast8' % prefix
    try:
        trans_col(filename).to_excel(filename.split('.')[0]+ '.result.xlsx', index=False)
    except:
        print(filename + ' is empty！')
