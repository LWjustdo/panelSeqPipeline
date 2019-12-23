#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/1 
@Description :比对，排序，去重，碱基质量控制
'''

import os
import time
import sys
from profile import Profile
import functools

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
    def __init__(self,prefix,threads,bed):
        self.var_path = Profile()
        self.prefix = prefix
        self.threads = threads
        self.bed = bed
        self.path = 'Mapping/InterFiles'


    @timefly
    def mapping(self):
        os.makedirs(self.path,exist_ok=True)
        os.system("{bwa} mem -t {th} -M -Y -R \'@RG\\tID:{T}\\tPL:ILLUMINA\\tLB:{T}\\tSM:{T}\'  "
                         "{refGenome} DataQC/InterFiles/{T}_R1.clean.fastq.gz DataQC/InterFiles/{T}_R2.clean.fastq.gz 1> Mapping/{T}.sam "
                  "2>> Mapping/{T}.map.log".format(T=self.prefix,th=self.threads,**self.var_path))

    @timefly
    def sort(self):
        os.system("{samtools} view -@ {th} -bS Mapping/{T}.sam -o Mapping/{T}.bam >>Mapping/{T}.map.log 2>&1".format(T=self.prefix,th=self.threads,**self.var_path))
        os.system("{samtools} sort -@ {th} Mapping/{T}.bam Mapping/{T}.sort >>Mapping/{T}.map.log 2>&1".format(T=self.prefix,th=self.threads,**self.var_path))

    @timefly
    def redup(self):
        if 'brca' in self.bed:
            inter_bam = "Mapping/{T}.sort.bam".format(T=self.prefix)
        else:
            os.system("{gatk} MarkDuplicates -I Mapping/{T}.sort.bam  -O {p}/{T}.sort.mrkdup.bam --REMOVE_DUPLICATES true "
                  "-M {p}/{T}.MarkDuplicates.xls >> Mapping/{T}.map.log 2>&1".format(T=self.prefix,p=self.path,**self.var_path))
            os.system("{samtools} index {p}/{T}.sort.mrkdup.bam".format(T=self.prefix,p=self.path,**self.var_path))
            inter_bam = "{p}/{T}.sort.mrkdup.bam".format(T=self.prefix,p=self.path)
        return inter_bam

    @timefly
    def bqsr(self):
        inter_bam = self.redup()
        bqsrKnownVcfs = self.var_path['bqsrKnownVcfs'].split(',')
        bvqsrKnownPath = self.var_path['bvqsrKnownPath']
        vcf_args = ' --known-sites '.join([bvqsrKnownPath+i for i in bqsrKnownVcfs])
        os.system("{gatk} BaseRecalibrator -I {ib} -R {refGenome} -L {b} "
                  #"--known-sites {known_dbsnp} --known-sites {known_1000G} "
                    "--known-sites {va} "
                  "-O {p}/{T}.recal.table >>Mapping/{T}.map.log 2>&1".format(ib=inter_bam,va=vcf_args,b=self.bed,T=self.prefix,p=self.path,**self.var_path))
        os.system("{gatk} ApplyBQSR -R {refGenome} -I {ib} -O Mapping/{T}.bqsr.final.bam "
                         "-bqsr {p}/{T}.recal.table >>Mapping/{T}.map.log 2>&1".format(ib=inter_bam,T=self.prefix,p=self.path,**self.var_path))
        os.system("{samtools} index Mapping/{T}.bqsr.final.bam".format(T=self.prefix,**self.var_path))  #MSI5points需要final.bam.bai而非final.bai格式的索引
        os.system("rm Mapping/{T}.sam Mapping/{T}.bam ".format(T=self.prefix,**self.var_path))

    def pipeline(self):
        self.mapping()
        self.sort()
        self.bqsr()


if __name__  == '__main__':
    if len(sys.argv) < 3:
        print('\nusage:  python {} [prefix] [threads] [bed] \n'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    threads = sys.argv[2]
    bed =  sys.argv[3]
    print(time.ctime(), 'MapSortRedupBQSR begin')
    MapSortRedupBQSR(prefix,threads,bed).pipeline()

