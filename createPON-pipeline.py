#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/11/7 
@Description :利用正常样本建立PON
'''
from multiprocessing import Process,Pool
from profile import Profile
import functools
import time
import os,sys

var_path = Profile()

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

@timefly
def QC(prefix, threads):
    path = 'DataQC/InterFiles'
    if not os.path.exists(path):
        os.makedirs(path)
    fastq = [ i for i in os.listdir() if prefix in i]
    fastq1 = ''.join([i for i in fastq if '1.f' in i])
    fastq2 = ''.join([i for i in fastq if '2.f' in i])
    os.system("{fastp} -w {th} --in1 {f1} --out1 {p}/{T}_R1.clean.fastq.gz "
              "--in2 {f2} --out2 {p}/{T}_R2.clean.fastq.gz "
              "--low_complexity_filter --correction --length_required=70 "
              "--html DataQC/{T}.QCReport.html --json {p}/{T}.json --report_title {p}/{T}.QCReport "
              " >{p}/{T}.fastp.log 2>&1".format(T=prefix, p=path,f1=fastq1, f2=fastq2, th=threads,**var_path))


class MapSortRedupBQSR():
    def __init__(self,prefix,threads):
        self.var_path = Profile()
        self.prefix = prefix
        self.threads = threads
        self.path = 'Mapping/InterFiles'

    def mapping(self):
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        os.system("{bwa} mem -t {th} -M -Y -R \'@RG\\tID:{T}\\tPL:ILLUMINA\\tLB:{T}\\tSM:{T}\'  "
                         "{refGenome} DataQC/InterFiles/{T}_R1.clean.fastq.gz DataQC/InterFiles/{T}_R2.clean.fastq.gz 1> Mapping/{T}.sam "
                  "2>> Mapping/{T}.map.log".format(T=self.prefix,th=self.threads,**self.var_path))

    def sort(self):
        os.system("{samtools} view -@ {th} -bS Mapping/{T}.sam -o Mapping/{T}.bam >>Mapping/{T}.map.log 2>&1".format(T=self.prefix,th=self.threads,**self.var_path))
        os.system("{samtools} sort -@ {th} Mapping/{T}.bam Mapping/{T}.sort >>Mapping/{T}.map.log 2>&1".format(T=self.prefix,th=self.threads,**self.var_path))

    def redup(self):
        os.system("{gatk} MarkDuplicates -I Mapping/{T}.sort.bam  -O {p}/{T}.sort.mrkdup.bam --REMOVE_DUPLICATES true "
              "-M {p}/{T}.MarkDuplicates.xls >> Mapping/{T}.map.log 2>&1".format(T=self.prefix,p=self.path,**self.var_path))
        inter_bam = "{p}/{T}.sort.mrkdup.bam".format(T=self.prefix,p=self.path)
        return inter_bam

    def bqsr(self):
        inter_bam = self.redup()
        # inter_bam = 'Mapping/InterFiles/{}.sort.mrkdup.bam'.format(self.prefix)
        os.system("{gatk} BaseRecalibrator -I {ib} -R {refGenome} "
                  "--known-sites {known_dbsnp} --known-sites {known_1000G} "
                  "-O {p}/{T}.recal.table >>Mapping/{T}.map.log 2>&1".format(ib=inter_bam,T=self.prefix,p=self.path,**self.var_path))
        os.system("{gatk} ApplyBQSR -R {refGenome} -I {ib} -O Mapping/{T}.bqsr.final.bam "
                         "-bqsr {p}/{T}.recal.table >>Mapping/{T}.map.log 2>&1".format(ib=inter_bam,T=self.prefix,p=self.path,**self.var_path))
        os.system("samtools index Mapping/{T}.bqsr.final.bam".format(T=self.prefix))  #MSI5points需要final.bam.bai而非final.bai格式的索引
        os.system("rm Mapping/{T}.sam Mapping/{T}.bam ".format(T=self.prefix,**self.var_path))

    @timefly
    def mapPipeline(self):
        self.mapping()
        self.sort()
        self.bqsr()

@timefly
def PonOfNormal(prefix,  panel):
    nbam = os.popen("ls Mapping/{0}*final.bam".format(prefix)).read().strip()
    path = 'VarientCalling/InterFiles_call'
    if not os.path.exists(path):
        os.makedirs(path)
    pon_path = '{0}_{1}'.format(var_path['PonPath'], panel)  # 不同panel基线目录
    if not os.path.exists(pon_path):
        os.makedirs(pon_path)
    os.system("{gatk} Mutect2 -R {refGenome} -I {nb} "
              "--max-mnp-distance 0 -O {pp}/{N}_for_pon.vcf.gz "
              ">>VarientCalling/{N}.call.log 2>&1".format(nb=nbam, N=prefix,  pp=pon_path, **var_path))


def basicAnalyze(prefix,panel,threads=16):
    print('{0} 开始对样本{1}建立PON......'.format(time.ctime(), prefix))
    QC(prefix, threads)
    MapSortRedupBQSR(prefix, threads).mapPipeline()
    PonOfNormal(prefix, panel)
    print('{0} 样本{1}建立PON已完成......'.format(time.ctime(), prefix))


if __name__ == '__main__':
    if len(sys.argv) <2:
        print("\n:Usage: python  {} [panel]  [prefix1,prefix2...]\n".format(sys.argv[0]))
        sys.exit(1)
    panel = sys.argv[1]
    prefix_list = sys.argv[2:]
    proc_list =[]
    for prefix in prefix_list:
        proc = Process(target=basicAnalyze,args=(prefix,panel,))
        proc_list.append(proc)
    for p in proc_list:
        p.start()
    for p in proc_list:
        p.join()
    print()
