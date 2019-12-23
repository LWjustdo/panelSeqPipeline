#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/6/24 
@Description :fastq文件质控
'''

import os
import time
import sys
import functools
from profile import Profile

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
    print('{0} 开始对样本{1}进行变异分析......'.format(time.ctime(), prefix))
    path = 'DataQC/InterFiles'
    os.makedirs(path,exist_ok=True)
    fastq = [ i for i in os.listdir() if prefix in i]
    fastq1 = ''.join([i for i in fastq if '1.f' in i])
    fastq2 = ''.join([i for i in fastq if '2.f' in i])
    os.system("{fastp} -w {th} --in1 {f1} --out1 {p}/{T}_R1.clean.fastq.gz "
              "--in2 {f2} --out2 {p}/{T}_R2.clean.fastq.gz "
              "--low_complexity_filter --correction --length_required=70 "
              "--html DataQC/{T}.QCReport.html --json {p}/{T}.json --report_title {p}/{T}.QCReport "
              " >{p}/{T}.fastp.log 2>&1".format(T=prefix, p=path,f1=fastq1, f2=fastq2, th=threads,**var_path))
    os.system("python {summary4fastp} {p}/{T}.json > DataQC/{T}.QCsummary.xls ".format(T=prefix,p=path,**var_path))
    if not os.path.exists('Result'):
        os.mkdir('Result')

if __name__  == '__main__':
    if len(sys.argv) < 2:
        print('\nusage:  python {} [prefix] [threads]\n'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    threads = sys.argv[2]
    QC(prefix, threads)
