#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/5 
@Description :
'''

import sys,os
import time
from profile import Profile
import functools

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
#只根据正常样本建立基线,在当前目录执行
def NormaBaselineCNV(bed):
    os.system("{cnvkit} target {b} --annotate {refFlat} --split -o {b}.target >baseline.log 2>&1".format(b=bed,**var_path))
    os.system("{cnvkit} antitarget {b} -g {access_bed} -o {b}.antitarget >>baseline.log 2>&1".format(b=bed,**var_path))
    os.system("{cnvkit} batch -n *bam --output-reference normal.cnn -t  {b}.target -a {b}.antitarget "
              "-f {refGenome} -g {access_bed} >>baseline.log 2>&1".format(b=bed,**var_path))

@timefly
#根据已建立基线分析肿瘤样本
def TumorOnlyCNV(Tprefix,bed):
    tbam = os.popen("ls Mapping/{0}*final.bam".format(Tprefix)).read().strip()
    path = 'CNVCalling/InterFiles_cnvkit'
    if not os.path.exists(path):
        os.makedirs(path)
    if 'exon' in bed:
        os.system("{cnvkit} batch {tb} -r {Agilent_V6_reference} -d {p} >CNVCalling/{t}.cnvkit.log 2>&1".format(tb=tbam,p=path,t=Tprefix,**var_path))
    cns = tbam.split('/')[-1].split('.bam')[0] + '.cns'
    os.system("{cnvkit} call {p}/{c} -o CNVCalling/{c}.result.txt >>CNVCalling/{t}.cnvkit.log 2>&1".format(c=cns,p=path, t=Tprefix,**var_path))

@timefly
#配对样本分析
def Cnvkit(Tprefix,Nprefix,bed):
    TNprefix = Tprefix + '_vs_' + Nprefix + '.somatic'
    tbam = os.popen("ls Mapping/{0}*final.bam".format(Tprefix)).read().strip()
    nbam = os.popen("ls Mapping/{0}*final.bam".format(Nprefix)).read().strip()
    path = 'CNVCalling/InterFiles_cnvkit'
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("{cnvkit} batch {tb} -n {nb} -t {b} --annotate {refFlat} "
              "-f {refGenome} -g {access_bed} --output-reference {p}/{tn}.cnn -d {p} "
              ">CNVCalling/{tn}.cnvkit.log 2>&1".format(tb=tbam,nb=nbam,b=bed,p=path, tn=TNprefix,**var_path))
    cns = tbam.split('/')[-1].split('.bam')[0]+'.cns'
    os.system("{cnvkit} call {p}/{c} -o CNVCalling/{tn}.CNV_cnvkit_result.xls >>CNVCalling/{tn}.cnvkit.log 2>&1".format(c=cns,p=path, tn=TNprefix,**var_path))
    os.system("cp CNVCalling/{}.CNV_cnvkit_result.xls Result/".format(TNprefix))

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('\nUsage: python {} [Tprefix] [Nprefix] [bed]'.format(sys.argv[0]))
        sys.exit(1)
    Tprefix =  sys.argv[1]
    Nprefix =  sys.argv[2]
    bed = sys.argv[3]
    print(time.ctime(),'cnv begin...')
    # NormaBaselineCNV(bed)
    # TumorOnlyCNV(Tprefix, bed)
    Cnvkit(Tprefix,Nprefix,bed)
