#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/6/24 
@Description :
'''
import os
import time
import sys
import functools
from profile import Profile
import pandas as pd

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
def CNV(Tprefix,Nprefix,threads,bed):
    TNprefix = Tprefix + '_vs_' + Nprefix + '.somatic'
    tbam = os.popen("ls Mapping/{0}*final.bam".format(Tprefix)).read().strip()
    nbam = os.popen("ls Mapping/{0}*final.bam".format(Nprefix)).read().strip()
    path = 'CNVCalling/InterFiles_freec'
    if not os.path.exists(path):
        os.makedirs(path)
    #新建配置文件
    configFile = open(var_path['CNVconfig_file'], 'r').read()
    newConfigFile = configFile.replace('$outDir', path).replace('$Threads', str(threads)).replace('$TumorBam', tbam).replace('$NormalBam', nbam).replace('$bed', bed)
    out = open('CNVCalling/InterFiles_freec/' +TNprefix + '_FREEC_config.txt', 'w')
    out.write(newConfigFile)
    out.close()
    os.system("{freec_path}/src/freec -conf {p}/{tn}_FREEC_config.txt >CNVCalling/{tn}.freec.log 2>&1".format(p=path,tn=TNprefix,**var_path))
    os.chdir(path)
    os.system("cat {freec_path}/scripts/assess_significance.R|R --slave --args {0}.bqsr.final.bam_CNVs {0}.bqsr.final.bam_ratio.txt "
              "2>> ../{1}.freec.log".format(Tprefix, TNprefix,**var_path))
    os.system("cat {freec_path}/scripts/makeGraph.R|R --slave --args 2 {0}.bqsr.final.bam_ratio.txt 2>> ../{1}.freec.log".format(Tprefix,TNprefix,**var_path))
    pvalue = '{0}.bqsr.final.bam_CNVs.p.value.txt'.format(Tprefix)
    df = pd.read_csv(pvalue, sep='\t')
    df2 = df[df.WilcoxonRankSumTestPvalue < 0.01]
    df3 = df2[df2.KolmogorovSmirnovPvalue < 0.01]
    df3.to_csv(pvalue + '2', sep='\t',index=None)
    os.system("{bedtools} intersect -a {pv} -b {genePosition} -wa -wb | {bedtools} groupby -i - -g 1,2,3,4,5 -c 11 -o collapse "
              "> ../{tn}.CNV_freec_result.xls ".format (pv=pvalue+'2',tn=TNprefix,**var_path))
    os.chdir('../../')
    os.system("cp CNVCalling/{}.CNV_freec_result.xls Result/".format(TNprefix))

if __name__  == '__main__':
    if len(sys.argv) < 5:
        print('\nusage:  python {} [Tprefix] [Nprefix] [threads] [bed] \n'.format(sys.argv[0]))
        sys.exit(1)
    Tprefix = sys.argv[1]
    Nprefix = sys.argv[2]
    threads = sys.argv[3]
    bed = sys.argv[4]
    print(time.ctime(),'CNV begin...')
    CNV(Tprefix, Nprefix, threads, bed)