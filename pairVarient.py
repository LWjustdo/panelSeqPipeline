#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/2
@Description :配对样本检测变异
'''

import os, sys
import time
from multiprocessing import Process, Pool
from profile import Profile
import functools


def timefly(func):
    @functools.wraps(func)
    def wrapper(*args, **kw):
        s = time.time()
        res = func(*args, **kw)
        e = time.time()
        print('{} runtime: {}'.format(func.__name__, TransTime(e - s)))
        return res

    return wrapper


def TransTime(seconds):
    h = seconds // 3600
    m = seconds % 3600 // 60
    s = seconds % 60
    return '{}h {}min {:.0f}s'.format(h, m, s)


var_path = Profile()


@timefly
# 利用正常样本入基线库，再汇总
def PonOfNormal(Tprefix, Nprefix, bed, panel):
    TNprefix = Tprefix + '_vs_' + Nprefix + '.somatic'
    nbam = os.popen("ls Mapping/{0}*final.bam".format(Nprefix)).read().strip()
    path = 'VarientCalling/InterFiles_call'
    if not os.path.exists(path):
        os.makedirs(path)
    pon_path = '{0}_{1}'.format(var_path['PonPath'], panel)  # 不同panel基线目录
    if not os.path.exists(pon_path):
        os.makedirs(pon_path)
    os.system("{gatk} Mutect2 -R {refGenome} -I {nb} "
              "--max-mnp-distance 0 -O {pp}/{N}_for_pon.vcf.gz "
              ">>VarientCalling/{TN}.call.log 2>&1".format(nb=nbam, N=Nprefix, TN=TNprefix, pp=pon_path, **var_path))
    vcf_files = os.popen(" ls %s/*pon.vcf.gz.tbi |while read i;do echo -V ${i/.tbi/};done" % pon_path).read().replace('\n', ' ')
    os.system("{gatk} GenomicsDBImport -R {refGenome} -L {b} --genomicsdb-workspace-path pon_db_{TN} {v} >>VarientCalling/{TN}.call.log 2>&1"
        .format(b=bed, v=vcf_files, TN=TNprefix, **var_path))
    os.system("{gatk} CreateSomaticPanelOfNormals -R {refGenome} -V gendb://pon_db_{TN} -O {p}/{N}.pon.vcf.gz >>VarientCalling/{TN}.call.log 2>&1"
        .format(v=vcf_files, N=Nprefix, p=path, TN=TNprefix, **var_path))
    os.system("rm -rf pon_db_{}".format(TNprefix))


@timefly
# 利用基线库对肿瘤样本变异检测
def Mutect2(Tprefix, Nprefix, bed, splitNum):
    path = 'VarientCalling/InterFiles_call'
    TNprefix = Tprefix + '_vs_' + Nprefix + '.somatic'
    path2 = '{0}/{1}'.format(path, TNprefix)
    tbam = os.popen("ls Mapping/{0}*final.bam".format(Tprefix)).read().strip()
    nbam = os.popen("ls Mapping/{0}*final.bam".format(Nprefix)).read().strip()
    os.system("{gatk} SplitIntervals -R {refGenome} -L {b} -scatter {num} -O {p2}_bedpath  >>VarientCalling/{TN}.call.log 2>&1"
              .format(TN=TNprefix, b=bed, p2=path2,num=splitNum, **var_path))
    bedpath = path2 + '_bedpath'
    vcfpath = path2 + '_vcfpath'
    if not os.path.exists(vcfpath):
        os.makedirs(vcfpath)
    proc_list = []
    for bed in os.listdir(bedpath):
        pre = bed.split('-')[0]
        bed2 = os.path.join(bedpath, bed)
        proc = Process(target=Split4Mutect2, args=(path,Nprefix,Tprefix, TNprefix, nbam,tbam, bed2, pre,))
        proc_list.append(proc)
    for p in proc_list:
        p.start()
    for p in proc_list:
        p.join()
    vcflist = [os.path.join(vcfpath,vcf) for vcf in os.listdir(vcfpath) if vcf.endswith('vcf.gz')]
    statslist = [os.path.join(vcfpath,vcf) for vcf in os.listdir(vcfpath) if vcf.endswith('stats')]
    vcflist2 = ' -I '.join(vcflist)
    statslist2 = ' -stats '.join(statslist)
    os.system("{gatk} MergeVcfs -I {vl} -O  {p}/{TN}.raw.vcf.gz  >>VarientCalling/{TN}.call.log 2>&1".format(TN=TNprefix, p=path,vl=vcflist2, **var_path))
    os.system("{gatk} MergeMutectStats -stats {stats} -O  {p}/{TN}.stats  >>VarientCalling/{TN}.call.log 2>&1".format(TN=TNprefix, p=path,stats=statslist2, **var_path))
    os.system("{gatk} FilterMutectCalls -V {p}/{TN}.raw.vcf.gz -O {p}/{TN}.vcf.gz -R {refGenome} -stats {p}/{TN}.stats "
              "--max-events-in-region 5 --max-alt-allele-count 5 --min-slippage-length 10 --long-indel-length 10 "
              ">>VarientCalling/{TN}.call.log 2>&1".format(TN=TNprefix, p=path, **var_path))
    os.system("{gatk} SelectVariants -V {p}/{TN}.vcf.gz -O VarientCalling/{TN}.PASS.vcf "
              "--exclude-filtered true >>VarientCalling/{TN}.call.log 2>&1".format(TN=TNprefix, p=path, **var_path))


def Split4Mutect2(path,Nprefix, Tprefix, TNprefix, nbam, tbam, bed, pre):
    vcfpath = '{0}/{1}_vcfpath'.format(path, TNprefix)
    os.system("{gatk} Mutect2 -R {refGenome} -I {nb} -normal {N} -I {tb} -L {b} --germline-resource {af-only-gnomad} --max-reads-per-alignment-start 0 "
        "--af-of-alleles-not-in-resource 0.0000025 -pon {p}/{N}.pon.vcf.gz -O {vp}/{n}.raw.vcf.gz >>VarientCalling/{TN}.call.log 2>&1"
        .format(nb=nbam, tb=tbam, T=Tprefix, p=path, vp=vcfpath, N=Nprefix, TN=TNprefix, b=bed, n=pre, **var_path)) #--max-reads-per-alignment-start 0


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('\nusage:  python {} [Tprefix] [Nprefix] [bed] [panel] [splitNum]\n'.format(sys.argv[0]))
        sys.exit(1)
    Tprefix = sys.argv[1]
    Nprefix = sys.argv[2]
    bed = sys.argv[3]
    panel = sys.argv[4]
    splitNum = sys.argv[5]
    PonOfNormal(Tprefix, Nprefix, bed, panel)
    Mutect2(Tprefix, Nprefix, bed, splitNum)
    print('Mutect2 end...')
