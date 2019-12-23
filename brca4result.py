#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/8/2 
@Description :
'''

import time
import os,sys
import re
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

@timefly
def BrcaResult(prefix):
    if not os.path.exists('Result'):
        os.mkdir('Result')
    brca_file = "VarientCalling/{0}.hg19_brca.txt".format(prefix)
    if os.path.exists(brca_file):
        out_raw = open( 'Result/{}.brca_raw_result.xls'.format(prefix), 'w')
        out_raw.write('Chr\tStart\tEnd\tRef\tAlt\tGene\tTrans\tCdna\tAA\tClinical_significance\tFre_EAS_ExAc\tVAF\tDP\tgenetype\n')
        with open(brca_file) as f:
            for line in f:
                lin = line.strip().split('\t')
                # 注释结果拆分
                gene = re.search('GENE=(.*?);', lin[1]).group(1)
                trans = re.search('Trans=(.*?);', lin[1]).group(1)
                cdna = re.search('cDNA=(.*?);', lin[1]).group(1)
                pro = re.search('Prot=(.*?);', lin[1]).group(1)
                Clinical_significance = re.search('Clinical_significance_ENIGMA=(.*?);', lin[1]).group(1)
                Fre_EAS = re.search('Fre_EAS_ExAc=(.*)', lin[1]).group(1)
                Fre_EAS_str = '-' if Fre_EAS == '-' else "%.2f" % (float(Fre_EAS) * 100)
                # vaf
                last_line = lin[-1].split(':')  # GT:AD:DP:GQ:PL
                genetype = last_line[0]
                if genetype == '0/0' or genetype == '1/1':
                    gt = '纯和'
                else:
                    gt = '杂合'
                alt_num = last_line[1].split(',')[1]
                dp = last_line[2]
                vaf = int(alt_num) / float(dp)
                vaf_str = "%.2f" % (vaf * 100)
                if vaf >= 0.01:
                    info = '\t'.join(lin[2:7]) + '\t' + gene + '\t' + trans + '\t' + cdna + '\t' + pro + '\t' + Clinical_significance + '\t' + Fre_EAS_str + '%\t' + vaf_str + '%\t' + dp + '\t' + gt + '\n'
                    out_raw.write(info)
        out_raw.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('\nUsage: python {} [prefix]'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    print(time.ctime(), 'BrcaResult begin...')
    BrcaResult(prefix)
