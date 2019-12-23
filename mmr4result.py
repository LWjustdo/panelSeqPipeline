#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/26 
@Description :此处Tprefix,Nprefix不要加germline/somatic。针对配对分析样本。
判断错配修复基因是正常（pMMR)还是缺陷（dMMR);
需要三个注释文件（somatic的multianno，germline的multianno，germline的hgmd）的raw_result结果文件存在；
'MLH1','PMS2','MSH2','MSH6' 4个基因中有一个含致病变异，则判断为dMMR;都没有致病变异，则判断为pMMR。
'''

import time
import os,sys
from collections import defaultdict
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

def Process_multi(file):
    mul_dict = defaultdict(list)
    with open(file) as multi:
        for mul in multi:
            if not mul.startswith('TMB'):
                mu = mul.strip().split('\t')
                gene = mu[5] ###! 5/6
                pathogenicity = mu[11:13]  ###! [11:13]/12:14
                clnsig,interVar = pathogenicity
                if gene in ['MLH1','PMS2','MSH2','MSH6'] and 'Pathogenic' in pathogenicity:
                    if clnsig == 'Pathogenic' and interVar == 'Pathogenic':
                        mul_dict[gene].append('Pathogenic')
                    else:
                        mul_dict[gene].append(','.join(pathogenicity))
    return mul_dict

def Process_hgmd(file):
    hgmd_dict={}
    with open(file) as hgmd:
        for hgm in hgmd:
            hg = hgm.strip().split('\t')
            gene = hg[5]
            pathogenicity = hg[10]
            if gene in ['MLH1','PMS2','MSH2','MSH6'] and pathogenicity == 'DM':
                hgmd_dict[gene] = pathogenicity
    return hgmd_dict

@timefly
def MMR(Tprefix,Nprefix):
    hgmd_file = 'Result/{}.germline.hgmd_raw_result.xls'.format(Nprefix)
    germ_mult_file =  'Result/{}.germline.raw_result.xls'.format(Nprefix)
    somat_mult_file = 'Result/{}_vs_{}.somatic.raw_result.xls'.format(Tprefix,Nprefix)
    for raw_file in [hgmd_file,germ_mult_file,somat_mult_file]:
        if not os.path.exists(raw_file):
            print('raw_result文件不存在！')
            sys.exit(1)
    hgmd_dict = Process_hgmd(hgmd_file)
    germ_mult_dict = Process_multi(germ_mult_file)
    somat_mult_dict = Process_multi(somat_mult_file)
    mmr_out = open('Result/{}_vs_{}.somatic.MMR_result.xls'.format(Tprefix,Nprefix),'w',encoding='utf-8')
    if hgmd_dict == germ_mult_dict == somat_mult_dict == {}:
        mmr_out.write('pMMR\n')
    else:
        mmr_out.write('dMMR\n')
        mmr_out.write('hgmd: '+ str(hgmd_dict)+'\n')
        mmr_out.write('germline_multi: '+ str(germ_mult_dict)+'\n')
        mmr_out.write('somatic_multi: '+ str(somat_mult_dict)+'\n')
    mmr_out.close()


if __name__ == '__main__':
    if len(sys.argv) <2:
        print('\nUsage: python {} [Tprefix] [Nprefix]'.format(sys.argv[0]))
        sys.exit(1)
    Tprefix = sys.argv[1]
    Nprefix = sys.argv[2]
    print(time.ctime(),'MMR begin...')
    MMR(Tprefix, Nprefix)

