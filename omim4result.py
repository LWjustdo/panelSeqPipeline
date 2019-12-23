#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/3 
@Description :用omimDB对hg19_multianno.txt文件注释，如果基因存在于omimDB文件，则将对应疾病写入hg19_multianno.txt
'''

import time
import os,sys
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
def Omim(prefix):
    multi_file = 'Result/' + prefix + '.hg19_multianno.txt'
    if not os.path.exists('Result'):
        os.mkdir('Result')
    if not os.path.exists(multi_file):
        os.system("cp VarientCalling/{0}.hg19_multianno.txt Result".format(prefix))
     #将omim文件中基因和疾病写入字典备用
    omim_dict = {}
    with open(var_path['omimDB']) as omim:
        for omi in omim:
            om = omi.strip().split('\t')
            omim_dict[om[0]] = om[1]

    out = open(multi_file.split('txt')[0]+'omim_result.xls','w')
    with open(multi_file) as multi:
        for mul in multi:
            mu = mul.strip().split('\t')
            if mu[6] in omim_dict:
                phe = omim_dict[mu[6]]
            else:
                phe = '.'
            out.write(phe+'\t'+mul)
    out.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('\nUsage: python {} [prefix]\n'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    Omim(prefix)