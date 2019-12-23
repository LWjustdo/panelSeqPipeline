#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/17 
@Description :健康位点筛查,输入健康位点列表文件healthRS.txt/bed文件/uncover.bed文件/变异注释文件hg19_multianno.txt
'''
import  os,sys
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
def Health(prefix,bed):
    multi_file = 'Result/' + prefix + '.hg19_multianno.txt'
    if not os.path.exists('Result'):
        os.mkdir('Result')
    if not os.path.exists(multi_file):
        os.system("cp VarientCalling/{0}.hg19_multianno.txt Result".format(prefix))
    #取bed文件与uncover.bed文件差集，得到cover.bed文件
    uncover_bed = 'StateMap/{t}/uncover.bed'.format(t=prefix.split('.')[0])
    cover_bed = os.popen("{bedtools} subtract -a {b} -b {unb} ".format(b=bed,unb=uncover_bed,**var_path)).read().strip().split('\n')
    rs_bed_dic = {}
    # out_tmp = open('rsbedtmp.txt','w')
    for be in cover_bed:
        b = be.split('\t')
        with open(var_path['healthRS']) as rsf:
            for rs in rsf:
                r = rs.strip('\n').split('\t')
                if r[0] == b[0] and r[1] >= b[1] and r[1] <= b[2]:
                    rs_bed_dic[r[4]] = r[2:4] + r[-2:]  # {'rs61816761':['T','C','ST14','皮肤保湿能力']}
                    # out_tmp.write(rs)
    # out_tmp.close()
    out = open(multi_file.split('hg19')[0] + 'health_result.xls', 'w')
    multi_dic = {}
    with open(multi_file) as multi:
        for mult in multi:
            if not mult.startswith('Chr'):
                mul = mult.strip().split('\t')
                ref, alt = mul[3], mul[4]
                rs_mul = mul[11]
                genetype = mul[-1].split(':')[0] #0/1或1/1
                multi_dic[rs_mul] = [ref, alt, genetype]  # {'rs61816761':['A','T','0/1']}
    for k, v in rs_bed_dic.items():
        if k in multi_dic:
            if multi_dic[k][-1] == '1/1':
                base = multi_dic[k][1] + multi_dic[k][1]
            else:
                base = multi_dic[k][0] + multi_dic[k][1]
        else:
            base = v[0]
        info = '\t'.join(v[-2:]) + '\t' + k + '\t' + base + '\n'
        out.write(info)
    out.close()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('\nUsage: python {} [prefix] [bed]\n'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    bed = sys.argv[2]
    Health(prefix,bed)