#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/30 
@Description :
需要查看（germline的multianno和germline的hgmd）两个文件，确定肿瘤相关致病的基因突变。乳腺癌还需要看BRCA的注释文件;
需要对应肿瘤的突变，如果突变不是关于肿瘤的，就过滤掉。即疾病名称（multianno中的CLNDN列或HGMD中PHEN=）中需要包含下列关键词之一：cancer、tumor、neoplasm，carcinoma;
仅列出致病风险为（pathogenic、likely pathogenic，DM）3种
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

# 变异类型英文转中文
def trans_varient_type(str):
    dic = {}
    with open(var_path['varientType'], encoding='utf-8') as varientType:
        for line in varientType:
            lin = line.strip().split('\t')
            dic[lin[0]] = lin[1]
    if str in dic:
        str = dic[str]
    return str

#germli_hgmd结果添加到表格
def process_hgmd(prefix):
    hgmd_file = 'Result/{}.hgmd_raw_result.xls'.format(prefix)
    if not os.path.exists(hgmd_file):
        print('hgmd_raw_result文件不存在！')
        sys.exit(1)
    # germli_tumor.writerow('基因,相关疾病,突变类型,cDNA改变,氨基酸改变,基因型,致病风险'.split(','))
    germli_tumor = open('Result/{}.germli_tumor_result.xls'.format(prefix), 'w', encoding='utf-8')
    with open(hgmd_file) as hgmd:
        for hgm in hgmd:
            if not hgm.startswith('Chr'):
                hg = hgm.strip().split('\t')
                gene = hg[5]
                cdna, AA = hg[7:9]
                pathogenic, disease = hg[10:12]
                genetype=hg[-1]
                if pathogenic == 'DM' and any(['cancer' in disease, 'tumor' in disease, 'neoplasm' in disease, 'carcinoma' in disease]):
                    info = [gene,disease,'-',cdna,AA,genetype,pathogenic]
                    germli_tumor.write('\t'.join(info)+'\n')

#germli_multi结果添加到表格,需要先运行process_hgmd
def process_multi(prefix):
    germ_mult_file = 'Result/{}.raw_result.xls'.format(prefix)
    if not os.path.exists(germ_mult_file):
        print('multi_raw_result文件不存在！')
        sys.exit(1)
    germli_tumor = open('Result/{}.germli_tumor_result.xls'.format(prefix), 'a', encoding='utf-8')
    with open(germ_mult_file) as multi:
        for mul in multi:
            if not (mul.startswith('TMB') or mul.startswith('Chr')):
                mu = mul.strip().split('\t')
                gene = mu[5]
                type = trans_varient_type(mu[6])
                cdna,AA = mu[7].split('|')[1:3]
                disease = mu[10]
                pathogenicity = mu[11:13]
                clnsig,interVar = pathogenicity
                genetype = mu[-1]
                if (clnsig or interVar) in ['pathogenic', 'likely pathogenic'] and any(['cancer' in disease, 'tumor' in disease, 'neoplasm' in disease]):
                    info = [gene,disease,type,cdna,AA,genetype,'|'.join(pathogenicity)]
                    # germli_tumor.writerow(info)
                    germli_tumor.write('\t'.join(info) + '\n')

def process_brca(prefix):
    brca_file = 'Result/{}.brca_raw_result.xls'.format(prefix)
    if not os.path.exists(brca_file):
        print('brca_raw_result文件不存在！')
        sys.exit(1)
    # germli_tumor = open('Result/{}.germli_tumor_result.csv'.format(prefix), 'a', encoding='gb18030')
    # germli_tumor = csv.writer(germli_tumor)
    germli_tumor = open('Result/{}.germli_tumor_result.xls'.format(prefix), 'a', encoding='utf-8')
    with open(brca_file) as brcaf:
        for brca in brcaf:
            if not brca.startswith('Chr'):
                br = brca.strip().split('\t')
                gene = br[5]
                cdna, AA,pathogenic = br[7:10]
                genetype = br[-1]
                if pathogenic == 'Pathogenic':
                    info = [gene, 'breast cancer', '-', cdna, AA, genetype, pathogenic]
                    # germli_tumor.writerow(info)
                    germli_tumor.write('\t'.join(info) + '\n')

@timefly
def Germline_tumor_pipe(prefix,panel):
    process_hgmd(prefix)
    process_multi(prefix)
    if panel == 'brca':
        process_brca(prefix)


if __name__ == '__main__':
    if len(sys.argv) <3:
        print('\nUsage: python {} [Tprefix] [panel]'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    panel = sys.argv[2]
    print(time.ctime(),'germli_tumor begin...')
    process_hgmd(prefix)
    process_multi(prefix)
    if panel == 'brca':
        process_brca(prefix)
