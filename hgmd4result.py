#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/30 
@Description :
生成hgmd的raw_result文件和hgmd_pro文件
'''

import os
import sys
import csv
import re
from profile import Profile
import functools
import time

var_parh = Profile()

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

# 1字母AA转换为3字母
def transAA(pro):
    AAdict = {'G': 'Gly', 'A': 'Ala', 'V': 'Val', 'L': 'Leu', 'I': 'Ile', 'P': 'Pro', 'F': 'Phe',
               'Y': 'Tyr', 'W': 'Trp', 'S': 'Ser', 'T': 'Thr', 'C': 'Cys', 'M': 'Met', 'N': 'Asn',
               'Q': 'Gln', 'D': 'Asp', 'E': 'Glu', 'K': 'Lys', 'R': 'Arg', 'H': 'His'}
    pro2 = re.findall('[A-Z]+', pro)
    for p in pro2:
        if p in AAdict:
            pro = re.sub(p, AAdict[p], pro)
    return pro

#hgmd生成raw_result文件
def hgmd2raw(prefix):
    if not os.path.exists('Result'):
        os.mkdir('Result')
    hgmd_file = "VarientCalling/{0}.hg19_hgmd.txt".format(prefix)
    if os.path.exists(hgmd_file):
        os.system("cp {0} Result".format(hgmd_file))
        out_raw = open('Result/{}.hgmd_raw_result.xls'.format(prefix), 'w')
        out_raw.write('Chr\tStart\tEnd\tRef\tAlt\tGene\tTrans\tCdna\tAA\tDBsnp\tClass\tPhen\tVAF\tDP\tgenetype\n')
        with open(hgmd_file) as f:
            for line in f:
                lin = line.strip().split('\t')
                # 注释结果拆分
                classes = re.search('CLASS=(.*?);', lin[1]).group(1)
                gene = re.search('GENE=(.*?);', lin[1]).group(1)
                dna = re.search('DNA=(.*?);', lin[1])
                if dna == None:
                    trans, cdna = '.', '.'
                else:
                    trans, cdna = dna.group(1).split(':')
                prot = re.search('PROT=(.*?);', lin[1])
                if prot == None:
                    pro = '.'
                else:
                    pro = prot.group(1).split(':')[-1]
                    pro = transAA(pro)
                db = re.search('DB=(.*?);', lin[1])
                if db == None:
                    rs = '.'
                else:
                    rs = db.group(1)
                phen = re.search('PHEN="(.*?)"', lin[1]).group(1)
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
                    info = '\t'.join(lin[2:7]) + '\t' + gene + '\t' + trans + '\t' + cdna + '\t' + pro + '\t' + rs + '\t' + classes + '\t' + phen + '\t' + vaf_str + '%\t' + dp + '\t' + gt + '\n'
                    out_raw.write(info)
        out_raw.close()

# hgmd_pro注释
def hgmd_pro(prefix):
    hgmd_file = "VarientCalling/{0}.hg19_hgmd.txt".format(prefix)
    out_pro = open('Result/{}.hgmd_pro_result.xls'.format(prefix), 'w')
    out_pro.write('chromosome\tstartCoord\tendCoord\tgene\tdescr\trefseq\thgvs\tdisease\tchrom'
                  '\tgenename\tgdbid\tomimid\tamino\tdeletion\tinsertion\tcodon\tcodonAff\thgvsAll'
                  '\tdbsnp\tinheritance\tgnomad_AC\tgnomad_AF\tgnomad_AN\ttag\tdmsupport\trankscore'
                  '\tmutype\tauthor\ttitle\tfullname\tallname\tvol\tpage\tyear\tpmid\tpmidAll\treftag\tcomments'
                  '\tacc_num\tnew_date\tbase\tclinvarID\tclinvar_clnsig\n')
    with open(hgmd_file) as f2:
        for line2 in f2:
            lin2 = line2.strip().split('\t')
            gene = re.search('GENE=(.*?);', lin2[1]).group(1)
            rs = re.search('DB=(.*?);', lin2[1]).group(1) if re.search('DB=(.*?);', lin2[1]) else '-'
            with open(var_parh['hgmd_pro_file']) as pro:
                for pr in pro:
                    p = pr.strip().split('\t')
                    if lin2[2] == p[0] and lin2[3] == p[1] and lin2[4] == p[2] and gene == p[3] and rs == p[18]:  # 染色体位置和基因名称
                        out_pro.write(pr)
    out_pro.close()

@timefly
def Hgmd_pipe(prefix):
    hgmd2raw(prefix)
    hgmd_pro(prefix)

if __name__ == '__main__':
    if len(sys.argv) <2:
        print('\nUsage: python {} [prefix]'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    print(time.ctime(),'hgmd2raw begin...')
    hgmd2raw(prefix)
    hgmd_pro(prefix)