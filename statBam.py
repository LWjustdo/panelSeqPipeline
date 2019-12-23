#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/1 
@Description :统计bam文件中测序深度，覆盖度，捕获效率等
'''

import time
import os
import sys
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
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
def StatBam(prefix,bed):
    path = 'StateMap/{0}/'.format(prefix)
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("{bamdst} -p {b} -q 30 -o {p} Mapping/{T}.sort.bam".format(p=path, T=prefix, b=bed, **var_path))
    # 重命名加上prefix
    for i in os.listdir(path):
        old_name = path + i
        new_name = path + prefix + '_sort_' + i
        os.rename(old_name, new_name)
    Stat2file(path + prefix + '_sort_coverage.report', prefix + '-sorted.bam',prefix)
    Draw_fig(prefix)
    os.system("{bamdst} -p {b} -q 30 -o {p} Mapping/{T}.bqsr.final.bam".format(p=path, T=prefix, b=bed,**var_path))
    Stat2file(path + 'coverage.report', prefix + '-final.bam',prefix)
    os.system("cp StateMap/stat_bam.xls Result/")

#根据bamdst结果生成统计文件
def Stat2file(report_file, sample,pre):
    # 提取统计信息到指定文件
    qc_file = 'DataQC/{}.QCsummary.xls'.format(pre)
    insert_peak = open(qc_file).readlines()[-1].split('\t')[-2] #插入片段长度
    report_dic = {}
    out_report = open('StateMap/stat_bam.xls', 'a')
    if open('StateMap/stat_bam.xls').read() == '':
        out_report.write('sample\tall_reads\tmapped_rads\tq30\ttarget\tAverage_depth\tCoverage_0x\tCoverage_30x\tbed_length\tinsert_peak\n')
    with open(report_file)  as report:
        for repo in report:
            rep = repo.strip().split('\t')
            if len(rep) == 2:
                report_dic[rep[0]] = rep[1]
    all_reads = report_dic['[Total] Raw Reads (All reads)']
    mapped_rads = report_dic['[Total] Fraction of Mapped Reads']
    q30 = report_dic['[Total] Fraction of MapQ reads in all reads']
    target = report_dic['[Target] Fraction of Target Reads in all reads']
    Average_depth = report_dic['[Target] Average depth']
    Coverage_0x = report_dic['[Target] Coverage (>0x)']
    Coverage_30x = report_dic['[Target] Coverage (>=30x)']
    bed_length = report_dic['[Target] Len of region']
    out_report.write(
        sample + '\t' + all_reads + '\t' + mapped_rads + '\t' + q30 + '\t' + target + '\t' + Average_depth + '\t' + Coverage_0x + '\t' + Coverage_30x + '\t' + bed_length +'\t'+insert_peak+ '\n')
    out_report.close()


#画图
def Draw_fig(prefix):
    # 深度作图
    depth_file = 'StateMap/{0}/{0}_sort_depth_distribution.plot'.format(prefix)
    os.system("sed '1i Depth\tNumber\tFraction\tcum_Number\tcum_Fraction' {0} > {0}2".format(depth_file))
    df = pd.read_csv(depth_file + '2', sep='\t', index_col=['Depth'])
    fig = plt.figure(figsize=(20, 6))
    ax1 = fig.add_subplot(1, 2, 1)  # 创建子图1
    ax1.set_xlabel('sequence depth')
    ax1.set_ylabel('Fraction of bashes(%)')
    plt.sca(ax1)  # 选择子图1
    plt.plot(df['Fraction'])
    ax2 = fig.add_subplot(1, 2, 2)  # 创建子图2
    ax2.set_xlabel('cumulative sequence depth')
    ax2.set_ylabel('Fraction of bashes(%)')
    plt.sca(ax2)
    plt.plot(df['cum_Fraction'])
    fig.savefig(depth_file.split('depth')[0] + 'depth.png')
    # 覆盖度作图
    cover_file = 'StateMap/{0}/{0}_sort_chromosomes.report'.format(prefix)
    df2 = pd.read_csv(cover_file, sep='\t')
    x = [str(i).strip() for i in df2['#Chromosome']]
    y1 = [i for i in df2[' Avg depth']]
    y2 = [i for i in df2[' Coverage%']]
    fig2 = plt.figure(figsize=(15, 8))  # 创建图形2
    chr = fig2.add_subplot(111)  # 创建子图
    chr.bar(x, y1, align='center')  # 柱状图
    chr.set_ylabel('Mean depth')
    chr2 = chr.twinx()  # 设置第2纵坐标，共享x轴
    chr2.plot(x, y2, c='r', lw=1, marker='o', mec='r', mfc='w')  # 折线图
    chr2.set_ylim([0, 105])
    chr2.set_yticks(np.arange(0, 101, 20))
    chr2.set_ylabel('Coverage%')
    fig2.savefig(cover_file.split('chromosome')[0] + 'coverage.png')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('\n Usage: python {} [prefix] [bed]\n'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    bed = sys.argv[2]
    print(time.ctime(),'StatBam begin...')
    StatBam(prefix, bed)
