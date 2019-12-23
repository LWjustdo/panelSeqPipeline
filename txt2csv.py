#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/9/2 
@Description :
txt2csv：将任意后缀的文本文件转换为csv文件；
csv2txt：将csv文件转换为文本文件；
'''

import csv
import sys
import os
import  time

def mkPath(result_path):
    if not os.path.exists(result_path):
        os.mkdir(result_path)

def txt2csv(txt_file,result_path):
    out_csv = open(txt_file.replace('Result',result_path).replace(file.split('.')[-1],'csv'), 'w', encoding='gb18030')
    out_csv = csv.writer(out_csv)
    with open(txt_file) as f:
        for line in f:
            lin = line.strip().split('\t')
            out_csv.writerow(lin)

def csv2txt(csv_file):
    out = open(file.replace('csv','txt'),'w',encoding='utf-8')
    with open(csv_file, encoding='gb18030') as csv_input:
        row = csv.reader(csv_input)
        for i in row:
            out.write('\t'.join(i)+'\n')
    out.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('\nUsage: {} [file]\n'.format(sys.argv[0]))
        sys.exit(1)
    file = sys.argv[1]
    #str_time = time.strftime('%Y%m%d-%H%M%S', time.localtime())
    result_path = 'Result-CSV'
    mkPath(result_path)
    txt2csv(file,result_path)
    # csv2txt(file)