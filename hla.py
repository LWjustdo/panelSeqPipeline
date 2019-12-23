#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/6/24 
@Description :HLA分型及预测新抗原
'''

import  os
import time
import sys
import functools
from profile import Profile

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

#HLA分型分析
@timefly
def HLA_type(Tprefix,threads):
    path = 'HLA/InterFiles_hla'
    os.makedirs(path,exist_ok=True)
    os.system("gunzip -c DataQC/InterFiles/{T}_R1.clean.fastq.gz >DataQC/{T}_R1.clean.fastq ".format(T=Tprefix))
    os.system("gunzip -c DataQC/InterFiles/{T}_R2.clean.fastq.gz >DataQC/{T}_R2.clean.fastq ".format(T=Tprefix))
    os.system("bash {hlahd_path}/bin/hlahd.sh -t {th} -m 100 -f {hlahd_path}/freq_data "
              "DataQC/{T}_R1.clean.fastq DataQC/{T}_R2.clean.fastq {hlahd_path}/HLA_gene.split.txt "
              "{hlahd_path}/dictionary {T} {p} >HLA/{T}.hla.log 2>&1".format(T=Tprefix,p=path,th=threads,**var_path))
    os.system("cp {p}/{T}/result/{T}_final.result.txt  {p}/{T}_final.result.txt".format(T=Tprefix,p=path))
    os.system("rm DataQC/InterFiles/{T}_R1.clean.fastq.gz DataQC/InterFiles/{T}_R2.clean.fastq.gz".format(T=Tprefix))

#利用vcf文件和hla分型文件预测新抗原
@timefly
def NeoPred(Tprefix,Nprefix):
    TNprefix = Tprefix + '_vs_' + Nprefix + '.somatic'
    #预处理hla分型文件格式
    path = 'HLA/InterFiles_neo'
    os.makedirs(path,exist_ok=True)
    hla_result = 'HLA/InterFiles_hla/{}_final.result.txt'.format(Tprefix)
    hla_result2 = open('Result/{}.hla_result.xls'.format(Tprefix),'w')
    hla_file = open('{p}/{tn}_hla.txt'.format(p=path,tn=TNprefix),'w') #neo输入文件
    hla_file.write("Patient\tHLA-A_1\tHLA-A_2\tHLA-B_1\tHLA-B_2\tHLA-C_1\tHLA-C_2\n")
    type_dict = {}
    with open(hla_result) as f:
        for line in f:
            lin = line.strip().split('\t')
            type_dict[lin[0]] = lin[1:]
    types = type_dict['A']
    types.extend(type_dict['B'])
    types.extend(type_dict['C'])
    hla_newtype = [types[n-1][4:11] if x=='-' else x[4:11] for n,x in enumerate(types)] #调整分型格式
    hla_result2.write('\t'.join(hla_newtype)+'\n')
    hla_result2.close()
    # 如果出现分型为‘-’，表示与前一个分型一样
    types = [types[n-1] if x=='-' else x for n,x in enumerate(types)]
    types2 = [i.replace('*','_').replace('-','_').replace(':','_') for i in types] #HLA-A*31:01:02 --> HLA_A_31_01_02
    hla_file.write(TNprefix+'\t'+'\t'.join(types2)+'\n')
    hla_file.close()
    #预处理vcf文件格式
    os.system("cp VarientCalling/{tn}.PASS.vcf {p}/".format(tn=TNprefix,p=path))
    vcf_file = "{p}/{tn}.PASS.vcf".format(p=path,tn=TNprefix)
    out_vcf = open("{p}/{tn}.vcf".format(p=path,tn=TNprefix),'w')
    with open(vcf_file) as vcf:
        flag = False
        for line in vcf:
            lin = line.strip().split('\t')
            if line.startswith('##'):
                line = line
            elif line.startswith('#CHROM'):
                if lin[-1] == Tprefix: #判断tumor/normal顺序；如果tumor在后则不需变动，反之则要调换顺序
                    line = line
                    flag = True
                else:
                    line = '\t'.join(lin[:-2])+"\t"+lin[-1]+"\t"+lin[-2]+'\n'
            else:
                if flag:
                    line = line
                else:
                    line = '\t'.join(lin[:-2]) + "\t" + lin[-1] + "\t" + lin[-2] + '\n'
            out_vcf.write(line)
    out_vcf.close()
    os.system("source /home/xueqiang.liu/anaconda3/bin/activate python27 && "
              "python {NeoPredPipe} -I {p} -H {p}/{T}_hla.txt -o {p} -n {T} -c 1 -E 8 9 10 >>HLA/{T}.neo.log 2>&1 && "
              "source /home/xueqiang.liu/anaconda3/bin/deactivate".format(p=path,T=TNprefix,**var_path))

    neo_file = 'HLA/InterFiles_neo/{}.neoantigens.txt'.format(TNprefix)
    neo_indel_file = 'HLA/InterFiles_neo/{}.neoantigens.Indels.txt'.format(TNprefix)
    rawresult_file = 'Result/{}.raw_result.xls'.format(TNprefix)
    if os.path.exists(neo_file) and os.path.exists(neo_indel_file):
        addPro4Neo(neo_file,rawresult_file)
        addPro4Neo(neo_indel_file,rawresult_file)
        ProcessNeo(neo_file+'_pro')
        ProcessNeo(neo_indel_file+'_pro')
        os.system("cat HLA/%s.neoantigens.result.xls HLA/%s.neoantigens.Indels.result.xls|sort -V -k5|awk '{if(NR<=10)print $0}'>Result/%s.neoantigens_result.xls" % (TNprefix,TNprefix,TNprefix))
    else:
        print('neoantigens文件不存在')

#neo结果取Binding Affinity<500 且Rank<0.5强结合的数据
def ProcessNeo(file):
    out_neo = open(file.replace('/InterFiles_neo','').replace('txt_pro','result.xls'),'w')
    # out_neo.write('geneAA\tHLA\tpeptide\tBinding Affinity\tRank\n')
    with open(file) as neofile:
        for neo in neofile:
            ne = neo.strip().split('\t')
            if float(ne[21]) < 500 and float(ne[22]) < 0.5:
                info = [ne[0]]+ne[10:12]+ne[21:23]
                out_neo.write('\t'.join(info)+'\n')
    out_neo.close()

def addPro4Neo(file,rawfile):
    out_pro = open(file+'_pro', 'w')
    with open(file) as neofile:
        for neo in neofile:
            if neo != '':
                ne = neo.strip().split('\t')
                # print(ne)
                with open(rawfile) as rawf:
                    for raw in rawf:
                        ra = raw.strip().split('\t')
                        # print(ra)
                        if len(ra) >4:
                            if all([ne[3] ==ra[0],ne[4]==ra[1],ne[5]==ra[3],ne[6]==ra[4]]):
                                gene = ne[7].split(':')[0]
                                pro = ra[7].split('|')[-1]
                                info = gene+'-'+pro+'\t'+ neo
                                out_pro.write(info)
    out_pro.close()


if __name__  == '__main__':
    if len(sys.argv) < 4:
        print('\nusage:  python {} [Tprefix] [Nprefix] [threads] \n'.format(sys.argv[0]))
        sys.exit(1)
    print(time.ctime(),'HLA begin')
    Tprefix = sys.argv[1]
    Nprefix = sys.argv[2]
    threads = sys.argv[3]
    HLA_type(Tprefix, threads)
    # NeoPred(Tprefix, Nprefix)