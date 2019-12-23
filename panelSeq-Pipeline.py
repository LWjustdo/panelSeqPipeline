#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/4/28 
@Description :本流程可处理single/paired数据
'''

from multiprocessing import Process,Pool
from argparse import ArgumentParser
import time
import sys,os
from qc import QC
from map import MapSortRedupBQSR
from statBam import StatBam
from singleVarient import SingleNormalVarient,SingleTumorVarient
from pairVarient import PonOfNormal,Mutect2
from cnv import CNV
from fusion import Fusion
from msi import MSI
from hla import HLA_type,NeoPred
from vcfAnnotate import VcfAnnotation
from brca4result import BrcaResult
from hgmd4result import Hgmd_pipe
from multi_drug4result import Multi2Raw,MultiDrug
from germ_tumor4result import Germline_tumor_pipe
from omim4result import Omim
from chemo4result import Chemo_drug_pipe
from cnv_cnvkit import Cnvkit,TumorOnlyCNV
from msi_5points import MSI5Points
from mmr4result import MMR
import requests,shutil

def TransTime(seconds):
    h = seconds//3600
    m = seconds%3600//60
    s = seconds%3600%60
    return '{}h {}min {:.0f}s'.format(h,m,s)

def Parser_opt():
    parser = ArgumentParser()
    parser.add_argument('--sample_num', dest='S_num', type=str, default='Single', help='单样本分析还是配对分析，[Single,Paired]，默认Single')
    parser.add_argument('--panel',dest='panel',type=str,default='',help='输入panel名称：[560gene,24gene,exon,brca]，必须！')
    parser.add_argument('--sample_type',dest='S_type',type=str,default='normal',help='样本类型，Single分析时需说明，[normal,tumor],默认normal！')
    parser.add_argument('--tumor_prefix', dest='T_prefix', type=str, default='', help='输入肿瘤样本名称，需与fastq文件前缀一致！必须！')
    parser.add_argument('--normal_prefix', dest='N_prefix', type=str, default='', help='输入正常样本名称，需与fastq文件前缀一致！Paired分析时必须！')
    parser.add_argument('--cancer',dest='cancer',type=str,default='Solid',help='输入患者癌症名称,'
                                        '[NSCLC(非小细胞肺癌),Breast(乳腺癌),Colorectal(结直肠癌),Gastric_Cancer(胃癌),Esophagus_Cancer(食管癌),Gastrointestinal_Stromal(胃肠道间质瘤)'
                                        'Thyroid_Cancer(甲状腺癌)]，默认实体瘤，必须！')
    parser.add_argument('--threads', dest='threads', type=int, default=16, help='输入配置cpu线程数，默认16')
    parser.add_argument('--cnv', dest='cnv', type=bool, default=True, help='是否进行CNV分析，适用于Paired及WES分析，[False,True]，默认True')
    parser.add_argument('--msi', dest='msi', type=bool, default=True, help='是否进行MSI分析，只适用于Paired分析，[False,True]，默认True')
    parser.add_argument('--fusion', dest='fusion', type=bool, default=True, help='是否进行Fusion分析，[False,True]，默认True')
    return parser

def basicAnalyze(T_prefix,threads, bed):
    MapSortRedupBQSR(T_prefix, threads, bed).pipeline()
    StatBam(T_prefix, bed)

def mkdir(path):
    # 去除首位空格
    path=path.strip()
    # 去除尾部 \ 符号
    path=path.rstrip("\\")
 
    # 判断路径是否存在
    # 存在     True
    # 不存在   False
    isExists=os.path.exists(path)
 
    # 判断结果
    if not isExists:
        os.makedirs(path) 
        return True
    else:
        return False


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('\nUsage:        python {} -h \n'.format(sys.argv[0]))
        sys.exit(1)
    parser = Parser_opt()
    args = parser.parse_args()

    if args.panel == '':
        print('\n请注明panel名称！\n')
        sys.exit(1)
    bed_path = '/home/longzhao/panelSel/bed'
    bed_list = os.listdir(bed_path)
    bed = [os.path.join(bed_path, i) for i in bed_list if args.panel in i][0]
    t0 = time.time()
    print('\n'+'#'*60)
    print('python ' + ' '.join(sys.argv))
    if args.S_num == 'Single': #单样本分析
        if args.S_type == 'tumor' :#tumor样本
            QC(args.T_prefix, args.threads)
            basicAnalyze(args.T_prefix, args.threads, bed)
            if args.fusion: #是否分析fusion
                Fusion(args.T_prefix,bed)
            SingleTumorVarient(args.T_prefix, bed,args.panel,args.threads)
            # 样本类型是normal时，文件名称加somatic；否则文件名称加germline
            base = '.somatic'
            vcf = VcfAnnotation(args.T_prefix + base, args.threads)
            vcf.tumor_annovar()
            vcf.hgmd_annovar()
            vcf.snpEff_Anno()
            if 'brca' in args.panel:  # 是否用brca数据库注释
                vcf.brca_annovar()
                BrcaResult(args.T_prefix + base)
            Hgmd_pipe(args.T_prefix + base)
            Multi2Raw(args.T_prefix + base, args.S_num,args.S_type).multi2raw()
            MultiDrug(args.T_prefix + base, args.cancer, args.S_num).multi_drug_pipe()
            Chemo_drug_pipe(args.T_prefix + base, bed, args.cancer)
            Omim(args.T_prefix + base)
        else: #normal样本
            QC(args.N_prefix, args.threads)
            basicAnalyze(args.N_prefix, args.threads, bed)
            pool1 = Pool(3)
            pool1.apply_async(SingleNormalVarient, args=(args.N_prefix, bed,args.panel, ))
            if 'exon' in args.panel:
                pool1.apply_async(TumorOnlyCNV, args=(args.N_prefix, bed,))
            if args.fusion: #是否分析fusion
                pool1.apply_async(Fusion, args=(args.N_prefix,bed,))
            pool1.close()
            pool1.join()
            base = '.germline'
            vcf = VcfAnnotation(args.N_prefix+base,args.threads)
            vcf.tumor_annovar()
            vcf.hgmd_annovar()
            vcf.snpEff_Anno()
            if 'brca' in args.panel: #是否用brca数据库注释
                vcf.brca_annovar()
                BrcaResult(args.N_prefix + base)
            Hgmd_pipe(args.N_prefix + base)
            Multi2Raw(args.N_prefix + base,args.S_num,args.S_type).multi2raw()
            Germline_tumor_pipe(args.N_prefix + base,args.panel)
            Chemo_drug_pipe(args.N_prefix +base, bed,args.cancer)
            Omim(args.N_prefix +base)
    else: # 配对样本分析
        p1 = Process(target=QC,args=(args.T_prefix, args.threads,))
        p2 = Process(target=QC, args=(args.N_prefix, args.threads,))
        p1.start()
        p2.start()
        p1.join()
        p2.join()
        p3 = Process(target=basicAnalyze,args=(args.T_prefix, args.threads, bed,))
        p4 = Process(target=basicAnalyze, args=(args.N_prefix, args.threads, bed,))
        p5 = Process(target=HLA_type,args=(args.T_prefix, args.threads,))
        p6 = Process(target=HLA_type, args=(args.N_prefix, args.threads,))
        p3.start()
        p4.start()
        p5.start()
        p6.start()
        p3.join()
        p4.join()
        pool2 = Pool(7)
        pool2.apply_async(PonOfNormal,args=(args.T_prefix, args.N_prefix, bed,args.panel))
        pool2.apply_async(SingleNormalVarient,args=(args.N_prefix, bed,args.panel,))
        if args.cnv:
            pool2.apply_async(CNV,args=(args.T_prefix, args.N_prefix, args.threads,bed,))
            pool2.apply_async(Cnvkit,args=(args.T_prefix, args.N_prefix,bed,))
        if args.msi:
            pool2.apply_async(MSI,args=(args.T_prefix, args.N_prefix, bed,))
            pool2.apply_async(MSI5Points,args=(args.T_prefix, args.N_prefix, args.threads,))
        if args.fusion:
            pool2.apply_async(Fusion,args=(args.T_prefix, bed,))
        pool2.close()
        pool2.join()

        Mutect2(args.T_prefix, args.N_prefix, bed,args.threads)
        VcfAnnotation(args.T_prefix+'_vs_'+args.N_prefix+'.somatic',args.threads).tumor_annovar()
        VcfAnnotation(args.T_prefix + '_vs_' + args.N_prefix + '.somatic', args.threads).snpEff_Anno()
        VcfAnnotation(args.N_prefix+'.germline', args.threads).tumor_annovar()
        VcfAnnotation(args.N_prefix+'.germline',  args.threads).hgmd_annovar()
        if 'brca' in args.panel:  # 是否用brca数据库注释
            VcfAnnotation(args.T_prefix + '_vs_' + args.N_prefix + '.somatic', args.threads).brca_annovar()
            BrcaResult(args.T_prefix + '_vs_' + args.N_prefix + '.somatic')
        Hgmd_pipe(args.N_prefix +'.germline')
        Multi2Raw(args.N_prefix +'.germline', 'Single','normal').multi2raw()
        Multi2Raw(args.T_prefix + '_vs_' + args.N_prefix+ '.somatic', args.S_num,args.S_type).multi2raw()
        MultiDrug(args.T_prefix + '_vs_' + args.N_prefix+ '.somatic', args.cancer, args.S_num).multi_drug_pipe()
        Germline_tumor_pipe(args.N_prefix +'.germline',args.panel)
        MMR(args.T_prefix,args.N_prefix)
        Chemo_drug_pipe(args.N_prefix + '.germline',  bed, args.cancer)
        Omim(args.T_prefix+'_vs_'+args.N_prefix+'.somatic')
        Omim(args.N_prefix + '.germline')
        p5.join()
        NeoPred(args.T_prefix, args.N_prefix)
    path_script = os.path.dirname(__file__)
    os.system("ls Result/*xls |while read i;do python "+path_script+"/txt2csv.py $i;done")
    print('分析结束')
    t1 = time.time()
    print('共运行时长：{}'.format(TransTime(t1-t0)))
    print('\n' + '#' * 60)
#
filePath= args.T_prefix + '-' + time.strftime("%Y%m%d")
mkpath= '/home/longzhao/panelSel/cancer_data/' + filePath
data= mkdir(mkpath)
currentPath=os.getcwd()
print(currentPath + '/Result-CSV')
current_folder = os.listdir(currentPath + '/Result-CSV')#current_foder是‘模拟’文件夹下所有子文件名组成的一个列表

# # 第二部分，将名称为file的文件复制到名为file_dir的文件夹中
for x in current_folder:
    #拼接出要存放的文件夹的路径
    file_dir = mkpath + '/' + x
    #将指定的文件file复制到file_dir的文件夹里面
    shutil.copy(currentPath + '/Result-CSV'+'/'+x,file_dir)

url='http://192.168.1.195:14001/library/report/save?sampleId=' + args.T_prefix + '&filePath=' + filePath
r=requests.get(url)
print("status code:",r.status_code)
response_dict=r.json()
print("Total repositories:",response_dict)
print('\n' + '#' * 60)
