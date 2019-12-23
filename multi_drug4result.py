#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/25 
@Description :见类中注释
'''

import time
import os,sys
import re
# import csv
from collections import defaultdict
import functools
from profile import Profile

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

class Multi2Raw():
    '''
    --!!!single/paired/UMI/noUMI不同 vaf和DP信息所在位置不同!!!
    SJR1-1.hg19_multianno.txt
    gatk处理，同一位置，不同突变类型,AD位置有2个变异，分别对应第一二行
    chr  start         end       ref alt         maf             vcf-alt            GT:AD:DP:GQ:PL
    6	31324666	31324666	A	T  ...     0.1521      ...   T,C   ...    1/2:0,166,209:375:99:16013,8726,8231,7002,0,6351
    6	31324666	31324666	A	C ...      0.6686      ...   T,C   ...    1/2:0,166,209:375:99:16013,8726,8231,7002,0,6351
    过滤条件：
    1.只关注外显子区，过滤掉无义突变
    2.肿瘤样本：过滤掉突变频率低于0.01，深度低于100，alt AD小于2的
    3.过滤掉人群中突变频率大于0.05的
    4.配对分析独有条件：过滤掉正常样本变异频率大于0.03的
    5.出现多个转录本信息时，若该基因收录在geneNM4resultProcess.txt*文件中（来自OncoKB的allAnnotatedVariants.txt文件），则选择相应基因的参考转录本，否则默认保留第一个。
    6. 同一位置出现2个基因时，选择bed文件中的那个，若均不在bed文件默认选第一个。
    7. 将1个字母的AA转换为3个字母的AA。
    '''
    def __init__(self,prefix,num,type):
        AA_dict = {'G': 'Gly', 'A': 'Ala', 'V': 'Val', 'L': 'Leu', 'I': 'Ile', 'P': 'Pro', 'F': 'Phe',
                   'Y': 'Tyr', 'W': 'Trp', 'S': 'Ser', 'T': 'Thr', 'C': 'Cys', 'M': 'Met', 'N': 'Asn',
                   'Q': 'Gln', 'D': 'Asp', 'E': 'Glu', 'K': 'Lys', 'R': 'Arg', 'H': 'His'}
        self.prefix = prefix
        self.AAdict = AA_dict
        self.num = num
        self.type = type
        self.var_parh = Profile()
        self.path =  'Result'

    #1字母AA转换为3字母
    def transAA(self,pro):
        pro2 = re.findall('[A-Z]+', pro)
        for p in pro2:
            if p in self.AAdict:
                pro = re.sub(p, self.AAdict[p], pro)
        return pro

    # 参考转录本文件生成字典
    def NM4gene(self):
        dic = {}
        with open(self.var_parh['NM_file']) as f:
            for i in f:
                ii = i.strip().split('\t')
                dic[ii[1]] = ii[0]
        return dic

    # 选择转录本，获取相应转录本信息
    def processAA(self,str): # CSF3R:NM_000760:exon10:c.T1260C:p.T420T,CSF3R:NM_156039:exon10:c.T1260C:p.T420T
        if str == '.' or str == 'UNKNOWN':
            trans, exon,cdna, pro,pro_drug = '-', '-', '-','-','-'
        else:
            trans_list = str.split(',')
            NM4gene_dic = self.NM4gene()
            x = [tran for tran in trans_list for gene in NM4gene_dic if NM4gene_dic[gene] in tran]
            if len(x) == 0:  # 如果不在参考转录本基因列表，默认选第一个转录本
                transs = trans_list[0].split(':')
            else:  # 否则，选择该基因的参考转录本
                transs = x[0].split(':')
            try:
                trans,exon,cdna = transs[1:4]
            except:
                trans,exon, cdna = '-','-','-'
            if len(transs) == 5:
                pro_raw = transs[4] #结果中原始蛋白
                pro = self.transAA(pro_raw) #转换为3字母蛋白
                pro_drug = pro_raw.replace('p.','')
            else:
                pro_drug,pro = '-','-'
        return trans,exon, cdna, pro,pro_drug

    # 同一位置出现2个基因名称时，选择genePosition文件中的基因名称
    def geneName(self,gene):
        num_gene = gene.split(';')
        if len(num_gene) == 2:
            gene_lis = list()
            with open(self.var_parh['genePosition']) as gene_pos:
                for line in gene_pos:
                    lin = line.strip().split('\t')[3]
                    gene_lis.append(lin)
            gene_lis2 = list(set(gene_lis))
            gene2 = [i if i in gene_lis2 else num_gene[0] for i in num_gene ][0]
        else:
            gene2 = gene
        return gene2

    #计算单样本分析的vaf/dp等值
    def vafSingleNormal(self,line,pos_dict):  # 顺序是 GT:AD:DP:GQ:PL
        col4 = '\t'.join(line[:4])  # 前4列：1	2488153	2488153	A
        pos_dict[col4] += 1  # 出现次数
        tumor_line = line[-1].split(':')  # tumor在倒数第1列
        genetype = tumor_line[0]
        if genetype == '0/0' or genetype == '1/1':
            gt = '纯和'
        else:
            gt = '杂合'
        if pos_dict[col4] == 1:  # 1/2:0,166
            tumor_alt_num = tumor_line[1].split(',')[1]
        else:  # 出现次数为2时，1/2:0,209
            tumor_alt_num = tumor_line[1].split(',')[-1]
        tumor_dp = tumor_line[2]
        try:
            tumor_vaf = int(tumor_alt_num) / float(tumor_dp)
        except:
            tumor_vaf = 0
        return tumor_vaf, tumor_dp,tumor_alt_num,gt

    def vafSingleTumor(self,line,pos_dict):  # 顺序是 GT:AD:AF:DP...
        col4 = '\t'.join(line[:4])  # 前4列：1	2488153	2488153	A
        pos_dict[col4] += 1  # 出现次数
        tumor_line = line[-1].split(':')  # tumor在倒数第1列
        if pos_dict[col4] == 1:  # 1/2:0,166
            tumor_alt = tumor_line[1].split(',')[1]
        else:  # 出现次数为2时，1/2:0,209
            tumor_alt = tumor_line[1].split(',')[-1]
        tumor_dp = tumor_line[3]
        tumor_vaf = int(tumor_alt) / int(tumor_dp)
        return tumor_vaf, tumor_dp, int(tumor_alt)


    # 计算配对样本分析的vaf/dp等值
    def vafPaired(self, line,pos_dict): #顺序是GT:AD:AF:DP:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB
        last_line = line[-1].split(':')
        last_2_line = line[-2].split(':')
        if last_line[0] in ['0/1','0|1'] and last_2_line[0] in ['0/0','0|0']:
            tumor_line = last_line
            normal_line = last_2_line
            tumor_ref, tumor_alt = tumor_line[1].split(',')
            normal_ref, normal_alt = normal_line[1].split(',')
            tumor_dp = int(tumor_ref) + int(tumor_alt)
            normal_dp = int(normal_ref) + int(normal_alt)
            tumor_vaf = int(tumor_alt) / tumor_dp
            normal_vaf = int(normal_alt) / normal_dp
        elif last_line[0] == '0/1/2' and last_2_line[0] == '0/0':
            pos = '\t'.join(line[-11:-9])  # 倒数10/11列：1	2488153
            pos_dict[pos] += 1  # 出现次数
            tumor_line = last_line  # tumor在倒数第2列
            normal_line = last_2_line
            if pos_dict[pos] == 1:  # 1/2:0,166
                tumor_alt = tumor_line[1].split(',')[1]
                normal_alt = normal_line[1].split(',')[1]
            else:  # 出现次数为2时，1/2:0,209
                tumor_alt = tumor_line[1].split(',')[-1]
                normal_alt = normal_line[1].split(',')[-1]
            tumor_dp = sum([int(i) for i in tumor_line[1].split(',')])
            normal_dp = sum([int(i) for i in normal_line[1].split(',')])
            tumor_vaf = int(tumor_alt) / tumor_dp
            normal_vaf = int(normal_alt) / normal_dp
        elif last_line[0] in ['0/0','0|0'] and last_2_line[0] in ['0/1','0|1']:
            tumor_line = last_2_line
            normal_line = last_line
            tumor_ref, tumor_alt = tumor_line[1].split(',')
            normal_ref, normal_alt = normal_line[1].split(',')
            tumor_dp = int(tumor_ref) + int(tumor_alt)
            normal_dp = int(normal_ref) + int(normal_alt)
            tumor_vaf = int(tumor_alt) / tumor_dp
            normal_vaf = int(normal_alt) / normal_dp
        elif last_line[0] == '0/0' and last_2_line[0] == '0/1/2':
            pos = '\t'.join(line[-11:-9])  # 倒数10/11列：1	2488153
            pos_dict[pos] += 1  # 出现次数
            tumor_line = last_2_line  # tumor在倒数第2列
            normal_line = last_line
            if pos_dict[pos] == 1:  # 1/2:0,166
                tumor_alt = tumor_line[1].split(',')[1]
                normal_alt = normal_line[1].split(',')[1]
            else:  # 出现次数为2时，1/2:0,209
                tumor_alt = tumor_line[1].split(',')[-1]
                normal_alt = normal_line[1].split(',')[-1]
            tumor_dp = sum([int(i) for i in tumor_line[1].split(',')])
            normal_dp = sum([int(i) for i in normal_line[1].split(',')])
            tumor_vaf = int(tumor_alt) / tumor_dp
            normal_vaf = int(normal_alt) / normal_dp
        else:
            print('数据异常')
            tumor_dp, normal_dp, tumor_vaf, normal_vaf, normal_alt,tumor_alt = 0, 0, 0, 0, 0,0
            err = open('Result/' + self.prefix + '_err.data', 'a')
            err.write('\t'.join(line)+'\n')
            err.close()
        return tumor_vaf, normal_vaf, tumor_dp, normal_dp, normal_alt,int(tumor_alt)

    #multi生成raw_result文件
    @timefly
    def multi2raw(self):
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        os.system("cp VarientCalling/{t}.hg19_multianno.txt {p}".format(t=self.prefix,p=self.path))
        multi_file = '{p}/{t}.hg19_multianno.txt'.format(t=self.prefix,p=self.path)
        out_raw = open('{p}/{t}.raw_result.xls'.format(t=self.prefix,p=self.path), 'w')
        out_raw.write('Chr\tStart\tEnd\tRef\tAlt\tGene.refGene\tExonicFunc.refGene\tTrans_cDNA_AA\tcosmic8N\tavsnp150\tCLNDN\tCLNSIG\tInterVar_automated\tExAC_EAS\tT_VAF\tT_DP\tN_VAF\tN_DP\n')
        row = 0
        pos_dict = defaultdict(int)
        with open(multi_file) as f:
            for line in f:
                if not line.startswith('Chr'):
                    lin = line.strip().split('\t')
                    gene_name = self.geneName(lin[6])
                    trans, exon, cdna, pro, pro_drug = self.processAA(lin[9])  # TNFRSF14:NM_001297605:exon1:c.A50G:p.K17R,TNFRSF14:NM_003820:exon1:c.A50G:p.K17R
                    Trans_cDNA_AA = '|'.join([trans, exon, cdna, pro, pro_drug])
                    if lin[5] == 'exonic' and lin[8] != 'synonymous SNV':  # 只关注外显子区，过滤掉无义突变
                        if lin[15] == '.' or float(lin[15]) <= 0.05:  # 过滤掉人群中突变频率大于0.05的
                            if self.num == 'Single':
                                if self.type == 'tumor' : # single tumor 顺序是 GT:AD:DP:GQ:PL
                                    tumor_vaf, tumor_dp,tumor_alt = self.vafSingleTumor(lin,pos_dict)
                                    tumor_vaf_str = "%.2f" % (100 * tumor_vaf)
                                    if  tumor_vaf >= 0.01 and int(tumor_dp) >= 100 and int(tumor_alt) >=2:
                                        row += 1
                                        info_raw = '\t'.join(lin[:5]) + '\t' +gene_name+ '\t' +lin[8]+'\t'+Trans_cDNA_AA+'\t'+'\t'.join(lin[10:12]) + '\t' +lin[22]+'\t' + '\t'.join(lin[25:27]) + '\t' + lin[15] + '\t' + tumor_vaf_str + '%\t' + tumor_dp + '\n'
                                        out_raw.write(info_raw)
                                else: #single normal 顺序是 GT:AD:DP:GQ:PL
                                    tumor_vaf, tumor_dp, tumor_alt, gt = self.vafSingleNormal(lin,pos_dict)
                                    tumor_vaf_str = "%.2f" % (100 * tumor_vaf)
                                    if tumor_vaf >= 0.01 and int(tumor_dp) >= 50 and int(tumor_alt) >= 2:
                                        row += 1
                                        info_raw = '\t'.join(lin[:5]) + '\t' + gene_name + '\t' + lin[8] + '\t' + Trans_cDNA_AA + '\t' + '\t'.join(lin[10:12]) + '\t' + lin[22] + '\t' + '\t'.join(lin[25:27]) + '\t' + lin[15] + '\t' + tumor_vaf_str + '%\t' + tumor_dp + '\t' + gt + '\n'
                                        out_raw.write(info_raw)
                            else: # paired  顺序是GT:AD:AF:DP:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB
                                tumor_vaf, normal_vaf, tumor_dp, normal_dp, normal_alt,tumor_alt = self.vafPaired(lin,pos_dict)
                                tumor_vaf_str = "%.2f" % (100 * tumor_vaf)
                                normal_vaf_str = "%.2f" % (100 * normal_vaf)
                                if  tumor_vaf >= 0.01 and tumor_dp >= 100 and tumor_alt >= 2:
                                    if normal_vaf <= 0.03:
                                        row += 1
                                        info_raw = '\t'.join(lin[:5]) + '\t' +gene_name+ '\t' +lin[8]+'\t'+Trans_cDNA_AA+'\t'+'\t'.join(lin[10:12])+ '\t' +lin[22]+'\t'  + '\t'.join(lin[25:27]) + '\t' + lin[15] + '\t' + tumor_vaf_str + '%\t' + str(tumor_dp) + '\t' + normal_vaf_str + '%\t' + str(normal_dp) + '\n'
                                        out_raw.write(info_raw)

        bed_length = open('StateMap/stat_bam.xls').readlines()[-1].strip().split('\t')[-2]
        bed_length = int(bed_length) / 1000000
        tmb = row / bed_length
        out_raw.write('TMB:\t' + str(tmb))
        out_raw.close()


class MultiDrug():
    '''
    1.如果“gene_drug_relation".txt表格中存在基因突变-药物关系（包括临床药物），就放入1/2类变异，没有推荐药物的基因突变放入3类。
    2.不同癌症的对应基因信息.XLS：这个表格中每个癌症包括24基因（包括I类\II类）。如果I类或II类基因，没有突变，需要以‘-’的形式写入，见示例报告。
    3.举个例子，如果BRCA突变，在卵巢癌中就是1类变异，在肺癌中就是2类变异。（我在卵巢癌中会标注成1类变异，肺癌中不会进行标注）
    4.Truncating Mutations=stopgain，无义突变
    5.Oncogenic Mutations：暂定为somatic的multianno.txt的CLNSIG/InterVar_automated注释为（pathogenic、likely pathogenic 和uncertain（.）；致病的、可能致病的、意义未明的。），对应的就是不包含（benign，likely benignb，良性、可能良性的）。因为需要看CLNSIG/InterVar_automated这两列：CLNSIG为CLINVAR的注释结果，InterVar_automated为ACMG的注释结果，因此只要其中一列为致病的、可能致病的，就可以进行保留）
    6.CNV和fusion的筛选规则：只要“gene_drug_relation.txt"表格内容的基因，其他基因暂时不需要提取。简而言之：在gene_drug_relation.txt该表格中就是有临床意义，不在该表格，就暂无临床意义。
    7.针对某癌症如肺癌，如果是肺癌/实体瘤/所有肿瘤类，变异类别/等级和表格中一致；如果是其他癌症，则类别/等级统一为 II类/C级
    '''
    def __init__(self,prefix,cancer,num):
        self.prefix = prefix
        self.cancer = cancer
        self.num = num
        self.var_parh = Profile()
        self.path =  'Result'

    #SNV/INDEL 基因-药物关系 字典
    def get_gene_drug(self):
        gene_drug_dict = {}
        with open(self.var_parh['geneDrug_file'], encoding='utf-8')  as genedrug:
            for line in genedrug:
                if line != '':
                    lin = line.strip().split('\t')
                    NM = lin[3].split('.')[0]
                    gene = lin[2]
                    AA = lin[4]  # 注意格式
                    drug = lin[6]
                    gene_info = '&'.join([gene, NM, AA])
                    if line.startswith(self.cancer) or line.startswith('Solid') or line.startswith('All') :
                        type = lin[5]
                        response_level = lin[7] +'/'+ lin[8]
                        gene_drug_dict[gene_info] = [type, drug, response_level]
                    else: #其他肿瘤均为II类/C级
                        type = 'II'
                        response_level = lin[7] + '/' + 'C级'
                        if not gene_info in gene_drug_dict:
                            gene_drug_dict[gene_info] = [type, drug, response_level]
        return gene_drug_dict

    #变异类型英文转中文
    def trans_varient_type(self,str):
        dic = {}
        with open(self.var_parh['varientType'],encoding='utf-8') as varientType:
            for line in varientType:
                lin = line.strip().split('\t')
                dic[lin[0]] = lin[1]
        if str in dic:
            str = dic[str]
        return str

    # multi SNV/INDEL结果添加到drug表格
    def raw2drug(self):
        rawresult_file = 'Result/{}.raw_result.xls'.format(self.prefix)
        if not os.path.exists(rawresult_file):
            print('raw_result文件不存在！')
            sys.exit(1)
        # out_simple = open('{p}/{t}.result.xls'.format(t=self.prefix, p=self.path), 'w', encoding='gb18030')
        # out_simple = csv.writer(out_simple)
        # out_simple.writerow('所属类别,基因,cDNA改变,氨基酸改变,变异类型,变异率,靶向药物,敏感性/证据等级'.split(','))
        out_simple = open('{p}/{t}.result.xls'.format(t=self.prefix, p=self.path), 'w', encoding='utf-8')
        with open(rawresult_file) as rawresult:
            for raw in rawresult:
                if not (raw.startswith('TMB') or raw.startswith('Chr')) :
                    ra = raw.strip().split('\t')
                    trans, exon, cdna, pro, pro_drug = ra[7].split('|')  # TNFRSF14:NM_001297605:exon1:c.A50G:p.K17R,TNFRSF14:NM_003820:exon1:c.A50G:p.K17R
                    gene_name = ra[5]
                    varient_type = self.trans_varient_type(ra[6])
                    tumor_vaf_str = ra[14]
                    gene_info = '&'.join([gene_name, trans, pro_drug])
                    gene_info2 = '&'.join([gene_name, trans, 'TruncatingMutations'])
                    gene_info3 = '&'.join([gene_name, trans, 'OncogenicMutations'])
                    gene_drug_dict = self.get_gene_drug()
                    if gene_info in gene_drug_dict:
                        type, drug, response_level = gene_drug_dict[gene_info]
                    elif gene_info2 in gene_drug_dict and ra[6] == 'stopgain':
                        type, drug, response_level = gene_drug_dict[gene_info2]
                    elif gene_info3 in gene_drug_dict and (ra[11] and ra[12]) in ['pathogenic', 'likely pathogenic','Uncertain_significance','Uncertain significance', 'Conflicting_interpretations_of_pathogenicity','not_provided','.']:
                        type, drug, response_level = gene_drug_dict[gene_info3]
                    else:
                        type, drug, response_level = 'III', '-', '-'
                    simple = [type, gene_name, cdna, pro,trans+':'+exon, varient_type, tumor_vaf_str, drug , response_level]
                    # out_simple.writerow(simple)
                    out_simple.write('\t'.join(simple)+'\n')

    #提取fusion结果添加到drug表格
    def fusion2drug(self):
        # out_simple = open('{p}/{t}.result.xls'.format(t=self.prefix, p=self.path), 'a', encoding='gb18030')
        # out_simple = csv.writer(out_simple)
        out_simple = open('{p}/{t}.result.xls'.format(t=self.prefix, p=self.path), 'a', encoding='utf-8')
        fusion_file ='Result/{t}.fusions_result.xls'.format(t=self.prefix.split('.')[0].split('_')[0])
        if not os.path.exists(fusion_file):
            print('缺少Fusion结果文件！')
            sys.exit(1)
        gene_drug_dict = {}
        with open(fusion_file, encoding='utf-8') as fusion:
            for line in fusion:
                if line != '' and not line.startswith('Est_Type'):
                    lin = line.strip().split('\t')
                    gene_fusion = lin[1:3] #两融合基因list
                    if (lin[11] and lin[16]) != '':
                        vaf = float(lin[11])/float(lin[16])
                        vaf_str = "%.2f" % (100 * vaf)+'%'
                    else:
                        vaf_str = '-'
                    with open(self.var_parh['geneDrug_file'], encoding='utf-8')  as genedrug:
                        for genedrg in genedrug:
                            genedr = genedrg.strip().split('\t')
                            gene_id = genedr[2]
                            if  gene_id in gene_fusion and 'Fusion' in genedr[4] :
                                cdna = gene_id + '-' + gene_fusion[1 - gene_fusion.index(gene_id)]
                                dic_key = '&'.join([gene_id, cdna, vaf_str, genedr[6], genedr[7]])  # gene-cdna-vaf-drug-response
                                if genedrg.startswith(self.cancer) or genedrg.startswith('实体瘤') or genedrg.startswith('所有肿瘤'):
                                    type = genedr[5]
                                    level = genedr[8]
                                    gene_drug_dict[dic_key] = [type,level]
                                else:
                                    type = 'II'
                                    level = 'C级'
                                    if not gene_id in gene_drug_dict: #如果该基因并不在I类变异
                                        gene_drug_dict[dic_key] = [type,level]
        for k, v in gene_drug_dict.items():
            type, level = v
            gene, cdna,vaf_str, drug, response = k.split('&')
            info = [type, gene, cdna, '-','-', '融合', vaf_str, drug, response + '/' + level]
            # out_simple.writerow(info)
            out_simple.write('\t'.join(info) + '\n')

    # 取cnvkit和freec交集,未使用，代码需更新
    def process_cnv_jiaoji(self):
        cnvkit_file = 'CNVCalling/{t}.cnvkit_result.xls'.format(t=self.prefix)
        freec_file = 'CNVCalling/{t}.freec_result.xls'.format(t=self.prefix)
        out_simple = open('{p}/{t}.result.xls'.format(t=self.prefix, p=self.path), 'a', encoding='gb18030')
        out_simple = csv.writer(out_simple)
        if os.path.exists(cnvkit_file) and os.path.exists(freec_file):
            common_file = '{p}/cnv_common.txt'.format(p=self.path)
            common_out = open(common_file,'w',encoding='utf-8')
            common = os.popen("{bedtools} intersect -a {c} -b {f} -wb".format(c=cnvkit_file,f=freec_file,**self.var_parh)).read()
            for line in common.split('\n'):
                if line != '':
                    lin = line.split('\t')
                    if lin[7] == lin[14] and int(lin[7]) >3: #CNV扩增,4个及以上
                        info = '\t'.join(lin[:3])+'\t'+lin[7]+'\n'
                        common_out.write(info)
            common_out.close()
            if os.path.getsize(common_file) : #判断两cnv软件共同结果文件是否为空
                #与genepos文件比较，提取基因名称
                cnv_result = os.popen("{bedtools} intersect -a {co} -b {genePosition} -wa -wb | {bedtools} groupby -i - -g 1,2,3,4 -c 8 -o collapse".format(co=common_file,**self.var_parh)).read()
                for cnvre in cnv_result.split('\n'):
                    if cnvre != '': # 3 gene1,gene2,gene3,,,
                        cnv = cnvre.split('\t')
                        cp = cnv[3]
                        gene_list = cnv[4].split(',')
                        with open(self.var_parh['geneDrug_file'], encoding='utf-8')  as genedrug:
                            for genedrg in genedrug:
                                genedr = genedrg.strip().split('\t')
                                pass #此处需再编写代码
            else:
                print('两cnv分析软件交集为空！')
        else:
            print("cnv 结果不全！")

    # 提取cnv结果添加到drug表格
    def cnv2drug(self):
        cnvkit_file = 'Result/{t}.CNV_cnvkit_result.xls'.format(t=self.prefix)
        # out_simple = open('{p}/{t}.result.xls'.format(t=self.prefix, p=self.path), 'a', encoding='gb18030')
        # out_simple = csv.writer(out_simple)
        out_simple = open('{p}/{t}.result.xls'.format(t=self.prefix, p=self.path), 'a', encoding='utf-8')
        if not os.path.exists(cnvkit_file):
            print('缺少cnvkit结果文件！')
            sys.exit(1)
        gene_drug_dict = {}
        with open(cnvkit_file) as cnvkit:
                for cnvk in cnvkit:
                    cnv = cnvk.strip().split('\t')
                    cp = cnv[7]
                    gene_list = cnv[3].split(',')
                    with open(self.var_parh['geneDrug_file'], encoding='utf-8')  as genedrug:
                        for genedrg in genedrug:
                            genedr = genedrg.strip().split('\t')
                            gene_id = genedr[2]
                            dic_key ='&'.join([gene_id,cp,genedr[6],genedr[7]])#gene-cdna-drug-response
                            if gene_id in gene_list and genedr[4] == 'Amplification' and int(cp) >3:
                                if genedrg.startswith(self.cancer) or genedrg.startswith('实体瘤') or genedrg.startswith('所有肿瘤'):
                                    type = genedr[5]
                                    level =  genedr[8]
                                    gene_drug_dict[dic_key] = [type, level]
                                else:
                                    type = 'II'
                                    level = 'C级'
                                    if not dic_key in gene_drug_dict:  # 如果该基因并不在I类变异
                                        gene_drug_dict[dic_key] = [type, level]
        for k,v in gene_drug_dict.items():
            type, level = v
            gene,cdna,drug,response = k.split('&')
            cdna2 = 'CN:' +cdna
            info = [type,gene,cdna2,'-','-','扩增','-',drug,response+'/'+level]
            # out_simple.writerow(info)
            out_simple.write('\t'.join(info) + '\n')

    #肿瘤变异I/II类空白信息添加到drug表格
    def get_class12_bank(self):#lung_class1 ALK	未检出	-	-	-	-	-	I
        dic=defaultdict(list)
        with open(self.var_parh['classGeneTumor'], encoding='utf-8')  as tumor_gene:
           for line in tumor_gene:
                lin = line.strip().split('\t')
                dic[lin[0]].append(lin[1:])
        # [['ALK', '未检出', '-', '-', '-', '-', '-', 'I'], ['BRAF', '未检出', '-', '-', '-', '-', '-', 'I'],,,]
        dic_key1 = self.cancer+'_class1'
        dic_key2 = self.cancer + '_class2'
        class1_cancer = dic[dic_key1]
        class2_cancer = dic[dic_key2]
        result_file = '{p}/{t}.result.xls'.format(t=self.prefix, p=self.path)
        with open(result_file, encoding='utf-8') as simple:
            # row = csv.reader(simple)
            for row in simple:
                if row != '':
                    i = row.strip().split('\t')
            # for i in row:
            #     if i != []:
                    gene_name = i[1]
                    if i[0] == 'I':
                        # gene出现在['ALK', '未检出', '-', '-', '-', '-', '-', 'I']中，移除该数据，否则打印
                        for class1 in class1_cancer:
                            if gene_name in class1:
                                class1_cancer.remove(class1)
                    if i[0] == 'II':
                        for class2 in class2_cancer:
                            if gene_name in class2:
                                class2_cancer.remove(class2)
        # out_simple = open('{p}/{t}.result.xls'.format(t=self.prefix, p=self.path), 'a', encoding='gb18030')
        # out_simple = csv.writer(out_simple)
        out_simple = open('{p}/{t}.result.xls'.format(t=self.prefix, p=self.path), 'a', encoding='utf-8')
        for cl1 in class1_cancer:
            # out_simple.writerow(cl1)
            out_simple.write('\t'.join(cl1)+'\n')
        for cl2 in class2_cancer:
            # out_simple.writerow(cl2)
            out_simple.write('\t'.join(cl2) + '\n')
        out_simple.close()

    #基因-功能对应结果
    def gene_function_result(self,genelist,name):
        genelist = list(set(genelist))
        geneExp = open('Result/{0}.{1}_geneInfo_result.xls'.format(self.prefix,name),'w',encoding='utf-8')
        with open(self.var_parh['geneInfo_file']) as info:
            for line in info:
                lin = line.strip().split('\t')
                for i in genelist:
                    if lin[0] == i:
                        geneExp.write(lin[0]+'\t'+lin[1]+'\n')
        geneExp.close()

    #药物及说明结果文件
    def drug_information_result(self,druglist,name):
        druglist = list(set(druglist))
        drugExp = open('Result/{0}.{1}_drugInfo_result.xls'.format(self.prefix,name),'w',encoding='utf-8')
        with open(self.var_parh['drugInfo_file']) as info :
            for line in info:
                lin = line.strip().split('\t')
                for i in druglist:
                    if lin[0] == i:
                        drugExp.write(lin[0]+'\t'+lin[-1]+'\n')
                    if lin[1] == i:
                        drugExp.write(lin[1] + '\t' + lin[-1] + '\n')
        drugExp.close()

    def split_file(self):
        result_file = '{p}/{t}.result.xls'.format(t=self.prefix, p=self.path)
        class_I_file =  open('{p}/{t}.result.class_I.xls'.format(t=self.prefix, p=self.path), 'w', encoding='utf-8')
        class_II_file = open('{p}/{t}.result.class_II.xls'.format(t=self.prefix, p=self.path), 'w', encoding='utf-8')
        class_III_file = open('{p}/{t}.result.class_III.xls'.format(t=self.prefix, p=self.path), 'w', encoding='utf-8')
        genelist_I,genelist_II,genelist_III = [],[],[]
        genelist_drug_I, genelist_drug_II= [], []
        # with open(result_file, encoding='gb18030') as simple:
        #     row = csv.reader(simple)
        #     for i in row:
        with open(result_file, encoding='utf-8') as simple:
            for row in simple:
                i = row.strip().split('\t')
                if i[0] == 'I' :
                    class_I_file.write('\t'.join(i[1:])+'\n')
                    if  i[2] != '未检出':
                        genelist_I.append(i[1])
                        drug = i[6].replace('+', ',').replace('|', ',').split(',')  #将药物按‘+’和‘|’分割
                        genelist_drug_I.extend(drug)
                if i[0] == 'II' :
                    class_II_file.write('\t'.join(i[1:])+'\n')
                    if i[2] != '未检出':
                        genelist_II.append(i[1])
                        drug = i[6].replace('+', ',').replace('|', ',').split(',')  # 将药物按‘+’和‘|’分割
                        genelist_drug_II.extend(drug)
                if i[0] == 'III' :
                    class_III_file.write('\t'.join(i[1:])+'\n')
                    if i[2] != '未检出':
                        genelist_III.append(i[1])
        class_I_file.close()
        class_II_file.close()
        class_III_file.close()
        self.gene_function_result(genelist_I,'class_I')
        self.gene_function_result(genelist_II, 'class_II')
        self.gene_function_result(genelist_III, 'class_III')
        self.drug_information_result(genelist_drug_I,'class_I')
        self.drug_information_result(genelist_drug_II, 'class_II')

    @timefly
    def multi_drug_pipe(self):
        self.raw2drug()
        self.fusion2drug()
        if self.num == 'Paired':
            self.cnv2drug()
        self.get_class12_bank()
        self.split_file()


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('\nUsage:  python {} [prefix] [Single/Paired] [tumor/normal] [NSCLC,Breast,Colorectal,Gastric_Cancer,Esophagus_Cancer,Gastrointestinal_Stromal]\n'.format(sys.argv[0]))
        sys.exit(1)
    print(time.ctime(), 'MultiDrug begin')
    prefix = sys.argv[1]
    sampleNum = sys.argv[2]
    sampleType = sys.argv[3]
    cancer = sys.argv[4]
    Multi2Raw(prefix,sampleNum,sampleType).multi2raw()
    MultiDrug(prefix,cancer,sampleNum).multi_drug_pipe()
