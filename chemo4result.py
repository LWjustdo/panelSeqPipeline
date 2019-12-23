#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/3 
@Description :
针对germline样本,化疗位点注释；输入化疗位点列表文件chemoRS.txt/bed文件/uncover.bed文件/变异注释文件hg19_multianno.txt;
再结合化疗药物文件chemoDrug4result.txt，提取相应信息；按高等级（1A/1B/2A/2B），判断增加/降低/-
'''
import time
import sys,os
import functools
from collections import defaultdict
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

#得到化疗位点
def chemo(prefix,bed):
    """
    input：化疗位点文件chemoRS.txt：1	97981343	97981343	rs55886062	A
                                        bed文件：1	97981343	97981343
                                 uncover.bed: 1 97981343	97981343
      注释结果*.hg19_multianno.txt：1	97981343	97981343	C	A
    :return:
    """
    multi_file = 'Result/' + prefix + '.hg19_multianno.txt'
    if not os.path.exists('Result'):
        os.mkdir('Result')
    if not os.path.exists(multi_file):
        os.system("cp VarientCalling/{0}.hg19_multianno.txt Result".format(prefix))
    chem_out = open(multi_file.split('hg19')[0] + 'chemo_raw_result.xls', 'w')
    #取bed文件与uncover.bed文件差集，得到cover.bed文件
    uncover_bed = 'StateMap/{t}/uncover.bed'.format(t=prefix.split('.')[0])
    cover_bed = os.popen("{bedtools} subtract -a {b} -b {unb} ".format(b=bed,unb=uncover_bed,**var_path)).read().strip().split('\n')
    dic_che, dic_mut = {}, {}
    for be in cover_bed:
        b = be.split('\t')
        with open(var_path['chemoRS']) as chemo:
            for chem in chemo:
                che = chem.strip().split('\t')
                if b[0] == che[0] and b[1] <= che[1] and b[2] >= che[2]:
                    pos = '-'.join(che[:3])  # 1-123-456
                    mut = che[3:5]  # rs16947 A
                    dic_che[pos] = mut
    with open(multi_file) as multi:
        for mu in multi:
            m = mu.strip().split('\t')
            pos2 = '-'.join(m[:3])  # 1-123-456
            mut2 = m[3:5]  # A T
            genetype = m[-1].split(':')[0]  # 0/1 或1/1
            mut2.append(genetype)
            dic_mut[pos2] = mut2
    for k, v in dic_che.items():
        if k in dic_mut:
            if dic_mut[k][2] == '1/1':  # 纯合变异
                base = dic_mut[k][1]*2
            else:  # 杂合变异
                base = dic_mut[k][0] + dic_mut[k][1]
            line = prefix.split('.')[0] + '\t' + dic_che[k][0] + '\t' + base + '\n'
        else:  # 无变异
            line = prefix.split('.')[0] + '\t' + dic_che[k][0] + '\t' + dic_che[k][1]*2 + '\n'
        chem_out.write(line)
    chem_out.close()

#指定肿瘤对应的化疗药物信息
def drug4cancer():
    cancer_dict = defaultdict(list)
    with open(var_path['chemoCancerDrug'],encoding='utf-8') as cancer:  # cancer-癌种-drug
        for canc in cancer:
            can = canc.strip().split('\t')
            cancer_dict[can[0]].append(can[-1])
    return cancer_dict

#化疗位点对应药物信息
def chemo2drug(prefix,cancer):
    out_chem_drug = open('Result/{}.chemo_drug_result.xls'.format(prefix),'w',encoding='utf-8')
    out_chem_drug.write('药物\t基因\t位点\t结果\t风险用药提示\t研究人群\t毒副作用\t疗效\t剂量\t等级\n')
    chemo_file =  'Result/{}.chemo_raw_result.xls'.format(prefix)
    with open(var_path['chemoDrug'],encoding='utf-8') as chemo_drug:
        for drug in chemo_drug: #drug-gene-rs-genetype
            dru = drug.strip().split('\t')
            if dru[0] in drug4cancer()[cancer]: #提取指定肿瘤对应药物
                with open(chemo_file) as chemo: #sample-rs-genetype
                    for chem in chemo:
                        che = chem.strip().split('\t')
                        if che[1] == dru[2] and (che[2] == dru[3] or che[2] == dru[3][::-1]):
                            out_chem_drug.write(drug)
        out_chem_drug.close()

#给出化疗位点药物提示
def drug_tips(prefix):
    drug_tip_file = open('Result/{}.chemo_drug_tip_result.xls'.format(prefix),'w',encoding='utf-8')
    chem_drug = 'Result/{}.chemo_drug_result.xls'.format(prefix)
    drug_dic = defaultdict(list)
    dic = {'增加':1,'降低':-1,'-':0}
    with open(chem_drug,encoding='utf-8') as chemdrug:
        for chdr in chemdrug:
            if not chdr.startswith('药物'):
                chd = chdr.strip().split('\t')
                drug_dic[chd[0]].append(chd[-4:]) #药物:[毒副作用,疗效,剂量,等级]
    for k,v in drug_dic.items():
        level_list = [i[-1] for i in v]
        top_level = sorted(level_list)[0]
        toxic = sum([dic[i[-4]] for i in v if top_level in i]) #最高级别对应毒副作用总和
        curative = sum( [dic[i[-3]] for i in v if top_level in i])#最高级别对应疗效总和
        dosage =  sum([dic[i[-2]] for i in v if top_level in i])#最高级别对应剂量总和
        if toxic > 0:
            toxic2 = '毒副作用增加'
        elif toxic < 0:
            toxic2 = '毒副作用降低'
        else:
            toxic2 = ''
        if curative > 0:
            curative2 = '疗效增加'
        elif curative < 0:
            curative2 = '疗效降低'
        else:
            curative2 = ''
        if dosage > 0:
            dosage2 = '剂量增加'
        elif dosage < 0:
            dosage2 = '剂量降低'
        else:
            dosage2 = ''
        drug_tip_file.write(' '.join([k,toxic2,curative2,dosage2])+'\n')
    drug_tip_file.close()

@timefly
def Chemo_drug_pipe(prefix,bed,cancer):
    chemo(prefix, bed)
    chemo2drug(prefix,cancer)
    drug_tips(prefix)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('\nUsage: python {} [prefix] [bed] [NSCLC,Breast,Colorectal,Gastric_cancer,Esophagus_cancer,Gastrointestinal_Stromal]\n'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    bed = sys.argv[2]
    cancer = sys.argv[3]
    print(time.ctime(),'Chemo begin...')
    chemo(prefix, bed)
    chemo2drug(prefix,cancer)
    drug_tips(prefix)
