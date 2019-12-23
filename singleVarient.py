#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/1 
@Description :单样本检测变异
'''

import  os,sys
import time
from profile import Profile
import functools
from multiprocessing import Process,Pool

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

var_path = Profile()

#------------单样本胚系变异检测--------------------
@timefly
def HaplotypeCaller(path,prefix,prefix2,bed):
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("{gatk} HaplotypeCaller -R {refGenome} -I Mapping/{T}.bqsr.final.bam -ERC GVCF -L {b}  "
              "--native-pair-hmm-threads 16  -O {p}/{T2}.g.vcf.gz "
              ">>VarientCalling/{T2}.call.log 2>&1".format(b=bed,T=prefix,p=path,T2=prefix2,**var_path))
    os.system("{gatk} GenotypeGVCFs -R {refGenome}  -V {p}/{T2}.g.vcf.gz  -O {p}/{T2}.raw.vcf.gz "
              ">>VarientCalling/{T2}.call.log 2>&1".format(T2=prefix2,p=path,**var_path))

    #硬过滤
@timefly
def HardFilter(path,prefix2):
    #分别提取snp/indel子集
    os.system("{gatk} SelectVariants -V {p}/{T2}.raw.vcf.gz -O {p}/{T2}.raw.snps.vcf.gz -select-type SNP "
              ">>VarientCalling/{T2}.call.log 2>&1".format(T2=prefix2, p=path, **var_path))
    os.system("{gatk} SelectVariants -V {p}/{T2}.raw.vcf.gz -O {p}/{T2}.raw.indels.vcf.gz -select-type INDEL -select-type MIXED "
              ">>VarientCalling/{T2}.call.log 2>&1".format(T2=prefix2, p=path, **var_path))
    #分别对snp/indel硬过滤
    os.system("{gatk} VariantFiltration  -R {refGenome} -O {p}/{T2}.filtered.snps.vcf.gz  -V {p}/{T2}.raw.snps.vcf.gz "
              "--filter-name FilterQual --filter-expression \'QUAL < 30.0\' "
              " --filter-name FilterSOR --filter-expression \'SOR > 3.0\'  "
              "--filter-name FilterQD --filter-expression \'QD < 2.0\' "
              "--filter-name FilterMQ --filter-expression \'MQ < 40.0\' "
              "--filter-name FilterFS --filter-expression \'FS > 60.0\' "
              "--filter-name FilterMQRankSum --filter-expression \'MQRankSum < -3.0\' "
              "--filter-name FilterReadPosRankSum --filter-expression \'ReadPosRankSum < -3.0\' "
              ">>VarientCalling/{T2}.call.log 2>&1".format(T2=prefix2,p=path,**var_path))
    os.system("{gatk} VariantFiltration  -R {refGenome} -O {p}/{T2}.filtered.indels.vcf.gz  -V {p}/{T2}.raw.indels.vcf.gz "
              "--filter-name FilterQual --filter-expression \'QUAL < 30.0\' "
              "--filter-name FilterQD --filter-expression \'QD < 2.0\' "
              "--filter-name FilterFS --filter-expression \'FS >200.0\' "
              "--filter-name FilterReadPosRankSum --filter-expression \'ReadPosRankSum < -20.0\' "
              ">>VarientCalling/{T2}.call.log 2>&1".format(T2=prefix2,p=path,**var_path))
    #合并snp/indel过滤结果
    os.system("{gatk} MergeVcfs -I {p}/{T2}.filtered.snps.vcf.gz -I {p}/{T2}.filtered.indels.vcf.gz -O {p}/{T2}.filtered.all.vcf.gz  "
              ">>VarientCalling/{T2}.call.log 2>&1".format(T2=prefix2,p=path,**var_path))
    os.system("{gatk} SelectVariants -V {p}/{T2}.filtered.all.vcf.gz -O VarientCalling/{T2}.PASS.vcf --exclude-filtered true "
              ">>VarientCalling/{T2}.call.log 2>&1".format(T2=prefix2, p=path, **var_path))

@timefly
def Vqsr(path,prefix2):
    #先对SNP建立模型并过滤
    os.system("{gatk} VariantRecalibrator  -V {p}/{T2}.raw.vcf.gz -R {refGenome} -mode SNP -O {p}/{T2}.vqsr_SNP.recal "
              "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {bvqsrKnownPath}hapmap_3.3.b37.vcf.gz "
              " -resource:omni,known=false,training=true,truth=true,prior=12.0 {bvqsrKnownPath}1000G_omni2.5.b37.vcf.gz "
              "-resource:1000G,known=false,training=true,truth=false,prior=10.0 {bvqsrKnownPath}1000G_phase1.snps.high_confidence.hg19.sites.nochr.vcf.gz "
              " -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {bvqsrKnownPath}dbsnp_138.hg19.nochr.vcf.gz "
              " -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR "
              "--tranches-file {p}/{T2}.vqsr_SNP.tranches --rscript-file {p}/{T2}.vqsr_SNP_plots.R"
              ">>VarientCalling/{T2}.call.log 2>&1".format(T2=prefix2, p=path, **var_path))
    os.system("{gatk} ApplyVQSR  -V {p}/{T2}.raw.vcf.gz  -R {refGenome} -mode SNP -O {p}/{T2}.vqsr_snps_raw_indels.vcf "
              "-ts-filter-level 99.0 --recal-file {p}/{T2}.vqsr_SNP.recal --tranches-file {p}/{T2}.vqsr_SNP.tranches "
              ">>VarientCalling/{T2}.call.log 2>&1".format(T2=prefix2, p=path, **var_path))
    # 后对INDEL建立模型并过滤
    os.system("{gatk} VariantRecalibrator  -V {p}/{T2}.vqsr_snps_raw_indels.vcf -R {refGenome} -mode INDEL -O {p}/{T2}.vqsr_INDEL.recal "
        "-resource:mills,known=false,training=true,truth=true,prior=12.0 {bvqsrKnownPath}Mills_and_1000G_gold_standard.indels.GRCh37.vcf.gz "
        " -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {bvqsrKnownPath}dbsnp_138.hg19.nochr.vcf.gz "
        " -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum --max-gaussians 4 "
        "--tranches-file {p}/{T2}.vqsr_INDEL.tranches --rscript-file {p}/{T2}.vqsr_INDEL_plots.R"
        ">>VarientCalling/{T2}.call.log 2>&1".format(T2=prefix2, p=path, **var_path))
    os.system("{gatk} ApplyVQSR  -V {p}/{T2}.vqsr_snps_raw_indels.vcf  -R {refGenome} -mode INDEL -O {p}/{T2}.vqsr.vcf "
        "-ts-filter-level 99.0 --recal-file {p}/{T2}.vqsr_INDEL.recal --tranches-file {p}/{T2}.vqsr_INDEL.tranches "
        ">>VarientCalling/{T2}.call.log 2>&1".format(T2=prefix2, p=path, **var_path))
    os.system("{gatk} SelectVariants -V {p}/{T2}.vqsr.vcf -O VarientCalling/{T2}.PASS.vcf --exclude-filtered true "
        ">>VarientCalling/{T2}.call.log 2>&1".format(T2=prefix2, p=path, **var_path))

def SingleNormalVarient(prefix,bed,panel):
    prefix2 = prefix + '.germline'
    path = 'VarientCalling/InterFiles_call'
    HaplotypeCaller(path,prefix,prefix2,bed)
    if panel == 'exon':
        Vqsr(path,prefix2)
    else:
        HardFilter(path,prefix2)



#------------单样本体细胞变异检测--------------------
def PON(Tprefix, bed, panel):
    TNprefix = Tprefix + '.somatic'
    path = 'VarientCalling/InterFiles_call'
    if not os.path.exists(path):
        os.makedirs(path)
    pon_path = '{0}_{1}'.format(var_path['PonPath'], panel)  # 不同panel基线目录
    if not os.path.exists(pon_path):
        os.makedirs(pon_path)
    vcf_files = os.popen(" ls %s/*pon.vcf.gz.tbi |while read i;do echo -V ${i/.tbi/};done" % pon_path).read().replace('\n', ' ')
    os.system("{gatk} GenomicsDBImport -R {refGenome} -L {b} --genomicsdb-workspace-path pon_db {v} >>VarientCalling/{TN}.call.log 2>&1"
        .format(b=bed, v=vcf_files, TN=TNprefix, **var_path))
    os.system("{gatk} CreateSomaticPanelOfNormals -R {refGenome} -V gendb://pon_db -O {p}/pon.vcf.gz >>VarientCalling/{TN}.call.log 2>&1"
        .format(v=vcf_files,  p=path, TN=TNprefix, **var_path))
    os.system("rm -rf pon_db")

def Mutect2(Tprefix, bed, splitNum):
    path = 'VarientCalling/InterFiles_call'
    TNprefix = Tprefix  + '.somatic'
    path2 = '{0}/{1}'.format(path, TNprefix)
    tbam = os.popen("ls Mapping/{0}*final.bam".format(Tprefix)).read().strip()
    os.system("{gatk} SplitIntervals -R {refGenome} -L {b} -scatter {num} -O {p2}_bedpath  >>VarientCalling/{TN}.call.log 2>&1"
              .format(TN=TNprefix, b=bed, p2=path2,num=splitNum, **var_path))
    bedpath = path2 + '_bedpath'
    vcfpath = path2 + '_vcfpath'
    if not os.path.exists(vcfpath):
        os.makedirs(vcfpath)
    proc_list = []
    for bed in os.listdir(bedpath):
        pre = bed.split('-')[0]
        bed2 = os.path.join(bedpath, bed)
        proc = Process(target=Split4Mutect2,args=(path,Tprefix, TNprefix, tbam, bed2, pre,))
        proc_list.append(proc)
    for p in proc_list:
        p.start()
    for p in proc_list:
        p.join()
    vcflist = [os.path.join(vcfpath,vcf) for vcf in os.listdir(vcfpath) if vcf.endswith('vcf.gz')]
    statslist = [os.path.join(vcfpath,vcf) for vcf in os.listdir(vcfpath) if vcf.endswith('stats')]
    vcflist2 = ' -I '.join(vcflist)
    statslist2 = ' -stats '.join(statslist)
    os.system("{gatk} MergeVcfs -I {vl} -O  {p}/{TN}.raw.vcf.gz  >>VarientCalling/{TN}.call.log 2>&1".format(TN=TNprefix, p=path,vl=vcflist2, **var_path))
    os.system("{gatk} MergeMutectStats -stats {stats} -O  {p}/{TN}.stats  >>VarientCalling/{TN}.call.log 2>&1".format(TN=TNprefix, p=path,stats=statslist2, **var_path))
    os.system("{gatk} FilterMutectCalls -V {p}/{TN}.raw.vcf.gz -O {p}/{TN}.vcf.gz -R {refGenome} -stats {p}/{TN}.stats "
              #"--max-events-in-region 5 --max-alt-allele-count 5 --min-slippage-length 10 --long-indel-length 10 "
              ">>VarientCalling/{TN}.call.log 2>&1".format(TN=TNprefix, p=path, **var_path))
    os.system("{gatk} SelectVariants -V {p}/{TN}.vcf.gz -O VarientCalling/{TN}.PASS.vcf "
              "--exclude-filtered true >>VarientCalling/{TN}.call.log 2>&1".format(TN=TNprefix, p=path, **var_path))

def Split4Mutect2(path,Tprefix, TNprefix,tbam, bed, pre):
    vcfpath = '{0}/{1}_vcfpath'.format(path, TNprefix)
    os.system("{gatk} Mutect2 -R {refGenome}  -I {tb}  -L {b} --germline-resource {af-only-gnomad} --max-reads-per-alignment-start 0 "
        "--af-of-alleles-not-in-resource 0.0000025 -pon {p}/pon.vcf.gz -O {vp}/{n}.raw.vcf.gz >>VarientCalling/{TN}.call.log 2>&1"
        .format(tb=tbam, T=Tprefix, p=path, vp=vcfpath, TN=TNprefix, b=bed, n=pre, **var_path)) #--max-reads-per-alignment-start 0

@timefly
# 利用基线库对肿瘤样本变异检测
def SingleTumorVarient(Tprefix,bed, panel,splitNum):
    PON(Tprefix, bed, panel)
    Mutect2(Tprefix, bed, splitNum)



if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('\nusage:  python {} [prefix] [bed]  [panel] [splitNum]\n'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    bed = sys.argv[2]
    panel = sys.argv[3]
    splitNum = sys.argv[4]
    print(time.ctime(), 'SingleVarientCalling begin')
    SingleNormalVarient(prefix,bed,panel)
    # SingleTumorVarient(prefix, bed, panel, splitNum)

