#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/4/26 
@Description :根据注释工具（annovar/snpEFF)，注释数据库（annovar自带，HGMD，BRCA_Exchange）对vcf文件注释
"""
import os,sys
import time
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


class VcfAnnotation():
    def __init__(self,prefix,threads):
        self.prefix = prefix
        self.threads = threads
        self.path = 'VarientCalling/InterFiles_vcf'
        self.var_path = Profile()
    #annovar自带数据库注释
    @timefly
    def tumor_annovar(self):
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        os.system("perl /data/biosoft/annovar/convert2annovar.pl -format vcf4 `ls VarientCalling/{t}*PASS.vcf` "
                  "-outfile {p}/{t}.avinput -allsample -withfreq -includeinfo >>VarientCalling/{t}.vcf.log 2>&1".format(t=self.prefix,p=self.path))
        os.system("perl /data/biosoft/annovar/table_annovar.pl {p}/{t}.avinput /data/biosoft/annovar/humandb -buildver hg19 "
                  "-out {p}/{t} --maxgenethread 40 --thread {th} "
                  "-remove -protocol refGene,cosmic8N,avsnp150,exac03,EAS.sites.2015_08,clinvar_20190305,intervar_20180118,dbnsfp35a "
                  " -operation g,f,f,f,f,f,f,f -nastring . -otherinfo >>VarientCalling/{t}.vcf.log 2>&1".format(p=self.path,t=self.prefix,th=self.threads))
        os.system("cp {p}/{t}.hg19_multianno.txt VarientCalling/{t}.hg19_multianno.txt".format(p=self.path,t=self.prefix))


    #HGMD遗传数据库注释
    @timefly
    def hgmd_annovar(self):
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        os.system("perl /data/biosoft/annovar/convert2annovar.pl -format vcf4 `ls VarientCalling/{t}*PASS.vcf` "
                  "-outfile {p}/{t}.avinput -allsample -withfreq -includeinfo >>VarientCalling/{t}.vcf.log 2>&1".format(p=self.path,t=self.prefix))
        hgmd_path,hgmd_db_file = self.var_path['hgmd_db_path'].split(',')
        os.system("perl /data/biosoft/annovar/annotate_variation.pl --infoasscore --buildver hg19 --filter --thread {th} "
                  "--maxgenethread 40 --dbtype vcf --vcfdbfile {hf} {p}/{t}.avinput {hp} >>VarientCalling/{t}.vcf.log 2>&1"
                  .format(p=self.path,t=self.prefix,hf=hgmd_db_file,hp=hgmd_path,th=self.threads))
        os.system("cp {p}/{t}.avinput.hg19_vcf_dropped  VarientCalling/{t}.hg19_hgmd.txt ".format(p=self.path,t=self.prefix))

    # BRCA_Exchange数据库注释
    @timefly
    def brca_annovar(self):
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        os.system("perl /data/biosoft/annovar/annotate_variation.pl --infoasscore --buildver hg19 --filter "
                  "--thread {th} --maxgenethread 40 --dbtype vcf --vcfdbfile BRCA_Exchange_db.vcf {p}/{t}.avinput "
                  "/data/DataBase/BRCA_Exchange/ >>VarientCalling/{t}.vcf.log 2>&1".format(t=self.prefix, p=self.path,th=self.threads))
        os.system("mv {p}/{t}.avinput.hg19_vcf_dropped  VarientCalling/{t}.hg19_brca.txt".format(t=self.prefix,p=self.path))

    @timefly
    def snpEff_Anno(self):
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        os.system("java -Xmx10G -jar /data/biosoft/snpEff/snpEff.jar eff -ud 500 -s {p}/{t}.snpEff.html -csvStats {p}/{t}.snpEff.csv "
                  "-interval /data/DataBase/HGMD2019.1/hgmd_pro_2019.1_hg19.vcf -interval /data/DataBase/COSMIC/CosmicCodingMuts.vcf "
                  "-v human_hg19 `ls VarientCalling/{t}*PASS.vcf` 1> VarientCalling/{t}.snpEff.vcf  2>>VarientCalling/{t}.vcf.log ".format(p=self.path,t=self.prefix))


if __name__  == '__main__':
    if len(sys.argv) < 2:
        print('\nusage:  python {} [prefix] [threads]\n'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    threads = sys.argv[2]
    print(time.ctime(),'VcfAnnotation begin...')
    x = VcfAnnotation(prefix,threads)
    x.brca_annovar()