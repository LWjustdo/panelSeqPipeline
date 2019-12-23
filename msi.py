#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/6/24 
@Description :检测MSI
'''
import os
import time
import sys
from profile import  Profile
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
def MSI(Tprefix, Nprefix, bed):
    TNprefix = Tprefix + '_vs_' + Nprefix + '.somatic'
    tbam = os.popen("ls Mapping/{0}*final.bam".format(Tprefix)).read().strip()
    nbam = os.popen("ls Mapping/{0}*final.bam".format(Nprefix)).read().strip()
    path = 'MSI/InterFiles_all'
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("{msi_path} msi -d {MSIpos_list} -e {b} -n {nbam} -t {tbam} -o {p}/{tn} "
              ">MSI/{tn}.all.log 2>&1".format(b=bed,nbam=nbam,tbam=tbam,tn=TNprefix,p=path,**var_path))
    os.chdir(path)
    os.system("cp {0} ../../Result/{0}.MSI_all_result.xls".format(TNprefix))
    os.system("cut -f 7 {0}_somatic >{0}.msi.score.list && "
              "perl /home/longzhao/panelSel/src/convert4MSI.pl {0}.msi.score.list > {0}.msi.status.score.list && "
              "Rscript {msi_plotR} {0}.msi.score.list ../{0}_msiScore.pdf  2>> ../{0}.all.log && "
              "Rscript {msi_plotR} {0}.msi.status.score.list {0}_msiScoreStatus.pdf 2>> ../{0}.all.log".format(TNprefix,**var_path))
    os.chdir('../../')



if __name__  == '__main__':
    if len(sys.argv) < 4:
        print('\nusage:  python {} [Tprefix] [Nprefix] [bed] \n'.format(sys.argv[0]))
        sys.exit(1)
    print(time.ctime(),'MSI begin')
    Tprefix = sys.argv[1]
    Nprefix = sys.argv[2]
    bed = sys.argv[3]
    MSI(Tprefix, Nprefix, bed)