#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/11 
@Description :
'''
import sys,os
import time
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
def MSI5Points(Tprefix, Nprefix, threads):
    TNprefix = Tprefix + '_vs_' + Nprefix + '.somatic'
    tbam = os.popen("ls Mapping/{0}*final.bam".format(Tprefix)).read().strip()
    nbam = os.popen("ls Mapping/{0}*final.bam".format(Nprefix)).read().strip()
    path = 'MSI/InterFiles_5points'
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("python {mantis} -b {MSI_Five} --genome {refGenome} -n {nb} -t {tb} --threads {th} -o {p}/{tn}.5points "
              ">MSI/{tn}.5points.log 2>&1".format(nb=nbam,tb=tbam,th=threads,p=path,tn=TNprefix,**var_path))
    os.system("cp {p}/{tn}.5points.status Result/{tn}.MSI_5points_result.xls".format(p=path,tn=TNprefix))


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("\n Usage: python {} [Tprefix] [Nprefix]  [threads]\n".format(sys.argv[0]))
        sys.exit(1)
    print(time.ctime(), 'MSI5Points begin')
    Tprefix = sys.argv[1]
    Nprefix = sys.argv[2]
    threads = sys.argv[3]
    MSI5Points(Tprefix, Nprefix, threads)
