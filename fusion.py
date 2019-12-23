#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/6/24 
@Description :检测融合基因
'''

import os
import time
import sys
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
def Fusion( Tprefix, bed):
    tbam = os.popen("ls Mapping/{0}*final.bam".format(Tprefix)).read().strip()
    path = 'Fusion/InterFiles'
    if not os.path.exists(path):
        os.makedirs(path)
    os.system("perl {factera_path}/factera.pl -o {p} {tb} {factera_path}/exons.bed  {factera_path}/GRCh37.2bit {b} "
              ">Fusion/{T}.fusion.log 2>&1".format(p=path,tb=tbam, b=bed,T=Tprefix,**var_path))
    fusion = tbam.split('/')[-1].replace('bam','factera.fusions.txt')
    os.system("cp {p}/{f} Result/{t}.fusions_result.xls".format(p=path,f=fusion,t=Tprefix))


if __name__  == '__main__':
    if len(sys.argv) < 3:
        print('\nusage:  python {} [Tprefix] [bed] \n'.format(sys.argv[0]))
        sys.exit(1)
    print(time.ctime(),'FUSION begin')
    Tprefix = sys.argv[1]
    bed = sys.argv[2]
    Fusion(Tprefix, bed)