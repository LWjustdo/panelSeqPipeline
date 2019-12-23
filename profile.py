#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/1 
@Description : Profile.txt文件
'''
import os

#脚本所在路径
path = os.path.dirname(__file__)
file = path +'/Profile.txt'

def Profile():
    dic = {}
    with open(file) as f:
        for line in f:
            lin = line.strip().split(':')
            if len(lin) == 2:
                dic[lin[0]] = lin[1]
    return dic

if __name__ == '__main__':
    var = Profile()
    print(var)
