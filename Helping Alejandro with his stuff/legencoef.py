#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 13:40:04 2017

@author: DangoMelon0701
"""

import numpy as np

class Funciones(object):
    def __init__(self,nombre,apellido,edad):
        self.name = nombre
        self.lastname = apellido
        self.age = edad
        
    def puto(self):
        print("Sabias que {} es un reverendo puto".format(self.name))

def legencoef(n):
    p0 = np.array([1])
    p1 = np.array([1,0])
    if n==0:
        return p0
    elif n==1:
        return p1
    else:
        for i in range(2,n+1):
            pn = ((2*i-1)*np.append(p1,0)-(i-1)*np.append([0,0],p0))/i
            p0=p1
            p1=pn
        return pn

if __name__ == '__main__':
    a = Funciones('Alejandro','Condori Alv',22)
    a.puto()
    
    b = Funciones('Gerardo','Rivera',21)
    b.puto()