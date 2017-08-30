#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:57:09 2017

@author: DangoMelon0701
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
import os,time,math
#%%
def plot_data(energy,data,save_img=0,name='image'):
    fig,axs = plt.subplots()
    for row in data:
        axs.plot(energy,row)
    axs.grid(linestyle='--')
    if save_img == 1:
        fig.savefig("{}.png".format(name),dpi=1000,bbox_inches='tight')
        
def read_coefs(input_file):
    return pd.read_csv(input_file,header=None, sep='\s+',names=['Energia','Coef_tot','Coef_fot'])

def interpol(data_frame,value,column='Coef_fot'):
    if column == 2:
        column = 'Coef_tot'
    func = interpolate.interp1d(data_frame['Energia'],data_frame[column],fill_value='extrapolate')
    return func(value)

def trapecio(np_array,number):
    return number*(np_array.sum()-(np_array[0]+np_array[-1])/2.)

def read_data():
    for files in os.listdir(os.getcwd()):
        if 'coeal' in files:
            coeal = read_coefs(files)
        elif 'coecar' in files:
            coecar = read_coefs(files)
        elif 'coenai' in files:
            coenai = read_coefs(files)
    return coeal,coecar,coenai

def main(co,E,d,r,l,xv,denv,den,n):
    coeal,coecar,coenai = read_data()
    e = np.zeros([len(d),len(E)])
    eang = 1.0/(2*np.pi)
    start_time = time.time()
    for dnum,dist in enumerate(d.astype(np.float)):
        for Enum,Et in enumerate(E):
            uventana=interpol(coeal,Et,2)
            udetector=interpol(coenai,Et)
            ec = np.zeros([n+1])
            fi = np.arange(n+1)
            
            if co >= 0 and co <=r:
                hfi = np.pi/n
                fi = fi* hfi
                for num,q in enumerate(fi):
                    g = (co*math.cos(q)+(r**2-(co*math.sin(q))**2)**0.5)
                    a = math.atan(g/(dist+l))
                    if dist == 0:
                        b = math.atan(np.inf)
                    else:
                        b = math.atan(g/dist)
                    e1 = 0
                    e2 = 0
                    if 0 < a:
                        h1=a/n
                        te = np.arange(n+1)*h1
                        x = np.divide(l,np.cos(te))
                        f1 =(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(np.divide(-uventana*denv*xv,np.cos(te)))
                        e1 = trapecio(f1,h1)
                    if a < b:
                        h2=(b-a)/n
                        te = a+np.arange(n+1)*h2
                        x = np.divide(g,np.sin(te))-np.divide(dist,np.cos(te))
                        f2 =(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(np.divide(-uventana*denv*xv,np.cos(te)))
                        e2 = trapecio(f2,h2)
                    ec[num] = e1+e2
                e12 = trapecio(ec,hfi)
            
            else:
                hfi =math.acos(r/co)/n
                fi = fi * hfi
                for num,q in enumerate(fi):
                    g = (co*math.cos(q)+(r**2-(co*math.sin(q))**2)**0.5)
                    g2 = (co*math.cos(q)-(r**2-(co*np.sin(q))**2)**0.5)
                    a = math.atan(g2/dist)
                    b = math.atan(g/(dist+l))
                    c = math.atan(g/dist)
                    e1 = 0
                    e2 = 0
                    if a < b:
                        h1 = (b-a)/n
                        te = a+np.arange(n+1)*h1
                        x = np.divide(l,np.cos(te))
                        f1 =(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(np.divide(-uventana*denv*xv,np.cos(te)))
                        e1 = trapecio(f1,h1)
                    if b < c and a < b:
                        h2 = (c-b)/n
                        te = b+np.arange(n+1)*h2
                        x = np.divide(g,np.sin(te))-np.divide(dist,np.cos(te))
                        f2 =(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(np.divide(-uventana*denv*xv,np.cos(te)))
                        e2 = trapecio(f2,h2)
                    if a > b and a<c:
                        h2 = (c-a)/n
                        te = b+np.arange(n+1)*h2
                        x = np.divide(g,np.sin(te))-np.divide(dist,np.cos(te))
                        f2 =(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(np.divide(-uventana*denv*xv,np.cos(te)))
                        e2 = trapecio(f2,h2)
                    ec[num] = e1+e2
                e12 = trapecio(ec,hfi)
            e[dnum][Enum]= e12*eang*100
    print("--- {} seconds --- \n".format(round(time.time() - start_time,2)))
    return e

#%%
if __name__ == '__main__':
    E = np.arange(0.1,1.21,0.01)
    d = np.array([0,1,2,3])
    a = main(0,E,d,2.54,5.08,0.0508,2.6984,3.67,128)
    plot_data(E,a)