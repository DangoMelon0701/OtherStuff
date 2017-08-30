# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 09:18:38 2017

@author: aleja_blkf3w7
"""

import numpy as np
import pandas as pd
import math as mt
import matplotlib.pyplot as plt
from scipy import interpolate
import time
#%%
def plot_data(energy,data,save_img=0,name='image'):
    fig,axs = plt.subplots()
    for row in data:
        axs.plot(energy,row)
    axs.grid(linestyle='--')
    if save_img == 1:
        fig.savefig("{}.png".format(name),dpi=1000,bbox_inches='tight')

def read_coef(imput_file,):
    return pd.read_csv(imput_file, sep='\s+',header=None, names=['Energias','Coe_tot','Coe_fot'] )

def interpol(data_x,data_y,value):
    func = interpolate.interp1d(data_x,data_y,fill_value='extrapolate')
    return func(value)

def main(co,d,r,l,xv,denv,den,n):
    coeal=read_coef('coeal.txt')
    coenai=read_coef('coenai2.txt')#para tener hast 20 mv
    Eal=np.array(coeal['Energias'])
    Enai=np.array(coenai['Energias'])
    ual=np.array(coeal['Coe_tot'])
    unai=np.array(coenai['Coe_fot'])
    E=np.unique(np.concatenate((Eal,Enai)))
    e = np.zeros([len(d),len(E)])
    eang = 1.0/(2*np.pi)
    for dn,dist in enumerate(d):
        for En,Et in enumerate(E):
            uventana=interpol(Eal,ual,Et)
            udetector=interpol(Enai,unai,Et)
            ec=np.zeros(n+1)
            ar=np.arange(n+1)
            if co>=0 and co<=r:
                hfi=mt.pi/n
                fi=ar*hfi
                for finum,fival in enumerate(fi):
                    g=(co*mt.cos(fival)+mt.sqrt(r**2-(co*mt.sin(fival))**2))
                    a=mt.atan(g/(dist+l))
                    if dist==0:
                        b=mt.pi/2
                    else:
                        b=mt.atan(g/dist)
                    e1=0
                    e2=0
                    if 0<a:
                        h1=a/n
                        te=ar*h1
                        x=l/np.cos(te)
                        f1=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                        e1=h1*(f1.sum()-(f1[0]+f1[-1])/2)
                    if a<b:
                        h2=(b-a)/n
                        te=a+ar*h2
                        x=g/np.sin(te)-dist/np.cos(te)
                        f2=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                        e2=h2*(f2.sum()-(f2[0]+f2[-1])/2)
                    ec[finum]=e1+e2
                e12=hfi*(ec.sum()-(ec[0]+ec[-1])/2)   
            else:
                hfi=mt.asin(r/co)/n
                fi=ar*hfi
                for finum,fival in enumerate(fi):
                    g=(co*mt.cos(fival)+mt.sqrt(r**2-(co*mt.sin(fival))**2))
                    g2=(co*mt.cos(fival)-mt.sqrt(r**2-(co*mt.sin(fival))**2))
                    b=mt.atan(g/(dist+l))
                    if dist==0:
                        a=mt.pi/2
                        c=a
                    else:
                        a=mt.atan(g2/dist)
                        c=mt.atan(g/dist)
                    e1=0
                    e2=0
                    if a<b:
                        h1=(b-a)/n
                        te=a+ar*h1
                        x=l/np.cos(te)
                        f1=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                        e1=h1*(f1.sum()-(f1[0]+f1[-1])/2)
                    if b<c and a<b:
                        h2=(c-b)/n
                        te=b+ar*h2
                        x=g/np.sin(te)-dist/np.cos(te)
                        f2=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                        e2=h2*(f2.sum()-(f2[0]+f2[-1])/2)
                    if a>b and a<c:
                        h2=(c-a)/n
                        te=a+ar*h2
                        x=g/np.sin(te)-dist/np.cos(te)
                        f2=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                        e2=h2*(f2.sum()-(f2[0]+f2[-1])/2)
                    ec[finum]=e1+e2
                e12=hfi*(ec.sum()-(ec[0]+ec[-1])/2)
            e[dn,En]=e12*eang*100    
                
    return e,E

#%%
if __name__ == '__main__':
    Eped = np.arange(0.1,1.21,0.01)
    d = np.array([0,1,1.5,2])
    start_time = time.time()
    a,b = main(0,d,2.54,5.08,0.0508,2.6984,3.67,128)
    print("--- {} seconds --- \n".format(round(time.time() - start_time,2)))
    c=np.zeros((len(d),len(Eped)))
    for En,E in enumerate(Eped):
        for dn,di in enumerate(d):
            f=a[dn,:]
            c[dn,En]=interpol(b,f,E) 
    plot_data(b,a)
    plot_data(Eped,c)