#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Used for parsing Tecan Infinite M200 time series read
Requirement:
 - Download/unzip the fitderiv folder from "http://swainlab.bio.ed.ac.uk/software/fitderiv/"
   PS Swain, K Stevenson, A Leary, LF Montano-Gutierrez, IBN Clark, J Vogel, and T Pilizota. 
   Inferring time derivatives including growth rates using Gaussian processes Nat Commun 7 (2016) 13766
"""

Folder="C:...\\fitderiv\\"
file="%sGrowth_curves.xlsx" %(Folder)
sheets="20180410"

s="strains name"

G_C= ['B2','C2','D2'] # wells control condition - repressed
G_S= ['E2','F2','G2'] # wells stress condition - induced

Hour_start = 0
Hour_end = 40 # One point every 10 minutes = 6.6 hours

##############################################################################
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
#import fitderiv
from fitderiv import *

def DFt(filename, sheetname, wells):
    data_xls = pd.read_excel(filename, sheetname, header=35)
    data_xls = data_xls.set_index('Cycle Nr.').T
    data_xls['Time [s]'] = data_xls['Time [s]'].astype(float) /3600
    data_xls=data_xls.iloc[Hour_start:Hour_end,:]
    data_xls = data_xls.set_index('Time [s]')
    return data_xls[wells]

GC = pd.concat([DFt(file,sheets,G_C[0]), DFt(file,sheets,G_C[1]), DFt(file,sheets,G_C[2])], axis=1)
GS = pd.concat([DFt(file,sheets,G_S[0]), DFt(file,sheets,G_S[1]), DFt(file,sheets,G_S[2])], axis=1)

#GC.sort_index(inplace=True) ; GS.sort_index(inplace=True) 
GC1 = GC.mean(axis=1) ; GC_sd = GC.std(axis=1)
GS1 = GS.mean(axis=1) ; GS_sd = GS.std(axis=1)

GC_fit = fitderiv(GC1.index.values, GC1).ds 
GS_fit = fitderiv(GS1.index.values, GS1).ds
GC_umax = 'Control - umax = %s' %(GC_fit['max df']) ; GC_A = 'Control - Carrying capacity  =   %s' %(GC_fit['max y'])
GS_umax = 'Stress - umax = %s' %(GS_fit['max df']) ; GS_A = 'Stress - Carrying capacity  =   %s' %(GS_fit['max y'])

fig, ax = plt.subplots(1)
ax.plot(GC1.index, GC1, color='blue', label='Repressed')
ax.plot(GS1.index, GS1, color='red', label='Induced')
ax.fill_between(GC_sd.index, GC1+GC_sd, GC1-GC_sd, facecolor='blue', alpha=0.5)
ax.fill_between(GS_sd.index, GS1+GS_sd, GS1-GS_sd, facecolor='red', alpha=0.5)
ax.set_xlabel('Time(Hours)')
ax.set_ylabel('OD610 $\mu$ $\pm \sigma$')
ax.set_title('%s'%(Title))
ax.legend(loc='upper left')
ax.set_ylim([0,1.0])

print(s)
print(GC_umax.replace('.',',')) ; print(GC_A.replace('.',','))
print(GS_umax.replace('.',',')) ; print(GS_A.replace('.',','))
