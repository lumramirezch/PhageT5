# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 10:36:47 2017

@author: lm.ramirez-chamorro
"""



#Folder Name, remember to add always two slashes
Folder="C:\\Users\\leo\\Desktop\\Thèse\\Résultats\\TECAN\\"
#File Name, without the extension
filename1 = '230518_A2'
filename2 = None
filename3 = None
strain = 'BL21' #Construction
#Range to be visualized
Hour_start = 0
Hour_end = 40
#List of wells in the file 1, Repressed and Induced
file1_Repr = ['B2', 'B3', 'B4', 'B5', 'B6']  #Repressed
file1_Ind = ['C2', 'C3', 'C4', 'C5', 'C6']    #Induced
#Change if not the same arrangement		

##############################################################################
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def Data_frame_generator(filename, folder):
    data_xls = pd.read_excel('%s%s.xlsx' %(Folder, filename), 'Sheet0', header=35)
    data_xls = data_xls.set_index('Cycle Nr.').T
    data_xls['Time [s]'] = data_xls['Time [s]'].astype(float) /3600
    data_xls=data_xls.iloc[Hour_start:Hour_end,:]
    #data_xls = data_xls.set_index('Time [s]')
    return data_xls

csvT1 = Data_frame_generator(filename1, Folder)
if filename2 is not None:
	file2_Repr = file1_Repr
	file2_Ind = file1_Ind
	csvT2 = Data_frame_generator(filename2, Folder)
	if filename3 is not None:
		file3_Repr = file1_Repr
		file3_Ind = file1_Ind
		csvT3 = Data_frame_generator(filename3, Folder)
		Repr = pd.concat([csvT1[file1_Repr], csvT2[file2_Repr], csvT3[file3_Repr]], axis=1, ignore_index=True)
		Ind = pd.concat([csvT1[file1_Ind], csvT2[file2_Ind], csvT3[file3_Ind]], axis=1, ignore_index=True)
		Hoursi = pd.concat([csvT1['Time [s]'], csvT2['Time [s]'], csvT3['Time [s]']], axis=1, ignore_index=True)
	else:
		Repr = pd.concat([csvT1[file1_Repr], csvT2[file2_Repr]], axis=1, ignore_index=True)
		Ind = pd.concat([csvT1[file1_Ind], csvT2[file2_Ind]], axis=1, ignore_index=True)
		Hoursi = pd.concat([csvT1['Time [s]'], csvT2['Time [s]']], axis=1, ignore_index=True)
else:
	Repr = pd.concat([csvT1[file1_Repr]], axis=1, ignore_index=True)
	Ind = pd.concat([csvT1[file1_Ind]], axis=1, ignore_index=True)
	Hoursi = pd.concat([csvT1['Time [s]']], axis=1, ignore_index=True)	


Hoursi = Hoursi[Hour_start:Hour_end]
Repr = Repr[Hour_start:Hour_end]
Ind = Ind[Hour_start:Hour_end]

Hours = Hoursi.mean(axis=1) 
Repr_mean = Repr.mean(axis=1)
Ind_mean = Ind.mean(axis=1)
Repr_sd = Repr.std(axis=1)
Ind_sd = Ind.std(axis=1)

fig, ax = plt.subplots(1)
ax.plot(Hours, Repr_mean, color='blue')
ax.plot(Hours, Ind_mean, color='red')
ax.fill_between(Hours, Repr_mean+Repr_sd, Repr_mean-Repr_sd, facecolor='blue', alpha=0.5)
ax.fill_between(Hours, Ind_mean+Ind_sd, Ind_mean-Ind_sd, facecolor='red', alpha=0.5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

white = mpatches.Patch(color='white', label=strain)
blue = mpatches.Patch(color='blue', label='repressed')
red = mpatches.Patch(color='red', label='induced')
plt.legend(handles=[white,blue,red], loc=2, borderpad=None)
ax.set_xlabel('Time(Hours)')
ax.set_ylabel('OD600')
ax.set_ylim([0,0.8]) ; ax.set_xlim([0,6.5])

plt.show()
"""fig_output = "%s_%s_%s_%s-%s.pdf"%(Folder, filename1[:10], filename2[4:10], filename3[4:10], strain)
ax.get_figure().savefig(fig_output, format='pdf', bbox_inches='tight')"""
