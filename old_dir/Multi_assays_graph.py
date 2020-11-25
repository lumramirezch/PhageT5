# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 10:36:47 2017

@authors: lm.ramirez-chamorro and Léo Zangelmi
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import os
import sys

#Folder Name, remember to add always two slashes
Folder="P:\\Equipes\\T5PHAG\\Programmes Python\\Tecan\\"
subfolder = input("Subfolder name? ")+"\\"
#File Name, without the extension
filename1 = input("Excel File n°1 name? (without .xlsx) ")
#Check if the file exists
my_file = Path("%sexcel files\\%s.xlsx"%(Folder,filename1))
if my_file.is_file() == False:
	print("Error: the file does not exist. Please check the spelling.")
	os.system("pause")
	sys.exit()
filename2 = input("Excel File n°2 name? (leave blank if no file 2) ")
#Check if the file exists
my_file = Path("%sexcel files\\%s.xlsx"%(Folder,filename2))
if my_file.is_file() == False and filename2 is not "":
	print("Error: the file does not exist. Please check the spelling.")
	os.system("pause")
	sys.exit()
#We don't want the user to give a file n°3 if he didn't give a file n°2
if filename2 is not "":
	filename3 = input("Excel File n°3 name? (leave blank if no file 3) ")
	#Check if the file exists
	my_file = Path("%sexcel files\\%s.xlsx"%(Folder,filename3))
	if my_file.is_file() == False and filename3 is not "":
		print("Error: the file does not exist. Please check the spelling.")
		os.system("pause")
		sys.exit()
strain = input("Name of the main label? ") #Construction
file = input("Name of the final file? ")
#Range to be visualized
Hour_start = 0
Hour_end = int(input("How many time points for x axis? "))
#List of wells in the file 1, Repressed and Induced
file1_Repr = input("List of wells with repressor (All files will use the same wells. Example format: A1,A2,A3) ").split(",") #split(",") converts string separated by a comma to list
file1_Ind = input("List of wells with inductor (All files will use the same wells. Example format: A1,A2,A3) ").split(",")

##############################################################################

def Data_frame_generator(filename, folder):
    data_xls = pd.read_excel('%sexcel files\\%s.xlsx' %(Folder, filename), 'Sheet0', header=35)
    data_xls = data_xls.set_index('Cycle Nr.').T
    data_xls['Time [s]'] = data_xls['Time [s]'].astype(float) /3600
    data_xls=data_xls.iloc[Hour_start:Hour_end,:]
    #data_xls = data_xls.set_index('Time [s]')
    return data_xls

csvT1 = Data_frame_generator(filename1, Folder)
if filename2 is not "":
	file2_Repr = file1_Repr
	file2_Ind = file1_Ind
	csvT2 = Data_frame_generator(filename2, Folder)
	if filename3 is not "":
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
ax.set_ylim([0,0.8]) ; ax.set_xlim([0,Hour_end//6])

#plt.show()
fig_output = "%sresults\\%s%s.pdf"%(Folder, subfolder, file)

#If the folder does not exist, we create it
if not os.path.exists("%sresults\\%s"%(Folder, subfolder)):
    os.makedirs("%sresults\\%s"%(Folder, subfolder))

my_file = Path(fig_output)
if my_file.is_file():
    # exists
	overwrite = input("File already exists. Overwrite? (y/n). Warning: if the file is open, don't forget to close it! ")
	if overwrite == "y" or overwrite == "Y":
		ax.get_figure().savefig(fig_output, format='pdf', bbox_inches='tight')
else:
    # doesn't exist

