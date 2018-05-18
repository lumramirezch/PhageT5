#List all proteins from T5
filename = "T5_st0"
folder = " " #write your directory path here
gb_filename = "%s%s.gb" %(folder, filename)
faa_filename = "%s%s_proteins.faa" %(folder, filename)
output_prot = "%s%s.csv" %(folder, filename)
csv_file = "%s%s_statistics.csv" %(folder, filename)
fig_output1 = "%s%s_fig1.pdf" %(folder, filename)
fig_output2 = "%s%s_fig2.pdf" %(folder, filename)

Import Biopython module SeqIO
from Bio import SeqIO
input_handle = open(gb_filename, 'r')
output_handle = open(faa_filename, 'w')

#Parse the genbank file 
for seq_record in SeqIO.parse(input_handle, 'genbank'):
    for seq_features in seq_record.features:
        if seq_features.type == "CDS":
            assert len(seq_features.qualifiers['translation']) == 1
            if "locus_tag" not in seq_features.qualifiers:
                print(seq_features.qualifiers['protein_id'])
                output_handle.write(">%s\n%s\n" %(
                   seq_features.qualifiers['protein_id'][0],
                   seq_features.qualifiers['translation'][0]))
            elif "gene" in seq_features.qualifiers:
                a = "gene"
            else:
                a = "locus_tag"
            output_handle.write(">%s_%s\n%s\n" %(
                   seq_features.qualifiers[a][0],
                   seq_features.qualifiers['locus_tag'][0][-3:],
                   seq_features.qualifiers['translation'][0]))
            
input_handle.close()
output_handle.close()

#blast against the database by using a -remote flag
import os
blast = 'tblastn -query {} -db nr -remote -out {} -evalue 1e-3 -matrix BLOSUM80 -outfmt "10 qseqid sseqid evalue bitscore"'.format(faa_filename, output_prot)
os.system(blast)

import pandas as pd
import matplotlib.pyplot as plt

data_csv = pd.read_csv(output_prot, header=None) #read to csv
data_csv = data_csv.drop_duplicates(subset=[0, 1]) #Drop the duplicated data
FST_genes_counts = data_csv[0].value_counts().to_dict() #get the counts of each gene
FST_genes_list= [key for key in dict(data_csv[0].value_counts())] #get a sorted list of genes
FST_genes_list = sorted(FST_genes_list, key=lambda qey: qey.split('_')[1])
#labellist = [i[:-4] for i in FST_genes_list]
labellist = [x.split('_')[0][-3::] for x in FST_genes_list]

w = ((data_csv[0] == 'A1_004') | (data_csv[0] == 'A2_006')) #Entries matching genes A1 and A2

T5_like_dict = dict((data_csv[w])[1].value_counts()) #Separate the and save the names of those entries
T5_virus = [key for key in T5_like_dict if T5_like_dict[key]> 1] #list the subjects with A1 and A2
T5virus_genes = data_csv[(data_csv[1].isin(T5_virus))]

#T5virus_genes = data_csv[(data_csv[1].isin(T5_virus) & ~(w))]
#genes that appear altogether with A1 and A2
#T5virus_genes = data_csv[(data_csv[1].isin(T5_virus) & ~((data_csv[0] == 'A1_004') | (data_csv[0] == 'A2_006')))]
T5virus_genes_dict = dict(T5virus_genes[0].value_counts())

#Plot and save as pdf
y_pos = [y for y in range(len(FST_genes_list))]
x_pos = []
x_pos2 = []

#Individual Frequency
fig, ax = plt.subplots()
for item in FST_genes_list:
    x_pos2.append(FST_genes_counts[item])
d = ['red' if ((x=='A1_004')|(x=='A2_006')) else 'blue' for x in FST_genes_list]
ax.barh(y_pos, x_pos2, align='center', color = d)
ax.set_yticks(y_pos)
ax.set_yticklabels(labellist)
ax.invert_yaxis()
ax.set_title("FST genes frequency of matches")
plt.ylabel('FST Genes')
plt.show()
#ax.get_figure().savefig(fig_output2, format='pdf', bbox_inches='tight')

#Frequency altogether with A1 and A2
fig, ax = plt.subplots()
for item in FST_genes_list:
    x_pos.append(T5virus_genes_dict[item])
    
c = ['red' if ((x=='A1_004')|(x=='A2_006')) else ('green' if (((x=='T5p002_002')|(x=='T5p005_005')|(x=='T5p007_007')|(x=='T5p010_010')|(x=='T5p013_013')|(x=='T5p014_014'))) else 'blue' )for x in FST_genes_list]
#c = ['red' if (x==max(x_pos)) else ('green' if (x>20) else 'blue') for x in x_pos]
ax.barh(y_pos, x_pos, align='center', color= c)
ax.set_yticks(y_pos)
ax.set_yticklabels(labellist)
ax.invert_yaxis()
ax.set_title("FST genes frequency in T5 viruses")
plt.ylabel('FST Genes')
plt.show()
#ax.get_figure().savefig(fig_output1, format='pdf', bbox_inches='tight')
