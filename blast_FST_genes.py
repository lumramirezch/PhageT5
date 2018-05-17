#List all proteins from T5
filename = "T5_st0"
folder = "D:\\2016_2020-Doctorat-RAMIREZ_Luis\\001-Bioinformatic_analysis\\Blast\\"
gb_filename = "%s%s.gb" %(folder, filename)
faa_filename = "%s%s_proteins.faa" %(folder, filename)
output_prot = "%s%s.csv" %(folder, filename)
csv_file = "%s%s_statistics.csv" %(folder, filename)
fig_output = "%s%s_fig.pdf" %(folder, filename)

"""
#Import Biopython module SeqIO
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
blast = 'tblastn -query {} -db nr -remote -out {} -evalue \
1e-3 -matrix BLOSUM80 -outfmt "10 qseqid sseqid evalue bitscore"'.format(faa_filename, output_prot)
os.system(blast)

#To try directly on the command prompt

#tblastn -query D:\2016_2020-Doctorat-RAMIREZ_Luis\001-Bioinformatic_analysis\Blast\FST_proteins.faa -db nr -remote -out D:\2016_2020-Doctorat-RAMIREZ_Luis\001-Bioinformatic_analysis\Blast\FST_proteins.csv -evalue 1e-3 -matrix BLOSUM80 -outfmt "10 qseqid sseqid evalue bitscore"
"""
import pandas as pd

#read to csv
data_csv = pd.read_csv(output_prot, header=None)

#Drop the duplicated data
data_csv = data_csv.drop_duplicates(subset=[0, 1])

#get the counts of each gene
genes = data_csv[0].value_counts()

#Separate the ones that contain A1 or A2
T5_like = data_csv[(data_csv[0]=='A1_004') | (data_csv[0]=='A2_006')]

#save the names
T5_like_dict = dict(T5_like[1].value_counts())

#select the subjects that appear with both A1 and A2, and make a list
T5_like_list = [key for key in T5_like_dict if T5_like_dict[key]> 1]

#retrieve the number of matches in the dataframe for the names appearing in the list 
most_common_genes = data_csv[(data_csv[1].isin(T5_like_list))]

#Plot and save as pdf
ax = most_common_genes[0].value_counts().plot(kind='bar')
ax.get_figure().savefig(fig_output, format='pdf', bbox_inches='tight')
