
import pandas as pd
from Bio import SeqIO
import numpy as np
import plotly.graph_objects as go
import os.path
import plotly.offline as pio
import re
import sys

os.chdir(os.path.dirname(os.path.realpath(sys.argv[0])))

ref_taxid = "10760" #T5
#ref_taxid = "10726" #T5
#ref_taxid = "10665" #T4

filen = ""
file1 = "%s\\Orthogroups%s.tsv" % (os.getcwd(), filen)
file2 = "%s\\Orthogroups2.csv" % (os.getcwd())
genecounts = "%s\\Orthogroups.GeneCount%s.tsv" % (os.getcwd(), filen)
outputC="%s\\Orthogroups-Circular%s.html" % (os.getcwd(), filen)
outputB="%s\\Orthogroups-Bar%s.html" % (os.getcwd(), filen)
filegb = "%s\\reference.gb" % (os.getcwd())
essentials = '%s\\essential_genes.txt' % (os.getcwd())


df = pd.read_csv(file1, sep='\t', header=0)
df = df.rename(columns=lambda x: re.sub('\r','',x))

df2 = pd.read_csv(genecounts, sep='\t', header=0)
df2 = df2.rename(columns=lambda x: re.sub('\r','',x))
df2 = df2[['Orthogroup','Total']]

Refgenes = {}
Reflocus_tag = {}

def to_db(file_name):
    input_handle = open(file_name, 'r')
    for seq_record in SeqIO.parse(input_handle, 'genbank'):
        for seq_features in seq_record.features:
            if seq_features.type == "CDS":
                try:
                    assert len(seq_features.qualifiers['translation']) == 1
                    dictkey = seq_features.qualifiers['protein_id'][0]
                    if "locus_tag" in seq_features.qualifiers:
                        Reflocus_tag[dictkey] = seq_features.qualifiers['locus_tag'][0]
                        if "gene" in seq_features.qualifiers:
                            Refgenes[dictkey] = seq_features.qualifiers['gene'][0]
                    else:
                        Reflocus_tag[dictkey] = dictkey
                        Refgenes[dictkey] = dictkey
                except KeyError:
                    continue
    input_handle.close()

to_db(filegb)

pat = '|'.join(r"\b{}\b".format(x) for x in Reflocus_tag.keys())

df.insert(0, 'locus_tag', df[ref_taxid].str.extract('('+ pat + ')', expand=False).map(Reflocus_tag))
df.insert(1, 'gene', df[ref_taxid].str.extract('('+ pat + ')', expand=False).map(Refgenes))
df = df.sort_values(["locus_tag"])
df.reset_index(inplace=True, drop=True)

df3 = df.merge(df2, left_on='Orthogroup', right_on='Orthogroup')

move_cols = ['locus_tag', 'gene', 'Total', 'Orthogroup']
new_cols = np.hstack((move_cols, df.columns.difference(move_cols)))
df3 = df3.loc[:, new_cols]
df3['Total'] = df3.apply(lambda x: df3.shape[1]- 4 - x.isnull().sum(), axis='columns')
#df3['Total'] = df3['Total']/df3['Total'].max()
df3.to_csv(file2)

df3['gene'] = df3['gene'].fillna(df3['locus_tag'])
df4 =df3[['gene','Total']].dropna()
x=df4['gene'].values.tolist()
y=df4['Total'].values.tolist()

df5 = pd.read_csv(essentials, sep="\t")
df5_list = list(df5['Gene'])
#df4['color'] = df['gene'].apply(lambda x: 'green' if x in df5_list else 'gray')
df4['color'] = df['locus_tag'].apply(lambda x: 'green' if x in df5_list else 'gray')
color=df4['color'].values.tolist()


layout = go.Layout(yaxis=dict(type='category'), plot_bgcolor='white', autosize=False, width=1000, height=2000)
fig = go.Figure(data=[go.Bar(x=y, y=x, marker_color=color, orientation='h')], layout=layout)
pio.plot(fig, filename=outputB)

fig = go.Figure()
fig.add_trace(go.Barpolar(
    r=df4['Total'],
    theta=360*df4.index.values/len(df4['Total']),
    width=360/len(df4['Total']),
    marker_color=color,
    opacity=0.8,
    base=100, 
))

fig.update_layout(
    template=None,
    polar = dict(
        radialaxis = dict(showticklabels=False, ticks=''),
        angularaxis = dict(showticklabels=False, ticks='', rotation=90, direction = "clockwise")
    )
)
pio.plot(fig, filename=outputC)
