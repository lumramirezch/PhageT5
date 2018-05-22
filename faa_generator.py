#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
import pandas as pd
from Bio import Entrez, SeqIO
import argparse
import urllib

parser = argparse.ArgumentParser()
parser.add_argument("-g", dest = "fg") # gb file to parse and write many faa files
parser.add_argument("-s", dest = "fcsv") #
parser.add_argument("-gb", dest = "fgb")
parser.add_argument("-fa", dest ="ffaa")

args = parser.parse_args()

filegb = args.fg
filecsv = args.fcsv
folder_gb = args.fgb 
folder_faa = args.ffaa

def to_db(file_name):
    if filecsv is not None:
        output_handle = open("%s%s.faa" %(folder_faa, file_name), 'a+')
        input_handle = open("%s%s.gb" %(folder_gb, file_name), 'r')
    elif filecsv is None:
        input_handle = open("%s%s" %(folder_gb, file_name), 'r')
    for seq_record in SeqIO.parse(input_handle, 'genbank'):
        for seq_features in seq_record.features:
            if seq_features.type == "CDS":
                try:
                    assert len(seq_features.qualifiers['translation']) == 1
                    if filegb is not None:
                        output_handle = open("%s%s-%s.faa" %(folder_faa, file_name[:-3], seq_features.qualifiers['locus_tag'][0]), 'a+')
                        if "gene" in seq_features.qualifiers:
                            a = seq_features.qualifiers['gene'][0]
                        elif "locus_tag" in seq_features.qualifiers:
                            a = seq_features.qualifiers['locus_tag'][0]
                        else:
                            a = seq_features.qualifiers['protein_id'][0]
                        output_handle.write(">%s|[%s]\n%s\n" %(a, file_name, seq_features.qualifiers['translation'][0]))
                        output_handle.close()
                    elif filegb is None:
                        if "gene" in seq_features.qualifiers:
                            a = seq_features.qualifiers['gene'][0]
                        elif "locus_tag" in seq_features.qualifiers:
                            a = seq_features.qualifiers['locus_tag'][0]
                        else:
                            a = seq_features.qualifiers['protein_id'][0]
                        output_handle.write(">%s|[%s]\n%s\n" %(a, file_name, seq_features.qualifiers['translation'][0]))
                except KeyError:
                    continue
    input_handle.close()
    output_handle.close()
   

Entrez.email = "enteryourmail@here.com"
if filecsv is not None:
    simseqs = pd.read_csv(filecsv, header=None, names=['sseqid','qseqid','evalue','bitscore','qcovs','pident'])
    simseqs.head()
    simseqs = simseqs.drop_duplicates(subset=['sseqid'])
    simseqs = simseqs.reset_index(drop=True)
    simseqs_list = []
    for item in list(simseqs['sseqid']):
        try:
            simseqs_list.append(item.split("|")[3])
        except IndexError:
            simseqs_list.append(item)
        continue
    print("Genomes to be downloaded: %s"%(len(simseqs_list)))
    Entrez.email = "enteryourmail@here.com"
    for item in simseqs_list:
        try:
            handle = Entrez.efetch(db="nucleotide", id="%s" %(item), rettype="gb", retmode="text")
            data = handle.read()
            handle.close()
            out_handle = open("%s%s.gb" %(folder_gb , item), "w")
            out_handle.write(data)
            out_handle.close()
            to_db(item)
        except urllib.error.HTTPError:
            print(item)
        continue

if filecsv is None:
    to_db(filegb)
