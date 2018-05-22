# Linux commands to analyze T5 viruses
Download the [blast+ software](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and install it
then it will be used to align faa files against the database
```
$ export PATH=$PATH:~/ncbi-blast-2.7.1+/bin ; mkdir ~/Blast/T5st0_faa ~/Blast/T5st0_blast_out ; cd ~/Blast/
$ python faa_generator.py -g "T5st0.gb" -gb "/home/ramirez/Blast/" -fa "/home/ramirez/Blast/T5st0_faa/"
$ cd ~/Blast/T5st0_faa ; for f in ~/Blast/T5st0_faa/* ; do tblastn -query $f -evalue 1e-3 -remote -db nr -out ~/Blast/T5st0_blast_out/$(basename $f).csv -outfmt '10 sseqid qseqid evalue bitscore qcovs pident' ; echo $(basename $f) ; done
```

to analyze only the FST region, copy the files from T5p001 to T5p0016
```
$ mkdir ~/Blast/FST/Blast 
$ cd ~/Blast/FST/Blast ; for f in ~/Blast/FST/Blast/* ; do cat $f >> ~/Blast/FST/simseqs-FST.csv ; echo $(basename $f) ; done
```

Download genomes -> store protein sequences
```
$ mkdir ~/Blast/FST/T5viruses/ ~/Blast/FST/T5viruses_faa ; cd ~/Blast ; python faa_generator.py -s "/home/ramirez/Blast/FST/simseqs-FST.csv" -gb "/home/ramirez/Blast/FST/T5viruses/" -fa "/home/ramirez/Blast/FST/T5viruses_faa/" && cd ~/Blast/FST/T5viruses_faa/ ; find . -type f -empty -delete
```
Result: when used the T5st0 FST region as seed: 122 viruses containing at least one T5st0 homologous gene.
All genomes not belonging to phages, or only contigs were not considered. 

## Cluster homologous gene families
Install the [get_homologues script](https://github.com/eead-csic-compbio/get_homologues/) to find orthologous genes 
and [IQtree](http://www.iqtree.org/) to build the tree.

Specify the input directory (`-d`)
- Clustering algorithm: COGtriangles (`-G`)
- Reference genome (`-r`): T5st0-Escherichia.gb
- Report all the clusters (`-t 0`) and the average % sequence identity values among pairs of genomes (`-A`)
- Request a pan- and core-genome analysis of the input genomes (`-c`)
```
$ cd ~/get_homologues ; ./get_homologues.pl -d ~/Blast/FST/T5viruses -G -r "T5st0-Escherichia.gb" -t 0 -c -A
```

## Build a pangenome matrix
Specify the input (`-d`) and output (`-o`) directories
- Produce intersection pangenome matrices (`-m`)
```
$ ./compare_clusters.pl -d ~/get_homologues/T5viruses_homologues/T5st0-Escherichia_f0_0taxa_algCOG_e0_ -o ~/Blast/FST/T5viruses_Intersection -t 0 -m
```

## Build a ML Tree
- 1000 SH-aLRT test (`-alrt`) 
- 1000 bootstraps (`-bb`)
```
$ iqtree -s ~/Blast/FST/T5viruses_Intersection/pangenome_matrix_t0.fasta -alrt 1000 -bb 1000 
```
Finally, you can use [iTOL](https://itol.embl.de/) to visualize the tree
