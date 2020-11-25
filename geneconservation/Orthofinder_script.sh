#!/bin/bash
##------------------------------------------------------------------------------
## Orthofinder_script
## Author: Luis Ramirez
## Affiliation: I2BC
## Aim: A workflow for essential phage genes comparison.
## Date: 2019111401
##------------------------------------------------------------------------------

while getopts a: name
do
        case $name in
a) p1="$OPTARG";; 
        esac
done

echo "=============================================================="
echo "Installing requirements"
echo "=============================================================="

if [ -d ~/OrthoFinder ]; then
 echo "OrthoFinder installed"
else
 wget https://github.com/davidemms/OrthoFinder/releases/download/2.3.3/OrthoFinder-2.4.0.tar.gz -P ../
 mkdir ~/OrthoFinder ; tar xzf ../OrthoFinder-2.4.0.tar.gz -C ~/OrthoFinder --strip-components 1
fi

conda --version &> /dev/null
if [ $? -eq 0 ]; then
 echo "conda installed"
else
  export PATH=~/anaconda3/bin:$PATH
  conda --version &> /dev/null
  if [ $? -eq 0 ]; then
   echo "conda installed"
  else
   wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
   bash Anaconda3-2020.07-Linux-x86_64.sh; rm Anaconda3-2020.07-Linux-x86_64.sh
   export PATH=~/anaconda3/bin:$PATH
   conda --version
  fi
fi

echo "Installing requirements: fastme, diamond, mcl, blast, and entrez-direct"
conda install -y -c bioconda fastme
conda install -y -c bioconda diamond
conda install -y -c bioconda mcl
conda install -y -c bioconda blast
conda install -y -c bioconda entrez-direct
echo "=============================================================="
echo "Creation of tree structure"
echo "=============================================================="

mkdir Project
mkdir Project/genomes
mkdir Project/outcome

echo "=============================================================="
echo "Dowload Genomes"
echo "=============================================================="

cd Project
if [ -e ./outcome/downloaded_genomes.txt ]; then
 echo "faa files downloaded"
else
 IFS=$'\n'	# make newlines the only separator
 for i in $(tail -n +2 ../taxidlist.txt)	# iterate over all lines in taxidlist.txt
  do
   if [ -e ./genomes/$(echo "$i").faa ]	#verify if the file has not been already downloaded
    then
     echo "$(echo "$i") already downloaded"
    else
     echo "txid$i"
     esearch -db nuccore -query $(echo "txid$i") | elink -target protein | efetch -format fasta | awk '/^>/ {$0=$1} 1' > ./genomes/$(echo "$i.faa")
   fi
 done
 cd genomes ; find . -size 0 -delete
 for i in *.faa; do mv "$i" "`echo $i | sed 's/\r//'`"; done
 ls -l | awk '{print $5 "\t" $9}' > ../outcome/downloaded_genomes.txt ; cd ..
 cp ./outcome/downloaded_genomes.txt $p1
fi

echo "=============================================================="
echo "Orthogroups construction"
echo "=============================================================="

if [ -e ./outcome/Orthogroups.tsv ]; then
 echo "Orthofinder results already available"
else
 ~/OrthoFinder/orthofinder -f ./genomes -S diamond
 lastfol=$(ls -td -- ./genomes/*/* | head -n 1)
 cp "$lastfol"/Orthogroups/Orthogroups.tsv $p1
 cp "$lastfol"/Orthogroups/Orthogroups.GeneCount.tsv $p1
fi

echo "=============================================================="
echo "Done"
echo "=============================================================="

exit

