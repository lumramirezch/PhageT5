# Sanger_graph
Install the sangerseq library 
```
source("https://bioconductor.org/biocLite.R")
biocLite("sangerseqR", dependencies = TRUE)
```
if it lacks a library type: install.packages("package", repos = c("http://rstudio.org/_packages","http://cran.rstudio.com"))

Import the tools
```
library(sangerseqR)
library(ggplot2)
library(Biostrings)
start <- 103 #start of the probe sequence in the seq results
seqspan <- 20 #length of the probe sequence
```

See the chromatogram and adjust `start` and `seqspan` variables if needed 
```
seqscf <- read.scf("~/../file.scf")
screen <- seqscf@sample_points[start:(start+seqspan),]
chromatogram(sangerseq(seqscf), trim5 = start, trim3= nchar(seqscf@basecalls)-start-seqspan, width = seqspan)
```


Create a matrix of probabilities based on the raw sequencing data
```
screen <- (screen/max(screen)) + 1
colnames(screen) <- c("A", "C", "G", "T")
target <- readDNAStringSet(file = "~/../file2.fasta", format = "fasta") 
targetr <- reverseComplement(target)
prob_generator <- function(probe, target){
	a <- c()
	for (i in 1:(nchar(target)-seqspan)){
		a[i] <- 1
		for (j in  (1:seqspan)){
			a[i] <- a[i]*probe[j, substr(target, i+j, i+j)]
		}
	}
	c <- max(a)
	return(a/c)
}
prob_target = prob_generator(screen, target)
prob_targetr = prob_generator(screen, targetr)
```

Plot the results of the exploration with the matrix, to know where is your sequence in the genome
```
qplot(seq_along(prob_target), prob_target, geom=c("line"))
qplot(seq_along(prob_targetr), prob_targetr, geom=c("line"))
```
