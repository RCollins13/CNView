# CNView
Visualization and annotation of CNVs from population-scale whole-genome sequencing data.

**Contact:** Ryan Collins (rcollins@chgr.mgh.harvard.edu)

All code copyright (c) 2016 Ryan Collins and is distributed under terms of the MIT license.  

---  
![Example CNView plot of a de novo 51kb duplication of SAE1 in an autism proband](/ExamplePlotA.CNView.jpg?raw=true "Example CNView plot of a de novo 51kb duplication of SAE1 in an autism proband")
---  
## Table of Contents  
#### Script documentation  
- [bidirectionalEnrichment.sh](https://github.com/RCollins13/CNView#cnviewr)  
  
#### Example usage  
- [Example #1](https://github.com/RCollins13/CNView#example1)  
---  

##CNView.R  
Performs joint normalization of binned coverage values across a batch of WGS libraries and facilitates visualization. Also interfaces with UCSC Genome Browser to underlay several annotation tracks.  

```
CNView.R [options] chr start end samples.list covmatrix.bed outfile

Options:
	-c CHARACTER, --compression=CHARACTER
		compression scalar for rebinning, if desired [default 'optimize']

	-i CHARACTER, --highlight=CHARACTER
		tab-delimited list of coordinate pairs for intervals to highlight and color as third column; NULL disables highlighting [default NA]

	-w INTEGER, --window=INTEGER
		distance to append to both sides of input interval for viewing [default 0]

	--ymin=INTEGER
		minimum value for y axis [default NULL]

	--ymax=INTEGER
		maximum value for y axis [default NULL]

	-n INTEGER, --normDist=INTEGER
		distance outside region to use for normalization (both sides) [default 5000000]

	-l, --nolegend
		disable legend on plot [default TRUE]

	-h, --help
		Show this help message and exit
``` 
**Usage Notes:**  
1.  TBD
--- 
Example A  
```
Rscript CNView.R --title="Example Plot A: Canonical Deletion, Single Sample" 2 178714141 178760307 SFARI_d12529p1 /Users/collins/scratch/SFARI_DeepSeq.nucleotide.1kb.cov_matrix.bed /Users/collins/Desktop/RCollins/Talkowski_Local/code/CNView/ExamplePlots/CNView.ExamplePlotA.pdf
```
---  
Example B

