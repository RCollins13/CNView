# CNView
Visualization and annotation of CNVs from population-scale whole-genome sequencing data.

**Contact:** Ryan Collins (rcollins@chgr.mgh.harvard.edu)

All code copyright (c) 2016 Ryan Collins and is distributed under terms of the MIT license.  

---  
![Canonical Deletion Plotted in a Single Sample](/ExamplePlots/CNView.ExamplePlotA.jpg?raw=true "Canonical Deletion Plotted in a Single Sample")  
---  
## Table of Contents  
#### Script documentation  
- [CNView.R](https://github.com/RCollins13/CNView#cnviewr)  
  
#### Example usage  
- [Example A: Canonical Deletion, Single Sample](https://github.com/RCollins13/CNView#examplea)  
- [Example B: Multiple Samples on the Same x-Axis](https://github.com/RCollins13/CNView#exampleb)  
- [Example B: Multiple Highlighted Intervals](https://github.com/RCollins13/CNView#examplec)  
---  

##CNView.R  
Performs joint normalization of binned coverage values across a batch of WGS libraries and facilitates visualization. Also interfaces with UCSC Genome Browser to underlay several annotation tracks.  

```
Usage: ./CNView.R [options] chr start end samples.list covmatrix.bed outfile

Options:
  -c INTEGER, --compression=INTEGER
    compression scalar for rebinning, if desired [default 'NULL']

  -i CHARACTER, --highlight=CHARACTER
    tab-delimited list of coordinate pairs for intervals to highlight and color as third column; NULL disables highlighting [default NA]

  -w INTEGER, --window=INTEGER
    distance to append to both sides of input interval for viewing [default 60% of plot interval]

  --ymin=INTEGER
    minimum value for y axis [default NULL]

  --ymax=INTEGER
    maximum value for y axis [default NULL]

  -n INTEGER, --normDist=INTEGER
    distance outside region to use for normalization (both sides) [default 5000000]

  -t CHARACTER, --title=CHARACTER
    custom title for plot [default NULL]

  -q, --quiet
    disable verbose output [default FALSE]

  -l, --nolegend
    disable legend on plot [default TRUE]

  -h, --help
    Show this help message and exit
``` 
**Usage Notes:**  
1.  Several base-level options not wrapped by Rscript (such as ```returnData``` or ```plot```). To access those features, import the R function directly and run from the R command prompt.  
---  
##Example A 
####Canonical Deletion Plotted in a Single Sample
![Canonical Deletion Plotted in a Single Sample](/ExamplePlots/CNView.ExamplePlotA.pdf?raw=true "Canonical Deletion Plotted in a Single Sample")  
```
bash$ ./CNView.R 2 178714141 178760307 SFARI_d12529p1 \
                 ~/cov_matrix.bed \
                 ./ExamplePlots/CNView.ExamplePlotA.pdf \
                 --title "Example Plot A: Canonical Deletion, Single Sample"
```
---  
Example B

