# CNView
Visualization and annotation of CNVs from population-scale whole-genome sequencing data.

**Contact:** Ryan Collins (rcollins@chgr.mgh.harvard.edu)

All code copyright (c) 2016 Ryan Collins and is distributed under terms of the MIT license.  

---  
![CNView: Example CNV Visualization](/ExamplePlots/CNView.BannerExample.jpg?raw=true "CNView: Example CNV Visualization")  
---  
## Table of Contents  
#### Script documentation  
- [CNView.R](https://github.com/RCollins13/CNView#cnviewr)  
  
#### Example usage  
- [Example A: Canonical Deletion, Single Sample](https://github.com/RCollins13/CNView#example-a)  
- [Example B: Multiple Samples on the Same x-Axis](https://github.com/RCollins13/CNView#example-b)  
- [Example C: Multiple Highlighted Intervals](https://github.com/RCollins13/CNView#example-c)  
- [Example D: Disabled UCSC Annotations](https://github.com/RCollins13/CNView#example-d)  
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
- Several base-level options not wrapped by Rscript (such as ```returnData``` or ```plot```). To access those features, import the R function directly and run from the R command prompt.  

---  
##Example A 
####Canonical Deletion Plotted in a Single Sample
![Canonical Deletion Plotted in a Single Sample](/ExamplePlots/CNView.ExamplePlotA.jpg?raw=true "Canonical Deletion Plotted in a Single Sample")  
```
bash$ ./CNView.R 2 178714141 178760307 SFARI_d12529p1 \
                 ~/cov_matrix.bed \
                 ./ExamplePlots/CNView.ExamplePlotA.pdf \
                 --title "Example Plot A: Canonical Deletion, Single Sample"
```
---  
##Example B  
####Multiple Samples on the Same x-Axis  
![Multiple Samples on the Same x-Axis](/ExamplePlots/CNView.ExamplePlotB.jpg?raw=true "Multiple Samples on the Same x-Axis")  
```
bash$ ./CNView.R 19 47636953 47687179 \
                 ~/sampleIDs.txt
                 ~/cov_matrix.bed \
                 ./ExamplePlots/CNView.ExamplePlotB.pdf \
                 --title "Example Plot B: Multiple Samples on the Same x-Axis" \
                 -w 50000 \
                 --ymin -6 \
                 --ymax 10 
```
---  
##Example C  
####Multiple Highlighted Intervals  
![Multiple Highlighted Intervals](/ExamplePlots/CNView.ExamplePlotC.jpg?raw=true "Multiple Highlighted Intervals")  
```
bash$ ./CNView.R X 36605888 37148341 \
                 ~/sampleIDs_2.txt \
                 ~/cov_matrix.bed \
                 ./ExamplePlots/CNView.ExamplePlotC.pdf \
                 --highlight ~/highlights.txt \
                 --title "Example Plot C: Multiple Highlighted Intervals" \
                 -c 10 \
                 -w 300000 \
                 --ymin -10 \
                 --ymax 3 
```
---  
##Example D  
####Disabled UCSC Annotations  
![Disabled UCSC Annotations](/ExamplePlots/CNView.ExamplePlotD.jpg?raw=true "Disabled UCSC Annotations")  
```
bash$ ./CNView.R 2 178714141 178760307 SFARI_d12529p1 \
                 ~/cov_matrix.bed \
                 ./ExamplePlots/CNView.ExamplePlotD.pdf \
                 --title "Example Plot D: Disabled UCSC Annotations" \
                 --noUCSC
```
