# CNView
Visualization and annotation of CNVs from population-scale whole-genome sequencing data.

**Contact:** Ryan Collins (rcollins@chgr.mgh.harvard.edu)

All code copyright (c) 2016 Ryan Collins and is distributed under terms of the MIT license.  

---  
![CNView: Example CNV Visualization](/ExamplePlots/CNView.BannerExample.jpg?raw=true "CNView: Example CNV Visualization")  

---  
## Table of Contents  
#### General Information  
- [CNView Summary](https://github.com/RCollins13/CNView#cnview-summary)
- [Accessing Reference Libraries](https://github.com/RCollins13/CNView#accessing-reference-libraries)  
- [Citing CNView](https://github.com/RCollins13/CNView#citing-cnview)

#### Script documentation  
- [CNView.R](https://github.com/RCollins13/CNView#cnviewr)  
  
#### Example usage  
- [Example A: Canonical Deletion, Single Sample](https://github.com/RCollins13/CNView#example-a)  
- [Example B: Multiple Samples on the Same x-Axis](https://github.com/RCollins13/CNView#example-b)  
- [Example C: Multiple Highlighted Intervals](https://github.com/RCollins13/CNView#example-c)  
- [Example D: Disabled UCSC Annotations](https://github.com/RCollins13/CNView#example-d)  

---  
##General Information
###CNView Summary  
CNView is a low-profile visualization tool for read depth in batches of next-generation sequencing libraries and, more specifically, for visually inspecting sites of [copy-number variation (CNV)](https://en.wikipedia.org/wiki/Copy-number_variation).  

Beta renditions of CNView plots have been featured in several prior publications ([Brand & Pillalamarri et al., 2014](http://www.ncbi.nlm.nih.gov/pubmed/25279985); [Brand & Collins et al., 2015](http://www.ncbi.nlm.nih.gov/pubmed/26094575)). The developers are happy to work with groups to generate customized or "publication-quality" plots, too: contact us at [rcollins@chgr.mgh.harvard.edu](mailto:rcollins@chgr.mgh.harvard.edu) to discuss options.

###Accessing Reference Libraries  
Like many depth-based CNV callers, CNView does not work on individual libraries. Instead, CNView jointly models multiple libraries simultaneously, normalizing both within and across libraries to reduce systematic coverage biases. Thus, while technically CNView will run successfully from just two libraries, generally CNVs become clearly resolvable with at least 20 total samples jointly modelled.  

If you do not have access to an appropriate cohort of normative reference samples, thousands of standard Illumina WGS libraries are available from the [1,000 Genomes Project](http://www.1000genomes.org/) and many hundreds of other, more esoteric library types are available through other public repositories like the [European Nucleotide Archive (ENA)](http://www.ebi.ac.uk/ena) or the [Sequence Read Archive (SRA)](http://www.ncbi.nlm.nih.gov/sra).  

###Citing CNView  
At present, there is no citation specifically for CNView, although we hope to publish the method in the future. If you use CNView, please cite either [Brand & Pillalamarri et al., 2014](http://www.ncbi.nlm.nih.gov/pubmed/25279985) or [Brand & Collins et al., 2015](http://www.ncbi.nlm.nih.gov/pubmed/26094575) until there is a dedicated CNView citation, as the CNView plots themselves are described in both.  

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
1. Several base-level options not wrapped by Rscript (such as ```returnData``` or ```plot```). To access those features, import the R function directly and run from the R command prompt.  

---  
##Example A 
####Canonical Deletion Plotted in a Single Sample  
The basic use-case for CNView is to visualize a known CNV locus 
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
