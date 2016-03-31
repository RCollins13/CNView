# CNView
Visualization and annotation of CNVs from population-scale whole-genome sequencing data.

**Contact:** Ryan Collins (rcollins@chgr.mgh.harvard.edu)

All code copyright (c) 2016 Ryan Collins and is distributed under terms of the MIT license.  

---  
![CNView: Example CNV Visualization](/ExamplePlots/CNView.BannerExample.jpg?raw=true "CNView: Example CNV Visualization")  

---  
## Table of Contents  
####[General Information](https://github.com/RCollins13/CNView#general-information-1)  
- [CNView Summary](https://github.com/RCollins13/CNView#cnview-summary)
- [Accessing Reference Libraries](https://github.com/RCollins13/CNView#accessing-reference-libraries)  
- [Citing CNView](https://github.com/RCollins13/CNView#citing-cnview)

####[Code documentation](https://github.com/RCollins13/CNView#code-documentation-1)  
- [CNView.R](https://github.com/RCollins13/CNView#cnviewr)  
  
####[Example Usage](https://github.com/RCollins13/CNView#example-usage-1)  
- [Getting Started](https://github.com/RCollins13/CNView#getting-started-1)
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
##Code Documentation
###CNView.R  
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
##Example Usage
###Getting Started  
**Required Input Data**  
The input data for CNView is a tab-delimited [bed-style](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) matrix of binned coverage values, which can be generated by tools like ```bedtools coverage``` ([bedtools homepage](http://bedtools.readthedocs.org/en/latest/)). Generally, 100bp-1kb sequential bins provide a reasonable tradeoff between resolution and modeling speed.

An example coverage matrix at 100bp binsize for human libraries aligned to reference genome hg19 would look something like this:  
```
Chr  Start     End       SampleA  SampleB  SampleC  ...  SampleZ
1    1         100       89       56       217      ...  141
1    100       200       98       60       230      ...  132
1    200       300       102      59       202      ...  142
...  ...       ...       ...      ...      ...      ...  ...
Y    59373400  59373500  93       68       198      ...  129
```  

CNView has been tested on libraries ranging from 1X to >300X coverage simultaneously and appears to perform relatively consistently irrespective of the ranges of coverage between individual libraries in the same batch.  


###Example A  
**Canonical Deletion Plotted in a Single Sample**  
The basic use-case for CNView is to visualize a known CNV locus.  This can be performed right out of the box by invoking ```CNView.R``` with all default parameters.  
![Canonical Deletion Plotted in a Single Sample](/ExamplePlots/CNView.ExamplePlotA.jpg?raw=true "Canonical Deletion Plotted in a Single Sample")  
**Example code to generate the above plot:**
```
bash$ ./CNView.R 2 178714141 178760307 SFARI_d12529p1 \
                 ~/cov_matrix.bed \
                 ./ExamplePlots/CNView.ExamplePlotA.pdf \
                 --title "Example Plot A: Canonical Deletion, Single Sample"
```  
This example is visualizing a 46kb deletion of two exons from *PDE11A*. The first three parameters were the coordinates of the deletion (the highlighted region), while the other parameters were as follows:  
- ```SFARI_d12529p1``` is the ID of the sample being plotted, which has to match one of the names of the columns in the coverage matrix exactly.  
- ```~/cov_matrix.bed``` is the path to the input coverage matrix, like the [example](https://github.com/RCollins13/CNView#getting-started-1) provided above.  
- ```./ExamplePlots/CNView.ExamplePlotA.pdf``` is the path to the desired output file (always will be pdf).  
- ```--title``` overrides the default title with the string in quotes.  

Running the above code should generate the above plot as a pdf and will also print the following expected text to stdout, which can be silenced with ```-q```/```--quiet```:  
```
+-------------------+
| CNView Visualizer |
|     (c) 2016      |
+-------------------+
Sample ID file 'SFARI_d12529p1' not found, assuming single sample ID provided
Attempting to connect to UCSC Genome Browser... Complete
Filtering & loading coverage matrix... Complete
Compressing coverage matrix [1,000 bp bins]...  Complete
Performing intra-sample normalization... Complete
Performing inter-sample normalization... Complete
Plotting samples to ./ExamplePlots/CNView.ExamplePlotA.pdf... Complete
Appending UCSC tracks... Complete

** FINISHED ON Thu Mar 31 14:12:43 2016 **
```
---  
###Example B  
**Multiple Samples on the Same x-Axis**  
CNView can plot up to 300 samples stacked atop the same shared x-axis. This is accomplished by passing a text file with the IDs of all samples to be plotted in place of the sample ID.  
![Multiple Samples on the Same x-Axis](/ExamplePlots/CNView.ExamplePlotB.jpg?raw=true "Multiple Samples on the Same x-Axis")  
**Example code to generate the above plot:**  
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
In this example, we're visualizing a 51kb duplication of *SAE1* that occurred *de novo* in the child (top panel, sample SFARI_d13874p1) and is present in neither the father (middle panel, SFARI_d13874fa) nor the mother (bottom panel, SFARI_d13874mo). In this example, we're also using a few more options to illustrate how to control the plot parameters:  
- ```~/sampleIDs.tx``` is the path to the text file containing all sample IDs to be plotted (max of 300 samples/plot, which is an R system default, although this can be overwritten if needed).  
- ```-w 50000``` tells CNView to pad 50kb of additional space to the plotting window on both sides of the specified interval. By default, the plotting window will be padded with [61.8%](https://en.wikipedia.org/wiki/Golden_ratio) of the size of the interval.  
- ```--ymim``` and ```--ymax``` override the default y axes, which is most useful when plotting multiple samples simultaneously to get a relative.  

And here's the expected stdout text:  
```
+-------------------+
| CNView Visualizer |
|     (c) 2016      |
+-------------------+
Sample ID file 'SFARI_d12529p1' not found, assuming single sample ID provided
Attempting to connect to UCSC Genome Browser... Complete
Filtering & loading coverage matrix... Complete
Compressing coverage matrix [1,000 bp bins]...  Complete
Performing intra-sample normalization... Complete
Performing inter-sample normalization... Complete
Plotting samples to ./ExamplePlots/CNView.ExamplePlotA.pdf... Complete
Appending UCSC tracks... Complete

** FINISHED ON Thu Mar 31 14:12:43 2016 **
```
---  
###Example C  
**Multiple Highlighted Intervals**  
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
**Disabled UCSC Annotations**  
![Disabled UCSC Annotations](/ExamplePlots/CNView.ExamplePlotD.jpg?raw=true "Disabled UCSC Annotations")  
```
bash$ ./CNView.R 2 178714141 178760307 SFARI_d12529p1 \
                 ~/cov_matrix.bed \
                 ./ExamplePlots/CNView.ExamplePlotD.pdf \
                 --title "Example Plot D: Disabled UCSC Annotations" \
                 --noUCSC
```
