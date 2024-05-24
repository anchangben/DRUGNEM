# DRUG-NEM An optimized individualized drug combination strategy for intratumoral heterogeneity using single-cell data



Cancers are comprised of populations of cells with distinct molecular and phenotypic features, a
phenomenon termed intratumor heterogeneity. Cells in general within a tumor have biomarkers that
can be exposed to drug(s) and the response measured using high-throughput single-cell technologies like
CyTOF. Some of the markers denoted here as lineage biomarkers can be used to characterize cell types
eg. T or B cells including malignant cell types and the others used to de ne the functional states of a cell
eg. apoptotic or intracellular signaling states producing a complex viewpoint of intratumor heterogeneity.
Intratumor heterogeneity poses an important challenge to cancer therapy [1]. Drug resistance mechanisms
from targeted single drugs have been identified due to intratumor heterogeneity [2], motivating the need for
combination therapy. The increasing number of potential FDA approved targeted drugs and corresponding
resistance mechanisms identified due to intratumor heterogeneity warrants a systematic and optimal
combination therapy that accounts for intratumor heterogeneity.
We introduce a model framework called DRUG-Nested Effects Models (DRUG-NEM) for analyzing
heterogeneous drug response data across a population of single cells to optimize and prioritize(rank)
combination therapy particularly for single drugs and in general any perturbation response data. This
algorithm combines (1) Population identi cation, (2) Estimating of desired e ects (3) Nested e ects
modeling and (4) Score-Rank analysis to identify and rank multi-target drug combinations from single
drug effects measured at the level of single cells. The main function drugnemmain takes as input either
(1) a folder containing all drug response FCS  les or (2) A data matrix of pooled drug responses for
all proteins(columns) across all cells(rows) or (3) a folder with manually gated cells for drug responses
from mass cytometry comprising of n drugs and m target (marker) responses on t cells taken from single
or several patient samples and returns several pdf  files including both the drug target nested network
and the sorted list of all possible drug combinations and their scores with the best drug regimen at the
top derived from pre-specified desired e ects such as down regulation effects or effects associated with a
desired response of a desired marker e.g death marker or from a model with no effects. The algorithm
also varies the desired effects, computes the DRUG-NEM score and chooses the model with the desired
effects that produce the maximal score.
The specific output of the algorithm is a list of 33 items. Item 2 corresponds to the DRUG-NEM
ranking using "Any effect or no prior", item 10 corresponds to the DRUG-NEM ranking using prior "Up
regulation effects", item 12 correspond to the DRUG-NEM ranking using prior "Down regulation effects",
item 14 corresponds to the DRUG-NEM ranking using prior "No effects" and item 33 corresponds the
best ranking after comparing priors "Up effects" and "Down effects". The algorithm also generates several
pdf(Tiff)  files corresponding to various diagnostic plots including plots of initial and  final optimized
CCAST decision trees for all cell lineage sub groups, Heatmaps of fold changes for all functional or
intracellular markers across all sub populations and DRUG-NEM network and ranking list of all drug
combinations with corresponding score.

# Introduction
drugnem uses data derived from mass cytometry (MCM) (CyTOF), which is a recently developed high
throughput MCM technology for labeling single cells with metal-chelated antibodies that greatly reduce the
auto fluorescence effect of individual protein signals, thereby circumventing the need for robust compensation
after data generation. MCM has been adapted to generate high throughput drug screening single-cell
signaling data using mass-tag cellular barcoding (MCB). This provides the input data for drugnem
although other molecular datasets with a similar input data structure can be adapted as well.

# Installation
drugnem relies on R (>= 3.4.0) and other R libraries: 
FlowCore, party,nem, RColorBrewer, limma, RBGL,graph,Rgraphviz,statmod,
MASS, PerformanceAnalytics, fastcluster, cluster, ggplot2, car,stats4, mvtnorm, grid, modeltools, plottrix,
xts, sandwich, and zoo. Note nem is no more maintained in Bioconductor. You will need to manually install it from the included source file. 
The current drygnem implementation runs on Mac OS X 10.7 and above operating sys-
tems. R can be installed in these systems by simply downloading the most recent R-x.y.z.tar.gz  file from
http://www.r-project.org and following the system speci c instructions.To install the package using the
devtools package in R simply open an R console or click the Rstudio drugnem. Rproj object in the package
folder which automatically opens the drugnem libraries in a new stand alone Rstudio window.

# Installing and Loading the Library using R studio

We start by installing the drugnem library.
>library("devtools")
>library(roxygen2)
>setwd("path/to/drugnem folder")
install("drugnem")
We start by loading the drugnem library. Next we create a working directory (drugnem output path) to
save all drugnem output  gures and data.
>library(drugnem)
>dir.create(path/to/drugnem output)
A complete installation of drugnem requires that all dependencies are loaded.

# Installing and Loading the Library using tar or zip file
install.packages("path/to/drugnem tar or zip file", repos = NULL, type = "source")
> library(drugnem)

# Running drugnem

drugnem" makes use of four different modeling processes: (1) Population identification, (2) Estimating of
desired effects (3) Nested effects modeling and (4) Score-Rank analysis to optimize and prioritize possible
combination therapeutic strategies from experimentally derived single-drug responses on single cells using
CyTOF. The flowchart of algorithm is shown in Figure2 in the main text. CyTOF produces mass
cytometry (MCM) data typically stored in a flow cytometry standard (FCS)  file as a data frame with rows
representing the cells or events and the columns corresponding to the markers of interest. Each FCS  file or
data matrix would ideally contain lineage markers used to define baseline subpopulations and intracellular
markers including death e ector proteins such as caspases or CPARP which allow us to compare both
survival and intracellular signaling pre and post drug treatment. Briefly stated, in the absence of manually
gated subpopulations, we pool all population of single-cell expression data into a single data matrix for all
drugs. CCAST is next used to identify and match homogeneous population of cells using the lineage
markers across all treatment conditions to de ne subpopulations. A set of apoptotic marker(s) is selected
from data and then CCAST is applied again using these marker(s) to separate survival from death cells
for each subpopulation of cells identi ed from the previous step. We then pool the mixture of labeled
survival cells derived into a new data matrix and then estimate the drug e ects for each subpopulation.
With respect to the baseline treatment (C), we estimate the probability of a protein to be di erentially
expressed under each drug using linear models for all cell types. We next build an integrated drug-nested
network using nested effects models from a drug desired effects matrix using the weighttype parameter
defining the type of desired phenotype information to use for optimizing treatment effects. Values include:
(i) 'T-stat': Regulation using T-stat of intracellular signaling e ects (ii) 'FC'(default): Regulation based
on Foldchange of intracellular signaling e ects, and (iii) 'deathmarker': effects associated with regulation
of a death effector protein or any other desired marker. (See main text for description). Finally, we score
and rank all drug combinations using posterior weights of the target positions on the derived drug-target
network with the best regimens on the top. Drug combinations are scored both under the independence
and nested assumption of drug effects.

In the drugnem manual we provide key technical details and implementation steps required for running the
algorithm and further demonstrate the usage of the main function of drugnem called "drugnemmain" on
various data input types specifically, on CyTOF drug response data derived from an cervical cancer cell line
(HeLa) and Acute Lymphoblastic Leukemia (ALL) cells showing the expected results. 


# Publications

Anchang, B , et al. (2018) DRUG-NEM: Optimizing drug combinations using single-cell perturbation
response to account for intratumoral heterogeneity, Proc Natl Acad Sci U S A. 2018 May 1;115(18):E4294-E4303. PMID: 29654148; PMCID: PMC5939057.
Anchang, B, Do, M, Zhao, X, Plevritis, S. (2014) CCAST: a model-based gating strategy to isolate
homogeneous subpopulations in a heterogeneous population of single cells. PLoS Comput Biol. 10(7): p.
e1003664 .

# QUESTIONS
Any questions can be addressed to benedict.anchang@nih.gov
