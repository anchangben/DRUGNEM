# DRUG-NEM
An optimized individualized drug combination strategy for intratumoral heterogeneity using single-cell data


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
