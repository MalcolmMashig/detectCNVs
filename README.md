CBMFW4761 Final Project

Columbia University, Spring 2022

Malcolm Mashig (mjm2396) and Michael Fagan (mef2224)

### File/Folder Descriptions:

CNV-Sim: the CNV-Sim toolset for simulating CNVs based on the target regions (BED) and reference (FASTA) required input files. Some adaptations were made to ensure the simulator runs under python3. This folder also includes the minimap2 sequence alignment toolset, which is utilized by the simulator.

CNVbenchmarkeR: relevant files and scripts from https://github.com/TranslationalBioinformaticsIGTP/CNVbenchmarkeR which we explored, but did not use, as a tool for evaluation

computeAcc.R: the function computeAcc takes in the dataset of detected CNVs (result of the CNV detection function CallCNVs) and the ground truth CNVs data, and (optionally) weighted regions (based on start position). The function computes specificity and sensitivity in aggregate and for deletions and duplications respectively, and if weighted (target) regions are provided, the accuracy is computed only for these regions (since our goal is to increase accuracy in these regions).

data: (from ExomeDepth) R package data. ExomeCount.RData is a read depth dataset for four samples. This was utilized and resampled for evaluation purposes.

DESCRIPTION: (from ExomeDepth) R package metadata

gatk: the GATK sequencing toolset that is utilized by the SECNVs simulator

ground-truth.csv: the true CNVs for sample read depth counts provided by ExomeDepth (ExomeCount) which were used for a heuristic evaluation of accounting for characteristic regions

heuristic-eval.R: this script (takes a bit to run) runs a mini heuristic evaluation of potential benefits for accounting for characteristic regions with increased likelihood of CNV presence

heuristicEval-results.csv: the average accuracies achieved based on the script above

HMM_improvements: our changes to the ExomeDepth tool, thus much of this folder is identical to the home repository. Within simple_hmm is our implementation of the extra state architecture. The relevant files for CNV detection are R/class_definition.R which contains the CallCNVs function which calls the HMM function found in viterbi.hmm function in R/tools.R which utilizes the underlying Viterbi function, C_hmm (written in C++) found in src/hmm.cpp.

ICR96: the BED file (target regions) and validation data for the ICR9 exon CNV Validation data found at https://ega-archive.org/datasets/EGAD00001003335. The actual samples are not provided here because they are too large (and must be requested). We did not end up utilizing these NGS panel samples for evaluation.

inputs: a number of BED files for various target regions that we tried. chr1-exons.bed includes all exon ranges for chromosome 1, copynumber.bed includes example simulator-generated target regions. target-regions.bed are the target regions from the ICR96 dataset above, and the other BED files are subsets based on chromosome and region length.

install.R: this script can be run to install our revisions to ExomeDepth; the package will still be called ExomeDepth

installBaseline.R: this script can be run to install (or re-install) the ExomeDepth package (if testing the baseline)

Makefile: (from ExomeDepth) for building and installing the ExomeDepth package

NAMESPACE: (from ExomeDepth) namespace metadata for the ExomeDepth package

picard: a sequence manipulation toolset utilized by the SECNVs simulator

R: (from ExomeDepth) source code for R functions in ExomeDepth. The key function is CallCNVs within R/class_definition.R which takes in the result of getBamCounts found in R/countBamInGranges.R (which aggregates read depth data from an aligned sample), and returns the detected CNVs.

resampleBamCounts.R: the function resampleBamCounts takes in a number of distribution-related and CNV specification parameters for simulating new read depth data from that which is provided by the ExomeDepth package. The user can either generate data based on resampling the package data (ExomeCount) and assuming the same ground truth CNVs or generating read depth data of a similar structure with different ground truth CNVs. The former was utilized for evaluation.

samtools-1.14: the genomic sequence processing toolset that was utilized for a number of post-processing steps (converting SAM files to clean/indexed BAM files)

SECNVs: the SECNVs CNV simulator (alternative to CNV-Sim) that requires a target regions BED file and reference fasta file as input. Some adaptations were made to ensure the simulator runs under python3.

simulate.R: the function simulate is responsible for carrying out both resampleBamCounts and computeAcc (mentioned above). The function takes in a boolean variable for whether or not regions will be weighted, what those regions are (optional), what transition probability (weight) should be utilized in the HMM for those regions (optional), and lastly what transition probability should be utilized if weighting is false. This function carries out the typical ExomeDepth workflow of calling CNVs based on read depth data (which includes a number of different steps), and then essentially evaluates accuracy of the calls.

simulator-pipeline.R: This script shows the ordered command-line commands that need to be run in order to simulate CNVs (with either CNV-Sim or SECNVs), then post-process the simulation output, then summarize target regions' read depth for the simulated samples. Once read depth figures are obtained, CNVs can be detected. Please note that these commands will not work unless necessary user-based adjustments (marked with <>) are made. For example, absolute paths are sometimes required, and some inputs (i.e. fasta for chromosome one) are too large to make accessible in this repository. Additionally, all of the toolsets utilized (python3, CNV-Sim, SECNVs, Picard, GATK, Samtools, Minimap2, etc.) may require compilation.

src: (from ExomeDepth) source code for non-R functions. As mentioned above, the relevant source file is hmm.cpp which implements the Viterbi algorithm to determine the most likely sequence of hidden states (copy number characterizations) given the depth of coverage for a sample in contrast to a reference.

util: (from ExomeDepth) irrelevant

vignette: (from ExomeDepth) irrelevant

### Experiment Replication Directions:

Assuming you have cloned this directory locally...

To run the simulators (CNV-Sim or SECNVs), please reference the description of simulator-pipeline.R above. Python3 and all of the tools utilized may require compilation. To download the hg38 chromosome one fasta, please visit https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/. 

To run the heuristic evaluation, please reference the description of heuristic-eval.R above. Since this is an R script, R must be installed and a few R packages need to be installed including fitdistrplus, GenomicRanges, and tidyverse. These can be installed in R with install.packages("PACKAGE_NAME").

To install our revisions to the ExomeDepth package (overwrite the ExomeDepth namespace), please reference the description above for install.R. To return to the baseline package installation, reference installBaseline.R.
