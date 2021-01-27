[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4458917.svg)](https://doi.org/10.5281/zenodo.4458917)


# LTDE



Repository for code associated with the preprint: "Microbial population dynamics and evolutionary outcomes under extreme energy-limitation"

https://doi.org/10.1101/2021.01.25.428163



## Dependencies
All code was written in Python 3.8 or R v3.5.0.

The following Python packages are required: numpy, pandas, matplotlib, scipy, statsmodels, and Biopython.

The following R packages are required: ggplot2, latex2exp, ggpubr, vegan, ggrepel, phylolm, MuMIn, viridis, lme4, metaMS, XML, xcms, scales, ggthemes, seqinr, ggtree, treeio, gridExtra, pmc, ape, bbmle, devtools, plotrix, pracma


Set up the repo under a folder named Github: `~/GitHub/LTDE/`

## Getting the data

Assembled genomes are available on NCBI and accession numbers are available in `genomes_info.txt`. Raw reads are available on SRA under BioProject PRJNA561216. All other data is on Zenodo data repository doi:10.5281/zenodo.4458917.


## Running the analyses

Once you have downloaded the data, just run `run_everything.sh`. This script contains commands for all data processing and analyses, including mutation calling.
