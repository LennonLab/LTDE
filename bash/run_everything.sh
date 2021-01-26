#!/bin/bash

# rerun mutation calling pipeline
# clean reads
#sh ~/GitHub/LTDE/bash/reads_clean.sh

# run breseq
#sh r~/GitHub/LTDE/bash/run_breseq.sh

# get breseq
#sh ~/GitHub/LTDE/bash/get_breseq.sh

# run irep
#sh ~/GitHub/LTDE/bash/run_irep.sh

# make tree
#sh ~/GitHub/LTDE/bash/raxml.sh

# run analyses for parallelism, cleaning iRep
#python ~/GitHub/LTDE/Python/clean_data.py

# identify targets of convergent evolution
#python ~/GitHub/LTDE/Python/calculate_convergence_table.py

# run demographic analyses
#R  ~/GitHub/LTDE/R/weibull.R

# examine demographic model fits
#R ~/GitHub/LTDE/R/get_mixed_model_fits.R

# plot phylogeny
#R ~/GitHub/LTDE/R/plot_tree.R

# perform phylogenetic analyses
#R ~/GitHub/LTDE/R/pmc.R

# analyze/plot metabolomics
#R ~/GitHub/LTDE/R/MSPlot.R

# plot survivial curves for all taxa
# python ~/GitHub/LTDE/Python/plot_survival_all_taxa.py

# plot longevity
# python ~/GitHub/LTDE/Python/plot_longevity.py

# plot likeligood and example survival curves
# python ~/GitHub/LTDE/Python/plot_survival_likelihood.py

# plod pN/pS
#python ~/GitHub/LTDE/Python/plot_dnds.py

# plot bacillus amino acids
#python ~/GitHub/LTDE/Python/plot_bacillus_aa.py

# plot bacillus spoiie survival curve
#python ~/GitHub/LTDE/Python/plot_spoiie_survival.py

# plot tajimas D
#python ~/GitHub/LTDE/Python/plot_tajimas_d.py

# plot proportion dead cells
#python ~/GitHub/LTDE/Python/plot_proportion_dead_cells.py

# plot irep
# python ~/GitHub/LTDE/Python/plot_irep_shape.py

# plot multiplicity synonymous vs nonsynonymous
# python ~/GitHub/LTDE/Python/plot_mult_syn_nonsyn.py

# plot lag
# python ~/GitHub/LTDE/Python/plot_lag_shape.py

# plot multiplicity vs frequency
# python ~/GitHub/LTDE/Python/plot_mult_frequency.py


# plot site frequency spectra
# python ~/GitHub/LTDE/Python/plot_afs.py

# plot diversity regressions
# python ~/GitHub/LTDE/Python/plot_dn_ds_tajimas_d.py
