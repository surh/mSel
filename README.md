# Instructions

1. Basically you one file per bacterial species which has one row per SNP
per population (i.e. indivdual) indicating the direction of the change in
the allele frequencies.

    - I create this file with `midas2bern.r` which takes a directory produced
    by midas, where each file is a species and SNPs are in MIDAS format (for
    historical reasons, inStrain is probably better). You also need

    - `midas2bern.r` also needs a mapping file (tsv) that has the columns `pt`,
    `start`, and `end`, for patient (should be population or individual),
    sample id of the first timepoint, and sample id of the end timepoint. These
    sample IDs need to correspond to the sample ids file with the SNPs in MIDAS
    format

    - . The output `output/sites.tsv`is then the input to 'bern_mix.r` 

2. You pass that file to `bern_mix.r` which produces several outputs, but
the main one is the `output/p_directional.tsv.gz` which has the per-site SNPs
estimates of the probability of the site being under directional selection.

3. The nextflow `bern_mix.nf`is just a wrapper that automates sending many jobs
on a cluster. Current version is still using DSL1 which is
now deprecated by nextflow.

4. The nextflow `midas2bern.nf` is a wrapper that starts from MIDAS-formatted
SNPs mapping files and runs the 2 steps.  Current version is still using DSL1
which is now deprecated by nextflow.

# Troubleshooting & recommendations

Parameters --vp --vq and vq have limited effect on the model results once you
have **at least 10k SNPs**, but they do influence convergence of the HMC
approach.

It is not recommended to run analysis of a species if it is present in less
than 5 individuals and/or it has less than 10K SNPs.

Defaults work well, but ocassionally we've found that increasing the values
can lead to faster convergence. If fit fails, try changing them.

# Requirements

Need to add versions

* R
    - argparser
    - tidyverse
    - rstan
* Nextflow (still using DSL1)