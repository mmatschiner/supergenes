# Code repository for analyses conducted in Matschiner et al. (2021)

This code repository is associated with our manuscript entitled "Origin and fate of supergenes in Atlantic cod"<!--, which has been posted to BioRxiv ([]())-->.

More information on the individual analyses can be found in the `README` files provided in subdirectories (e.g. [`inversion_detection/README`](inversion_detection/README)).

Analyses were conducted in the following order:

1. [`mito_nuclear_comparison`](mito_nuclear_comparison)
2. [`inversion_detection`](inversion_detection)
3. [`three_taxon_alignments`](three_taxon_alignments)
4. [`gadidae_time_calibration`](gadidae_time_calibration)
5. [`gadidae_phylogenomics`](gadidae_phylogenomics)
6. [`cod_phylogenomics`](cod_phylogenomics)
7. [`demography`](demography)

Note that the files named `run_all.sh` that can be found in each `src` directory (e.g. [`inversion_detection/src/run_all.sh`](inversion_detection/src/run_all.sh)) are not intended to be executed, but they specify the order in which all other scripts should be executed.