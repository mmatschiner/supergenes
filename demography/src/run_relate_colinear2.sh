#!/bin/sh

# JMIB [jmi.barth@gmail.com], 2020-06-10
# Wrapper for demography estimation using Relate.

## Run Relate again to infer the demography.
sbatch run_relate_prep_demography.slurm
