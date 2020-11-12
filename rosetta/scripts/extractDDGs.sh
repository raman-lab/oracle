#!/bin/bash

#Define path to R statistical analysis package
#---------------------------------------------------------------------------------------------------------
rscript=/scratch/ameger/software/R-4.0.0/bin/Rscript
#---------------------------------------------------------------------------------------------------------

#This script compiles all of the Rosetta computed ddG values to the file "../analysis/ddg_breakdown.csv"
#---------------------------------------------------------------------------------------------------------
mkdir ../analysis/
${rscript} r/extractDDGs.R
#---------------------------------------------------------------------------------------------------------

