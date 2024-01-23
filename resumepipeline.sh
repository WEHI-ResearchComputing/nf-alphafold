#!/bin/bash
module purge
#Uncomment line below if using conda environments
#module load miniconda3/latest 
#Uncomment line below if using containers
#module load apptainer/1.0.0 
module load nextflow/23.10.0 
nextflow run main.nf -profile milton -resume