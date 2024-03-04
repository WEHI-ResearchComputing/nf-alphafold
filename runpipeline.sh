#!/bin/bash
module purge
module load nextflow/23.10.0 
rm -rf /vast/scratch/users/$USER/results
nextflow run main.nf -profile milton 
