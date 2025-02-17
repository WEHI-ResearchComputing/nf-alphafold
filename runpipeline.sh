#!/bin/bash
module purge
module load nextflow/24.04.2
nextflow run main.nf -profile milton debug
