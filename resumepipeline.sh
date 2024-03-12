#!/bin/bash
module purge
module load nextflow/23.10.0 
nextflow run main.nf -profile milton -resume -with-dag flowchart.html -preview