# AlphaFold Nextflow Pipeline to predict structures for further DL training

This pipeline is a simplified pipeline from the main developed one. All fasta files found in the input dir are processed using only `model 1` and then relaxed. Only the final relaxed pdb is moved to the output dir.

The following parameters are pre-set in [nextflow.config](https://github.com/WEHI-ResearchComputing/nf-alphafold/blob/af-for-ai/nextflow.config) and hidden from [seqera platform](https://seqera.services.biocommons.org.au/) launchpad.
  * model_indices      = "0"         
  * num_predictions    = 1
  * model_to_relax     = "best"        
  * max_template_date  = "2024-03-01" 

This is developed for Richard Birkenshaw, Czabotar Lab
