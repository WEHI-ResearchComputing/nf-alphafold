#!/usr/bin/env nextflow
import java.text.SimpleDateFormat
println "*****************************************************"
println "*  Nextflow <name> pipeline                         *"
println "*  A Nextflow wrapper pipeline                      *"
println "*  Written by Research Computing Paltform           *"
println "*  research.computing@wehi.edu.au                   *"
println "*                                                   *"
println "*****************************************************"
println " Required Pipeline parameters                        "
println "-----------------------------------------------------"
println "Input  Directory   : $params.inputdir                   "
println "Output Directory   : $params.outdir                  " 
println "*****************************************************"

// include modules


// step_1_output_dir = "$params.output_dir"
// input_files = "${params.input_dir}/*.fasta" // match everything else


process ALPHAFOLD_Feature{
    
    label 'Alphafold2'

    tag "${fasta}"

    input:
    path(fasta)

    output:
    path(fasta)

    module 'alphafold/2.3.2'
    script:
    
    """
    alphafold -f -o $params.outdir \
             -t $params.template_date $fasta
    """
}

process ALPHAFOLD_Inference{
    queue 'gpuq'
    clusterOptions '--gres=gpu:A30:1 --nice'
 
    label 'Alphafold2'
    tag "${fasta}"

    input:
    tuple val(model_index),path(fasta)

    output:
    path("*.pdb")

    module 'alphafold/2.3.2'
    script:
    """
    today=\$(date +%F)
    alphafold  -o $params.outdir -t $params.template_date \
               -g  true \
               -m $params.model_preset  \
               -n $model_index \
               -i $params.num_predictions \
               -r $params.model_to_relax \
               $fasta
    """
}


workflow {

    //Check max_template_date Param
    if (params.max_template_date==null || params.max_template_date.isEmpty())
        params.template_date = new SimpleDateFormat("yyyy-MM-dd").format(new Date()) 
    else
        params.template_date = params.max_template_date
    
    
    
    def query_ch = Channel.fromPath(params.inputdir+"/*.fasta",checkIfExists:true)
                          .ifEmpty {
                                    error("""
                                    No samples could be found! Please check whether your input directory
                                    is correct, and that your sample name matches *.fasta.
                                    """)
                          }
    

    fasta_ch=ALPHAFOLD_Feature(query_ch)

    Channel.from(params.model_indices.split(',').toList())
           .set { model_indicies_ch }

    ALPHAFOLD_Inference(model_indicies_ch.combine(fasta_ch))
    
}