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

include {  CountUniqueSequences } from './modules/alphafold.nf'
include {  ALPHAFOLD_Feature as Monomer_Feature } from './modules/alphafold.nf'
include {  ALPHAFOLD_Feature as Multimer_Feature } from './modules/alphafold.nf'

include {  ALPHAFOLD_Inference as Monomer_Inference } from './modules/alphafold.nf'
include {  ALPHAFOLD_Inference as Multimer_Inference } from './modules/alphafold.nf'



workflow {

    def query_ch = Channel.fromPath(params.inputdir+"/*.fasta",checkIfExists:true)
                          .ifEmpty {
                                    error("""
                                    No samples could be found! Please check whether your input directory
                                    is correct, and that your sample name matches *.fasta.
                                    """)
                          }
    
    Channel.from(params.model_indices.split(',').toList())
           .set { model_indicies_ch }
    
    count_ch=CountUniqueSequences(query_ch)

    count_ch.map{ name,file,count ->
            return tuple(name,file,count.splitText(limit:1).first().trim().toInteger())
            }
            .branch { name,file,count ->
                monomer  : count == 1 
                    return tuple(name,file,"monomer_ptm")
                multimer : count > 1 
                    return tuple(name,file,"multimer")
            }
            .set { inference_ch }
    
    Mutli_feature_ch=Multimer_Feature(inference_ch.multimer)
    
    Mono_feature_ch=Monomer_Feature(inference_ch.monomer)

    Monomer_Inference(Mono_feature_ch.combine(model_indicies_ch))
    Multimer_Inference(Mutli_feature_ch.combine(model_indicies_ch))
    
}