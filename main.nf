#!/usr/bin/env nextflow
import java.text.SimpleDateFormat

include {  CountUniqueSequences } from './modules/alphafold.nf'
include {  ALPHAFOLD_Feature as Monomer_Feature } from './modules/alphafold.nf'
include {  ALPHAFOLD_Feature as Multimer_Feature } from './modules/alphafold.nf'

include {  ALPHAFOLD_Inference as Monomer_Inference } from './modules/alphafold.nf'
include {  ALPHAFOLD_Inference as Multimer_Inference } from './modules/alphafold.nf'

include {  ALPHAFOLD_Relax_Only as Monomer_Relaxation } from './modules/alphafold.nf'
include {  ALPHAFOLD_Relax_Only as Multimer_Relaxation } from './modules/alphafold.nf'


include { Generate_Report } from './modules/generate_report'



workflow {
    
    println "*****************************************************"
    println "*  Nextflow <name> pipeline                         *"
    println "*  A Nextflow wrapper pipeline                      *"
    println "*  Written by Julie Iskander,                       *"
    println "*              Research Computing Platform          *"
    println "*  research.computing@wehi.edu.au                   *"
    println "*                                                   *"
    println "*****************************************************"
    println " Required Pipeline parameters                        "
    println "-----------------------------------------------------"
    println "Input  Directory   : $params.inputdir                "
    println "Output Directory   : $params.outdir                  "
    println "Use Calculated MSA : $params.msa_calculated          " 
    println "*****************************************************"

    def query_ch = Channel.fromPath(params.inputdir+"/*.fasta",checkIfExists:true)
                          .ifEmpty {
                                    error("""
                                    No samples could be found! Please check whether your input directory
                                    is correct, and that your sample name matches *.fasta.
                                    """)
                          }
    
    Channel.from(params.model_indices.split(',').toList())
           .set { model_indicies_ch }
    
    number_of_model_output= params.model_indices.split(',').size()
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
    
    

    if (params.msa_calculated == "False"){
        
        Multimer_Feature(inference_ch.multimer)
        Monomer_Feature(inference_ch.monomer)

        Monomer_Inference(Monomer_Feature.out.meta.combine(model_indicies_ch))
        Multimer_Inference(Multimer_Feature.out.meta.combine(model_indicies_ch))

    }
    else {
        Monomer_Inference(
            inference_ch.monomer
                .combine(model_indicies_ch)
        )
       
        Multimer_Inference(
            inference_ch.multimer
                .combine(model_indicies_ch)
        )
        
    }
    Monomer_Relaxation(
        Monomer_Inference.out.pdb
            .groupTuple(by:[0,1],size:number_of_model_output)
            .map{ fasta,preset, files ->
                return tuple(fasta,preset, files.flatten())
            }
    )
    Multimer_Relaxation(
        Multimer_Inference.out.pdb
            .groupTuple(by:[0,1],size:number_of_model_output)
            .map{ fasta,preset, files ->
                return tuple(fasta,preset, files.flatten())
            }
    )
    Generate_Report(Monomer_Relaxation
                        .out.pdb
                        .mix(Multimer_Relaxation.out.pdb)
                        .map{ fasta,preset, files ->
                                return tuple(fasta,preset, files.flatten())
                        }
                        .combine(Channel.fromPath("${projectDir}/assets/proteinfold_template.html"))
    )
}
