process CountUniqueSequences {
    label 'Count'

    tag "${fasta}"

    input:
    path(fasta)


    output:
    tuple val(fasta.baseName),path(fasta),path("${fasta.baseName}_count.txt")

    script:
    """
    grep ">" ${fasta} | sort | uniq | wc -l > ${fasta.baseName}_count.txt
    """
}


process ALPHAFOLD_Feature{
    
    label 'Alphafold2_cpu'

    tag "${fasta}"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/", mode: 'move', pattern: "${fasta}/*.pkl"
    publishDir "${params.outdir}", mode: 'move', pattern: "${fasta}/msas/*"

    input:
    tuple val(fasta),path(fasta_file),val(preset)

    output:
    tuple   val(fasta),
            path(fasta_file),val(preset), emit: meta
    tuple   val(fasta),
            path("${fasta}"), emit:feature

    script:
    
    """
    alphafold -f -o ./  -m $preset \
            -i $params.num_predictions \
            -t $params.max_template_date $fasta_file
    """
}

process ALPHAFOLD_Inference 
{
    label 'Alphafold2'

    tag "${fasta}"

    publishDir "${params.outdir}/", mode: 'move', pattern: "${fasta}/*.pdb"
    publishDir "${params.outdir}/", mode: 'move', pattern: "${fasta}/*.json"
    publishDir "${params.outdir}/", mode: 'move', pattern: "${fasta}/*.pkl"
    publishDir "${params.outdir}", mode: 'move', pattern: "${fasta}/plots/*.pdf"

    
    input:
        tuple   val(fasta),
                path(fasta_file),val(preset),val(model_index),
                path(msa)

    output:
        tuple val(fasta),val(preset),path("${fasta}/*model*.pdb"),path("${fasta}/*model*.pkl"), emit:pdb
        tuple val(fasta),val(preset),path(fasta_file), emit:pdb_meta
        path("${fasta}/*.pdb")
        path("${fasta}/*.json")
        path("${fasta}/plots/*.pdf")

    script:
   
    """
    #!/bin/bash    
    alphafold  -o ./ -t $params.max_template_date \
               -u \
               -g  true \
               -m $preset  \
               -n $model_index \
               -i $params.num_predictions \
               -r none \
               $fasta_file
    """
}

process ALPHAFOLD_Relax_Only{
    label 'Alphafold2'
    tag "${fasta}"
    stageInMode 'copy'
    
    publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/ranked*.pdb"
    publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/relaxed*.pdb"
    
    input:
    tuple val(fasta),val(preset),path(unrelaxed),path(unrelaxed_results),path(fasta_file),path(msa)

    output:
    tuple val(fasta),val(preset),path("${fasta}/unrelaxed*.pdb"),path(unrelaxed_results),path("${fasta}/features.pkl"), emit:pdb
    tuple val(fasta), path("${fasta}/relaxed*.pdb") , emit:pdb_relaxed

    script:
    """
    
    cp *model*.pdb ${fasta}/
    cp *model*.pkl ${fasta}/
    
    alphafold  -o ./ -t $params.max_template_date \
               -g  true \
               -m $preset  \
               -n 0,1,2,3,4 \
               -j \
               -i ${params.num_predictions} \
               -r ${params.model_to_relax} \
               ${fasta}.fasta
    """
}