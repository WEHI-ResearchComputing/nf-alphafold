process CountUniqueSequences {
    label 'Alphafold2'

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

    publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/*.pkl"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${fasta}/msas/*"

    input:
    tuple val(fasta),path(fasta_file),val(preset)

    output:
    tuple val(fasta),path(fasta_file),val(preset), emit:meta
    path("${fasta}/*.pkl")
    path("${fasta}/msas/*")

    script:
    
    """
    alphafold -f -o ./  -m $preset \
            -i $params.num_predictions \
            -t $params.max_template_date $fasta_file
    #mkdir -p ${params.outdir}/${fasta}/msas
    #cp -r ${fasta}/msas ${params.outdir}/${fasta}/msas
    #cp ${fasta}/*.pkl ${params.outdir}/${fasta}/
    """
}

process ALPHAFOLD_Inference 
{

    label 'Alphafold2'

    tag "${fasta}"

    publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/*.pdb"
    publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/*.json"
    publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/*.pkl"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${fasta}/plots/*.pdf"

    input:
        tuple val(fasta),path(fasta_file),val(preset),val(model_index)

    output:
        tuple val(fasta),val(preset),path("${fasta}/*model*.p*"), emit:pdb
        tuple val(fasta),path(fasta_file),val(preset),val(model_index), emit:pdb_meta
        path("${fasta}/*.json")
        path("${fasta}/plots/*.pdf")

    script:
   
    """
    #!/bin/bash
    mkdir -p ${fasta}/msas
    cp -r ${params.outdir}/${fasta}/msas ${fasta}
    cp ${params.outdir}/${fasta}/features.pkl ${fasta}/
    
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

    publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/ranked*.pdb"
    publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/relaxed*.pdb"
    
    input:
    tuple val(fasta),val(preset),path(unrelaxed)

    output:
    tuple val(fasta),val(preset),path("${fasta}/ranked*.pdb"), emit:pdb
    tuple val(fasta), path("${fasta}/relaxed*.pdb") , emit:pdb_relaxed

    script:
    """
    mkdir -p ${fasta}/msas
    cp -r ${params.outdir}/${fasta}/msas ${fasta}
    cp ${params.outdir}/${fasta}/features.pkl ${fasta}/
    cp *.p* ${fasta}/
    touch ${fasta}.fasta
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