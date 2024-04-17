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
    
    label 'Alphafold2'

    tag "${fasta}"

    //publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/*.pkl"
    //publishDir "${params.outdir}", mode: 'copy', pattern: "${fasta}/msas/*"

    input:
    tuple val(fasta),path(fasta_file),val(preset)

    output:
    tuple val(fasta),path(fasta_file),val(preset),path("${fasta}/*.pkl"),path("${fasta}/msas/*"),  emit:output
    
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

process ALPHAFOLD_Inference{
    queue 'gpuq'
    clusterOptions '--gres=gpu:A30:1 --nice'
    errorStrategy 'ignore'
    label 'Alphafold2'
    tag "${fasta}"

    publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/*.pdb", saveAs: { filename -> "${fasta}" }
    //publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/*.json"
    //publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/*.pkl"
    //publishDir "${params.outdir}", mode: 'copy', pattern: "${fasta}/plots/*.pdf"

    output:
    path("${fasta}/*.pdb")
    path("${fasta}/*.json")
    path("${fasta}/*.pkl")
    path("${fasta}/plots/*.pdf")


    input:
    tuple val(fasta),path(fasta_file),val(preset),path(pkl),path(msas),val(model_index)

    script:
    """
    mkdir -p ${fasta}/msas
    mv  bfd_uniref_hits.a3m ${fasta}/msas/
    mv  mgnify_hits.sto ${fasta}/msas/
    mv  uniref90_hits.sto ${fasta}/msas/
    mv  pdb_hits.hhr ${fasta}/msas/
    mv features.pkl ${fasta}/
    alphafold  -o ./ -t $params.max_template_date \
               -g  true \
               -m $preset  \
               -n $model_index \
               -i $params.num_predictions \
               -r $params.model_to_relax \
               $fasta_file
    """
}