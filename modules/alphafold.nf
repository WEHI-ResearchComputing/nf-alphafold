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
    
    label 'Alphafold2'

    tag "${fasta}"

    input:
    tuple val(fasta),path(fasta_file),val(preset)

    output:
    tuple val(fasta),path(fasta_file),val(preset), emit:output
    path("${fasta}/*.pkl")



    module 'alphafold/2.3.2'
    script:
    
    """
    alphafold -f -o ./  -m $preset \
            -i $params.num_predictions \
            -t $params.max_template_date $fasta_file
    mkdir -p ${params.outdir}/${fasta}/msas
    cp -r ${fasta}/msas ${params.outdir}/${fasta}/msas
    cp ${fasta}/*.pkl ${params.outdir}/${fasta}/
    """
}

process ALPHAFOLD_Inference{
    queue 'gpuq'
    clusterOptions '--gres=gpu:A30:1 --nice'
    errorStrategy 'ignore'
    label 'Alphafold2'
    tag "${fasta}"

    publishDir "${params.outdir}/${fasta}/", mode: 'copy', pattern: "*.pdb"
    publishDir "${params.outdir}/${fasta}/", mode: 'copy', pattern: "*.json"
    publishDir "${params.outdir}/${fasta}/plots", mode: 'copy', pattern: "*.pdf"

    output:
    path("${fasta}/*.pdb")
    path("${fasta}/*.json")
    path("${fasta}/*.pkl")
    path("${fasta}/plots")


    input:
    tuple val(fasta),path(fasta_file),val(preset),val(model_index)

    module 'alphafold/2.3.2'
    script:
    """
    mkdir -p ${fasta}/msas
    cp -r ${params.outdir}/${fasta}/msas ${fasta}/msas
    cp ${params.outdir}/${fasta}/*.pkl ${fasta}/
    alphafold  -o ./ -t $params.max_template_date \
               -g  true \
               -m $preset  \
               -n $model_index \
               -i $params.num_predictions \
               -r $params.model_to_relax \
               $fasta_file
    """
}