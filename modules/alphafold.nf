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
    tuple val(fasta), path(fasta_file),val(preset)


    module 'alphafold/2.3.2'
    script:
    
    """
    alphafold -f -o $params.outdir  -m $preset \
            -i $params.num_predictions \
            -t $params.max_template_date $fasta_file
    """
}

process ALPHAFOLD_Inference{
    queue 'gpuq'
    clusterOptions '--gres=gpu:A30:1 --nice'
 
    label 'Alphafold2'
    tag "${fasta}"

    input:
    tuple val(fasta),path(fasta_file),val(preset),val(model_index)

    module 'alphafold/2.3.2'
    script:
    """
    alphafold  -o $params.outdir -t $params.max_template_date \
               -g  true \
               -m $preset  \
               -n $model_index \
               -i $params.num_predictions \
               -r $params.model_to_relax \
               $fasta_file
    """
}