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

    input:
    tuple val(fasta),path(fasta_file),val(preset)

    output:
    tuple val(fasta),path(fasta_file),val(preset),path("${fasta}/*.pkl"),path("${fasta}/msas/*"),  emit:output
    
    script:
    
    """
    mkdir -p /vast/scratch/users/$USER/tmp
    export SINGULARITY_BINDPATH="/vast/scratch/users/$USER/tmp:/tmp"
    alphafold -f -o ./  -m $preset \
            -i $params.num_predictions \
            -t $params.max_template_date $fasta_file
    """
}

process ALPHAFOLD_Inference{
    queue 'gpuq'
    clusterOptions '--gres=gpu:A30:1 --nice'
    errorStrategy 'ignore'
    label 'Alphafold2'
    tag "${fasta}"

    publishDir "${params.outdir}/", mode: 'copy', pattern: "${fasta}/ranked_0.pdb", saveAs: { filename -> "${fasta}.pdb" }

    output:
    path("${fasta}/*.pdb")

    input:
    tuple val(fasta),path(fasta_file),val(preset),path(pkl),path(msas),val(model_index)

    script:
    """
    mkdir -p ${fasta}
    cp features.pkl ${fasta}/
    alphafold  -u -o ./ -t $params.max_template_date \
               -g  true \
               -m $preset  \
               -n $model_index \
               -i $params.num_predictions \
               -r $params.model_to_relax \
               $fasta_file
    """
}
