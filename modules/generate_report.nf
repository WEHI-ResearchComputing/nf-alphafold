process Generate_Report {
    label 'Report'
    
    tag "${fasta}"
    

    publishDir "${params.outdir}/", mode: 'copy', pattern: "*.html"

    input:
    tuple val(fasta),val(model_preset), path(pdb), path(pkl), path(features),path(template)

    output:
    tuple val(fasta), path ("*report.html"), emit: report


    
    script:
    """
    generate_plots.py --html_template ${template} --pdb ${pdb.join(' ')} \
                      --name ${fasta} --model_preset ${model_preset}
    """
}
