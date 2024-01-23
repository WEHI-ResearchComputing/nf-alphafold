process Test {

    label 'Test'
    
    tag "${filename}"
    
    input:
    tuple val(filename),path(file)
    
    output:
    stdout

    script:
    """
    #!/usr/bin/env python
    print("File input is ${filename}")

    """
}

SPLIT = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

process SplitSequences {

    publishDir "${params.outdir}", mode: 'copy'
    input:
    path('input.fa')

    output:
    path('seq_*')
    """
    $SPLIT input.fa '%^>%' '/^>/' '{*}' -f seq_
    """

}
process Reverse {

    input:
    path(x) 
    
    output:
    stdout 

    """
    cat $x | rev
    """
}