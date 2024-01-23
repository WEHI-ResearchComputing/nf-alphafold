println "*****************************************************"
println "*  Nextflow <name> pipeline                         *"
println "*  A Nextflow wrapper pipeline                      *"
println "*  Written by <------------------------>            *"
println "*  <-------------Email---------->                   *"
println "*                                                   *"
println "*****************************************************"
println " Required Pipeline parameters                        "
println "-----------------------------------------------------"
println "Input  Directory   : $params.inputdir                   "
println "Output Directory   : $params.outdir                  " 
println "*****************************************************"

// include modules
include {  Test           } from './modules/test.nf'
include {  SplitSequences } from './modules/test.nf'
include {  Reverse        } from './modules/test.nf'


//Start main workflow
workflow {

    Channel.fromPath(params.inputdir+"/*.fa",checkIfExists:true).set{input_ch}
    
    split_ch = SplitSequences(input_ch) 
    Reverse(split_ch).view()


    input_ch.map{ file -> tuple( file.baseName , file)}.set{mapped_input_ch}
    Test(mapped_input_ch).view()

}
