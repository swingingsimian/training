#!/usr/bin/env nextflow

//params.greeting = "Holà mundo!"
// You can put the parameter declaration inside the workflow block if you prefer. Whatever you choose, try to group similar things 
// in the same place so you don't end up with declarations all over the place.
params.input_file = "data/greetings.csv"


/*
 * Use echo to print 'Hello World!' to standard out
 */
process sayHello {
    publishDir 'results', mode: 'copy' // By default these are symlinks. Also overwrite mode available, this seems to be the default!
    // Files are copied async and may not be immediately available after process, hence downstream processes should not access directly but via channels 
    // Note: Newer workflow level public API here: https://www.nextflow.io/docs/latest/workflow.html#publishing-outputs


    input:
        val greeting // See above for default

    output:
        //path 'output.txt' // This isn't required to actually write the output, but will be required for wiring into another process?
        // The warning in the tutorial is wrong, it does not break in this context
        // https://training.nextflow.io/hello_nextflow/02_hello_world/#33-run-the-workflow-again
        path "${greeting}-output.txt"
        // This cascades to following processes


    script:
    // """
    // echo '$greeting' > output.txt
    // """
        """
        echo '$greeting' > '$greeting-output.txt'
        """
}

process convertToUpper {
    publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}" // single quotes will not interpolate here, obvs!
        
    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
    """
}

// Nice! Preview DAG with nf lang server
workflow {

    // emit a greeting
    //greeting_ch = Channel.of('Hello World!')
    // greeting_ch = Channel.of(params.greeting) // nextflow ./hello-world.nf --greeting "Hello from the cmdline"
    // greeting_ch = Channel.of('Hello', 'Bonjour', 'Hola') // nextflow ./hello-world.nf --greeting "Hello from the cmdline"
    // Note -single-dash-style-params are nextflow builtins/settings whilst --double-dash-style-params are pipeline params 
    greeting_ch = Channel.fromPath(params.input_file)
        .splitCsv() //.view{ "After splitCsv: $it" }
        .flatten() //.view{ "After flatten: $it" }
    // fromPath takes varargs rather than a list, hence the flatten
    // view is optional, but useful

    // This on it's own writes to a single dir, hence overwriting the output we need to use:
    //  ./hello-world.nf -ansi-log false 
    // But this only applied to the default workdirs, the publishDir appears to suffer from the same issue
    // Once unique filenames are implemented no real need for -ansi-log false
    sayHello(greeting_ch)
    convertToUpper(sayHello.out)
}
