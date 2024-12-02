Pull an nf-core pipeline without cloning:

nextflow pull nf-core/demo

Pulled pipelines are stored in a hidden assets folder. By default, this folder is $HOME/.nextflow/assets, but in this training environment the folder has been set to $NXF_HOME/assets.

$NXF_HOME was not set in my local env, was in $HOME/.nextflow/assets.

nextflow list

Now let's run the test profile (see conf test.config for references to nf-core/test-datasets), for this pipeline:

nextflow run nf-core/demo -profile docker,test --outdir results

We need to manually specify an output dir. Also unclear how this test.config maps to a profile, this is not how we set up a profile previously??

See this in the nextflow.config profile spec:

    test      { includeConfig 'conf/test.config'      }

If there are nextflow version mismatches you can run with a specific version using the NXF_VER env var prefix e.g.

NXF_VER=24.09.2-edge nextflow run nf-core/demo -profile docker,test --outdir results


### Creating a pipeline

nf-core was not installed in my local env:

    conda install nf-core

Sigh:

Error while loading conda entry point: conda-libmamba-solver (dlopen(/Users/nathanjohnson/miniconda3/lib/python3.12/site-packages/libmambapy/bindings.cpython-312-darwin.so, 0x0002): Library not loaded: @rpath/libarchive.20.dylib                                                                                                                                                              
  Referenced from: <496442DC-0EDE-3705-A2B5-401A4FC0D733> /Users/nathanjohnson/miniconda3/lib/libmamba.2.0.0.dylib                                                                               
  Reason: tried: '/Users/nathanjohnson/miniconda3/lib/libarchive.20.dylib' (no such file), '/Users/nathanjohnson/miniconda3/lib/python3.12/site-packages/libmambapy/../../../libarchive.20.dylib' (no such file), '/Users/nathanjohnson/miniconda3/lib/python3.12/site-packages/libmambapy/../../../libarchive.20.dylib' (no such file), '/Users/nathanjohnson/miniconda3/bin/../lib/libarchive.20.dylib' (no such file), '/Users/nathanjohnson/miniconda3/bin/../lib/libarchive.20.dylib' (no such file), '/usr/local/lib/libarchive.20.dylib' (no such file), '/usr/lib/libarchive.20.dylib' (no such file, not in dyld cache))                                                                  

Attempting to address the underlying mamba reqs by ensuring they are from the same channel eventaully worked, but it took a while:

    conda install --solver=classic conda-forge::conda-libmamba-solver conda-forge::libmamba conda-forge::libmambapy conda-forge::libarchive 

Reported: https://github.com/nf-core/tools/issues/3303

No need to mess with pipx/pyenv i.e. export PIPX_DEFAULT_PYTHON="$HOME/.pyenv/versions/3.12/bin/python"

conda borks the pyenv anyway :/

✗ which python
/Users/nathanjohnson/miniconda3/bin/python


nf-core --help

Create pipeline from base template:
nf-core pipelines create

Nice G/TUI!
- Custom


This section probably need more work:
https://training.nextflow.io/hello_nextflow/09_hello_nf-core/#check-the-input-data

The concept of `take` hasn't been introduced, nor is it clear where `ch_samplesheet` is coming from or where it is being used.

https://www.nextflow.io/docs/latest/workflow.html#workflow-inputs-take

The inline comment does explain this, but may a reference to the docs from the tutorial might help here?

Let's add a module:

    nf-core modules list remote

Ensure we are in the right dir:

    ✗ pwd
    /Users/nathanjohnson/src/nextflow-io/training/hello-nextflow/hello-nf-core/nf-core-myfirstpipeline

nf-core modules install seqtk/trim

    INFO     Installing 'seqtk/trim'          
    INFO     Use the following statement to include this module:                          
    
    include { SEQTK_TRIM } from '../modules/nf-core/seqtk/trim/main'     

Note: can also omit the tool and follow prompts.

You can check the hashes for installed modules here:

    ✗ cat modules.json
    {
        "name": "nf-core/myfirstpipeline",
        "homePage": "https://github.com/nf-core/myfirstpipeline",
        "repos": {
            "https://github.com/nf-core/modules.git": {
                "modules": {
                    "nf-core": {
                        "multiqc": {
                            "branch": "master",
                            "git_sha": "cf17ca47590cc578dfb47db1c2a44ef86f89976d",
                            "installed_by": ["modules"]
                        },
                        "seqtk/trim": {
                            "branch": "master",
                            "git_sha": "666652151335353eef2fcd58880bcef5bc2928e1",
                            "installed_by": ["modules"]
                        }
                    }
                },
                "subworkflows": {
                    "nf-core": {
    ...

Look in the main.nf for seqtk to see the inputs or more helpful:

    nf-core modules info seqtk/trim


               ╷                                                                                                                                                                               ╷             
 📥 Inputs     │Description                                                                                                                                                                    │     Pattern 
╺━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
 input[0]      │                                                                                                                                                                               │             
╶──────────────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼────────────╴
  meta  (map)  │Groovy Map containing sample information e.g. [ id:'test', single_end:false ]                                                                                                  │             
╶──────────────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼────────────╴
  reads  (file)│List of input FastQ files                                                                                                                                                      │*.{fastq.gz} 
     

The channel ch_samplesheet used by the FASTQC process can be used as the reads input.

Added some module.config to set a hardcode param for seqtk

Now let's expose a param to enable skipping that step via nextflow.config.

This should now work:
    nextflow run . -profile test,docker --outdir results --skip_trim

Weird:

    WARN: The following invalid input values have been detected:

    * --skip_trim: true


    executor >  local (1)
    [38/207a7b] process > NFCORE_MYFIRSTPIPELINE:MYFIRSTPIPELINE:MULTIQC [100%] 1 of 1 ✔

It warns, but clearly skips??

This is also in the tutorial but not explained.

Ah..the next bit...


This is due to lack of representation in the nextflow_schema.json.

Suggestion is to use a web UI tool do this rather than manually editing due to complexity.

    nf-core pipelines schema build

Whoa, nice :)

Help test modal initilal focussed in text area input, but clicking on it removes focus which can't be returned unless closing/saving and re-opening.

Now let's create a custom module for FastQE

    nf-core modules create


After updating the template and plumbing in myfirstpipeline.nf, re-run the pipeline and look in results/fastqe for some smiley action.
