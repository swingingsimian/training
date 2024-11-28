nf-test needs to be installed to run this locally:
- https://www.nf-test.com/installation/
- TLDR; `curl -fsSL https://get.nf-test.com | bash`
- `sudo mv nf-test /usr/local/bin`



nt-test init
nf-test generate process modules/local/samtools/index/main.nf

Follow convention:
mv tests/modules/local/samtools/index/main.nf.test modules/local/samtools/index/tests/

& update script path in test.

nf-test test --profile docker_on modules/local/samtools/index/tests/main.nf.test

Use --update-snapshot to overwrite invalid/old snapshots. Careful with this though!

    Warning: every snapshot that fails during this test run is re-record.

TODO Can we run with individual tests to avoid corrupting other snapshots?
Maybe? There seems to be support for running using a hash, but not the test name?
https://www.nf-test.com/docs/running-tests/


A test with input from a previous step:

nf-test generate process modules/local/gatk/haplotypecaller/main.nf

best practise location:
mkdir -p modules/local/gatk/haplotypecaller/tests
mv tests/modules/local/gatk/haplotypecaller/main.nf.test modules/local/gatk/haplotypecaller/tests/

and update the script location.

nf-test generate process modules/local/gatk/jointgenotyping/main.nf

Let's use pre-generated data as the inputs this time.


Now let's look at a workflow level test.


nf-test generate pipeline main.nf

Add pipeline config to test dir, update to use sample_bams.local.txt and run...

nf-test test --profile docker_on tests/main.nf.test

Now update the testDir nt-test.config to allow running all the tests:

nf-test test --profile docker_on