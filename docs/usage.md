# nf-core/taxtriage: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/taxtriage/usage](https://nf-co.re/taxtriage/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```console
--input '[path to samplesheet file]'
```

### Multiple Samples 

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1,fastq_2,platform,from,trim,sequencing_summary,single_end,barcode
Sample_1,AEG588A1_S1_L001_R1_001.fastq.gz,AEG588A1_S1_L001_R2_001.fastq.gz,ILLUMINA,,FALSE,,FALSE,FALSE
Sample_2,AEG588A1_S2_L001_R1_001.fastq.gz,AEG588A1_S2_L001_R2_001.fastq.gz,ILLUMINA,,FALSE,,FALSE,FALSE
Sample_3,AEG588A1_S3_L001_R1_001.fastq.gz,AEG588A1_S3_L001_R2_001.fastq.gz,ILLUMINA,,FALSE,,FALSE,FALSE
```




### Multiple Samples AND Platforms 

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:



```console
single_end,sample,from,platform,barcode,fastq_1,fastq_2,sequencing_summary,trim
TRUE,NB01,data/test-run/fastq_demux/NB01,OXFORD,TRUE,,,data/test-run/sequencing_summarycovid.txt,
TRUE,NB03,data/test-run/fastq_demux/NB03,OXFORD,TRUE,,,,
TRUE,NB11,data/test-run/fastq_demux/NB11,OXFORD,TRUE,,,,
TRUE,Sample_1,,OXFORD,FALSE,data/test-run/sample_metagenome.fastq.gz,,,
FALSE,ERR6913101,,ILLUMINA,FALSE,data/test-run/ERR6913101_1.fastq.gz,data/test-run/ERR6913101_2.fastq.gz,,
```

### Samplesheet Information

| Column     | Description                                                                                                                                                                            |
| ---------  | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`   | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`  | OPTIONAL (if not using 'from'). Full path to FastQ file for Illumina short reads 1 OR OXFORD reads. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                             |
| `fastq_2`  | OPTIONAL. Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `barcode`  | OPTIONAL. TRUE/FALSE, is the row attributed to a demultiplexed barcode folder of 1 or more fastq files or is it a single file that is .gz?                                                       |
| `from`     | OPTIONAL (if not using fastq1/2) Directory path of the barcode, only used with the column being set as TRUE in the barcode column                                                                                       |
| `platform` | Platform used, [ILLUMINA, OXFORD]                                                            |
| `trim` | TRUE/FALSE, do you want to run trimming on the sample?                                       |
| `sequencing_summary` | OPTIONAL. If detected, output plots based on the the sequencing summary file for that sample | 

An [example samplesheet](../examples/Samplesheet.csv) has been provided with the pipeline alongside some demo data.

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/taxtriage --input samplesheet.csv --outdir <OUTDIR> --demux --asssembly <assembly file> -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work                # Directory containing the nextflow working files
<OUTIDR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### CLI Parameters Possible and Explained

| Parameter     | Description                                                                                                                                                                            |
| ---------  | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input <samplesheet.csv_path>`   | CSV file containing samplesheet information for each sample. See [Samplesheet section](#samplesheet-information)  |
| `--outdir <path_to_output_to_be_made>`  | Output folder.                                             |
| `--skip_assembly`  | Skip the assembly process. Currently this is recommmended as spades/flye are not functioning properly |
| `--skip_fastp`  | TRUE/FALSE,  do not filter with fastp |
| `--skip_plots`  | TRUE/FALSE,  do not make any plots  |
| `--remove_taxids  <numbers[]>`  | "taxidA taxidB..." a list of one or more taxids to remove from the kraken report prior to downstream analysis. Use "'9606'" for human reads |
| `--skip_realignment`  | TRUE/FALSE, Skip realignment step. You will not get a metrics report as a result                    |            
| `--minq <number>`     | What minimum quality would you want in your samples. Disabled if you run --skip_fastp. Default is 7 for Oxford Nanopore, 20 for Illumina         |
| `--trim`     | Remove adapters from data prior to qc filtering. Trimgalore for Illumina, Porechop for Illumina (SLOW)  |
| `--subsample <number>`     | Take a subsample of n reads from each sample. Useful if your data size is very large and you want a quick triage analysis      |
| `--db <path_to_kraken2_database>` | Database to be used. IF `--low_memory` is called it will read the database from the fileystem. If not called, it will load it all into memory first so ensure that the memory available (limited as well by `--max_memory` is enough to hold the database). If using with --download-db, choose from download options {minikraken2, flukraken2} instead of using a path |
| `--download_db` | Download the preset database indicated in `--db` to `--outdir` |
| `--max_memory <number>GB` | Max RAM you want to dedicate to the pipeline. Most of it will be used during Kraken2 steps so ensure you have more memory than the size of the `--db` called |
| `-latest` | Use the latest revision (-r). If you want to autopull the latest commit for a branch from https://github.com/jhuapl-bio/taxtriage, specify this. Used in Basestack and the Cloud (default toggle) |
| `--low_memory <number>` | If you don't have enough memory to load the kraken2 db, call this to read it from the filesystem for each read. THIS IS MUCH SLOWER THAN THE DEFAULT METHOD |
| `--max_cpus <number>` | Max CPUs you want to dedicate to the pipeline | 
| `--demux` | If your Samplesheet contains a folder (rather than 1-2 fastq files), you MUST call this flag | 
| `-resume` | Resume the run from where it left off. IF not called the pipeline will restart from the Samplesheet check each time | 
| `-r [main, stable, etc.]` | Specify the branch/revision name to use if pulling from github (not local main.nf file) | 
| `-profile [docker,singularity,conda]` | Conda, Singularity, or Docker | 


## AWS with Nextflow Tower

While we support Taxtriage in both Basestack and native and local CLI deployment, you can also import and run code from Nextflow Tower. This process can be convoluted in order to get it applicable for cloud environments but, once fully setup, becauses very easy to reproduce at low cost on AWS. Please be aware of sensitivity of data when working in the cloud

1. First, you must create a [Nextflow Tower](https://cloud.tower.nf/) account
    - Further documentation for the following steps can be found [here](https://help.tower.nf/22.2/)
2. Request a user account to brian.merritt@jhuapl.edu
    - Information for accessing the S3 Buckets and creating a credential-specific Compute environment will be included in the email
At this point follow all steps for setting up AWS. You should see some information like setting up an IAM user like so:


![CloudDeck](images/Cloud_AWS_NFTower_only/Slide2.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide3.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide4.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide5.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide6.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide7.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide8.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide9.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide10.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide11.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide12.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide13.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide14.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide15.png)
![CloudDeck](images/Cloud_AWS_NFTower_only/Slide16.png)






<!-- ![IAM Examples](images/iam_example.png) -->
<img src="images/iam_example.png" alt="drawing" height="300"/>

When all necessary items have been setup, you'll need to load your data you want to run into an S3 bucket. Be aware that all filesystem refrences to files/directories need to be relative to the S3 bucket. However, with things like the Samplesheet, the paths can be relative only within the S3 bucket so you don't need to use things like `S3://bucketname/path/to/file `

<img src="images/s3bucket.png" alt="drawing" height="300"/>
<img src="images/s3bucket2.png" alt="drawing" height="300"/>

In the above example, I've made an s3 bucket with necessary permissions for nextflow-tower to run. I have a directory with the same test-data you can find in the `examples` directory in the root level of this bit of code [here](../examples)


Once the AWS system is setup, let's head back to Nextflow Tower. On the left, you can select a drop-down and open up your own personal launchpad (and other features)

<img src="images/cloud_showcase_nftower.png" alt="drawing" height="300"/>

<img src="images/addpipeline.png" alt="drawing" />

<!-- <img src="images/taxtriagelaunchpad1.png" alt="drawing" height="300"/> -->

![Launchpad1](images/taxtriagelaunchpad1.png)

MAKE SURE that the compute environment matches the one you set up when you set your credentials in AWS with NF tower

If you expand the pipeline parameters, you can mimic what I've written for my example with your own paths for the S3 bucket and example data. Note that these are going to be identical to the parameters available at [here](#cli-parameters-possible-and-explained)



![Launchpad2](images/taxtriagelaunchpad2.png)

Next, Run the pipeline from the Launchpad. Be aware that all inputs to the S3 bucket are tailored for my example. You own inputs will vary including things like `low_memory` or database name used

[cloud_run_2.webm](https://user-images.githubusercontent.com/50592701/192596313-7e30f285-dc1d-4c62-99d2-5791a5d8c0e9.webm)

[cloud_run_3.webm](https://user-images.githubusercontent.com/50592701/192596272-46007980-cc07-46c3-978f-e1846adbfffb.webm)

[cloud_run_1.webm](https://user-images.githubusercontent.com/50592701/192596324-57162d50-2738-4b7f-ba8c-2fa473ca1433.webm)

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/taxtriage releases page](https://github.com/nf-core/taxtriage/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```


> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
