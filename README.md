# Feature Barcode QC for SC 

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.com/robinmeyers/featurebarcode-qc-snakemake.svg?branch=master)](https://travis-ci.com/robinmeyers/featurebarcode-qc-snakemake)

This is a workflow built in snakemake to process a feature barcode library from a 10X scRNA-seq experiment. Starting from raw sequencing reads, the pipeline will map to reference of feature barcodes, quantify cell barcodes and UMIs, and generate an analysis report.

## Authors

* Robin Meyers (@robinmeyers)

## Usage

### Simple

#### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/robinmeyers/featurebarcode-qc-snakemake/releases).
If you intend to modify and further extend this workflow or want to work under version control, fork this repository as outlined in [Advanced](#advanced). The latter way is recommended.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).


Install a conda distribution, such as miniconda. On a shared computing cluster using modules, you may be able to run

```$ module load miniconda```

Clone this repositiory

```$ git clone git@github.com:robinmeyers/featurebarcode-qc-snakemake.git```

Create and activate the conda environment

```
$ conda env create -q -f=envs/conda.yaml -n featurebarcode-qc
$ conda activate featurebarcode-qc
```

Run on test data

```$ snakemake --directory .test```


```

#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`. This includes the path to the fasta file of reference feature barcodes, the directory of fastq files, location of plasmid library sequencing, and how targets are specified in the names of the features.

#### Step 3: Execute workflow

Ensure the correct conda environment is active with

```$ conda activate featurebarcode-qc```

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster qsub --jobs 100


# Step 4: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.

### Updating the workflow

If you installed the workflow by cloning the github repo, you can pull latest updates to workflow with 

```$ git pull --rebase```

This will require you to first commit any changes you made to your configuration file before pulling new updates.

Then simply rerun the `snakemake` command.

### Advanced

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for the concrete project/run on your machine.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).
4. Modify the config, and any necessary sheets (and probably the workflow) as needed.
5. Commit any changes and push the project-branch to your fork on github.
6. Run the analysis.
7. Optional: Merge back any valuable and generalizable changes to the [upstream repo](https://github.com/robinmeyers/featurebarcode-qc-snakemake) via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
8. Optional: Push results (plots/tables) to the remote branch on your fork.
9. Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
10. Optional: Delete the local clone/workdir to free space.


## Testing

Tests cases are in the subfolder `.test`. They are automtically executed via continuous integration with Travis CI.
