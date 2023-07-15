# dadi-cli-analysis

## Introduction

This repo contains [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines for replicating the analysis in our `dadi-cli` paper. These pipelines were tested on Linux operating systems (Oracle Linux 8).

To replicate our analysis, users should install [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) first, then use the following commands to create a virtual environment for the analysis.

	conda install mamba -n base -c conda-forge
	mamba env create -f environment.yml
	conda activate dadi-cli-analysis

Besides, users should download `annovar` based on the instruction from https://annovar.openbioinformatics.org/en/latest/ and put it into the directory `ext`.

## Running the pipelines

To replicate our analysis, users should run the four snakemake files in the `workflows` directory step by step.

        snakemake -s workflows/step1_download.smk --profile slurm/
        snakemake -s workflows/step2_annovar.smk --profile slurm/
        snakemake -s workflows/step3_dfes.smk --profile slurm/
        snakemake -s workflows/step4_plots.smk --profile slurm/

In the `slurm` directory, a `config.yaml` file is provided. Users can adjust the parameters in this file based on the settings of their clusters.

Log information can be found in the `logs` directory when running the pipelines. After running the pipelines, users can find the figures of our results in the `results/plots` directory.
