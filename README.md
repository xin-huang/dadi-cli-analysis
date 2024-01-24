# dadi-cli-analysis

## Introduction

This repo contains [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines for replicating the analysis in our `dadi-cli` paper. These pipelines were tested on Linux operating systems (Oracle Linux 8).

To replicate our analysis, users should install `mamba` first. According to the latest [guidelines](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html), it is recommended to install `mamba` via `mambaforge` or `miniforge`, then use the following commands to create a virtual environment for the analysis.

```
mamba env create -f environment.yml
conda activate dadi-cli-analysis
```

Besides, users should download `annovar` based on the instruction from https://annovar.openbioinformatics.org/en/latest/ and put it into the directory `ext`.

## Running the pipelines

To replicate our analysis, users should run the four snakemake files in the `workflows` directory step by step.

```
snakemake -s workflows/step1_download.smk --profile slurm/
snakemake -s workflows/step2_annovar.smk --profile slurm/
snakemake -s workflows/step3_dfes.smk --profile slurm/
snakemake -s workflows/step4_plots.smk --profile slurm/
```

In the `slurm` directory, a `config.yaml` file is provided. Users can adjust the parameters in this file based on the settings of their clusters.

If users want to run the pipelines locally, please use the `-c` option, for example,

```
snakemake -s workflows/step1_download.smk -c 1
snakemake -s workflows/step2_annovar.smk -c 1
snakemake -s workflows/step3_dfes.smk -c 1
snakemake -s workflows/step4_plots.smk -c 1
```

`-c` specifies the number of threads and `snakemake` could run jobs parallelly as many as possible with the given number of threads.

Log information can be found in the `logs` directory when running the pipelines. After running the pipelines, users can find the figures of our results in the `results/plots` directory.
