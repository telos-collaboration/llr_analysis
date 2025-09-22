# Finite-temperature Yang-Mills theories with the density of states method: towards the continuum limit - Analysis workflow

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16579683.svg)](https://doi.org/10.5281/zenodo.16579683)

The workflow in this repository performs the analyses presented in the paper
[Finite-temperature Yang-Mills theories with the density of states method: towards the continuum limit]().

## Requirements

- Conda, for example, installed from [Miniforge][miniforge]
- [Snakemake][snakemake], which may be installed using Conda

## Setup

1. Install the dependencies above.
2. Clone this repository including submodules
   (or download its Zenodo release and `unzip` it)
   and `cd` into it:

   ```shellsession
   git clone --recurse-submodules https://github.com/telos-collaboration/llr_analysis
   cd llr_analysis
   ```

3. The raw data and metadata can be downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16580109.svg)](https://doi.org/10.5281/zenodo.16580109).
Download `metadata.zip`, `raw_data.zip` and `decompress_raw_data.sh` from Zenodo and decompress the
archives in this directory by invoking `decompress_raw_data.sh`. On slow and/or unstable connections consider using
`wget -c` to download large files. After decompression the raw data takes
up roughly 90GB of space.

## Running the workflow

The workflow is run using Snakemake:

``` shellsession
snakemake --cores 1 --use-conda
```

where the number `1`
may be replaced by
the number of CPU cores you wish to allocate to the computation.

Snakemake will automatically download and install
all required Python packages.
This requires an Internet connection;
if you are running in an HPC environment where you would need
to run the workflow without Internet access,
details on how to preinstall the environment
can be found in the [Snakemake documentation][snakemake-conda].

Using all 16 cores on an AMD 5950x CPU the analysis takes roughly 6 minutes.

``` shellsession
time snakemake --use-conda --cores 16 --forceall
real    6m11.190s
user    46m39.241s
sys     2m39.236s
```

## Output

Output plots, tables, and definitions
are placed in the `assets/plots`, `assets/tables`, and `assets/definitions` directories.

Output data assets are placed into the `data_assets` directory.

Intermediary data are placed in the `intermediary_data` directory.

## Reusability

This workflow is relatively tailored to the data
which it was originally written to analyse.
Additional ensembles may be added to the analysis
by adding relevant files to the `raw_data` directory,
and adding corresponding entries to the files in the `metadata` directory.
However,
extending the analysis in this way
has not been as fully tested as the rest of the workflow,
and is not guaranteed to be trivial for someone not already familiar with the code.

[miniforge]: https://github.com/conda-forge/miniforge
[snakemake]: https://snakemake.github.io
[snakemake-conda]: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html
