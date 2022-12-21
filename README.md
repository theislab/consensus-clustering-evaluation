# Multi-resolution consensus clustering evaluation

## Introduction

This repository contains a pipeline for comparing single-cell RNA sequencing (scRNA-seq) clustering methods used to evaluation the [multi-resolution consensus clustering method](https://github.com/theislab/multires-consensus-clustering) (MRCC).
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It was based on the [nf-core](https://nf-co.re/) template (although most of the bells and whistles have been removed).

### Directory structure

- `bin/` - Scripts for stages of the pipeline
- `conf/` - Configuration files
- `datasets/` - Scripts for creating dataset input files. The created files should also be placed here for the pipeline to run.
- `envs/` - Conda environment YAML files for pipeline stages
- `workflows/` - NextFlow files specifying the pipeline
- `LICENSE` - MIT license
- `main.nf` - The main NextFlow pipeline file
- `nextflow.config` - The main NextFlow configuration file
- `README.md` - This README

## Pipeline summary

The pipeline follows a standard benchmarking format of datasets, methods and metrics.
Each dataset is prepared and passed to each method and the results of the method are scored using the metrics.
This is done in a combinational way so the number of stages is roughly `number of datasets` x `number of methods` x `number of metrics`.

### Datasets

The pipeline uses a selection of real and simulated datasets for evaluation.
Scripts for creating the datasets are available in `datasets/` but are not part of the pipeline and need to be run in advance.
#### Real datasets

The pipeline uses several real datasets for evaluation, each of these has pre-existing curated labels.
We have selected a single sample from each dataset in order to avoid needing to consider batch effects when running the methods.

- [Azimuth kidney (CZIKidney763280)](https://azimuth.hubmapconsortium.org/references/#Human%20-%20Kidney) - CZIKidney763280 sample from the Azimuth example kidney query dataset
- [Azimuth lung (Dropseq-2)](https://azimuth.hubmapconsortium.org/references/#Human%20-%20Lung%20v1) - Dropseq-2 sample from the Azimuth example lung query dataset
- [Azimuth bone marrow (batch2)](https://azimuth.hubmapconsortium.org/references/#Human%20-%20Bone%20Marrow) - batch2 sample from the Azimuth example bone marrow query dataset
- [Azimuth mouse motor cortex (352357)](https://azimuth.hubmapconsortium.org/references/#Human%20-%20Bone%20Marrow) - 352357 sample from the Azimuth example mouse motor cortex query dataset
- [COMBAT (S00052-Ja005E-PBCa)](https://zenodo.org/record/6120249) - S00052-Ja005E-PBCa sample from the COvid-19 Multi-omics Blood ATlas (COMBAT) study
- [GCA (F73-FPIL-0-SC-1)](https://gutcellatlas.cellgeni.sanger.ac.uk/) - F73-FPIL-0-SC-1 sample from the Gut Cell Atlas (GCA)
- [NeurIPS (site2-donor1)](https://openproblems.bio/neurips_docs/data/dataset/) - site2-donor1 sample from the Open Problems in Single-cell Analysis NeurIPS 2021 multimodal integration challenge

#### Simulated datasets

Simulated datasets were generated using the [{splatter} package](https://bioconductor.org/packages/splatter/) and the scripts in `datasets/`.

- Simulation (blob) - Simulation with one cell group (i.e. a single cluster)
- Simulation (groups) - Simulation with multiple cell groups (i.e. multiple clusters)
- Simulation (path) - Simulation with a continuous transition between two cell types
- Simulation (rare) - Simulation with multiple cell groups where some of those groups have low occurrences (1-5%)

### Methods

Both standard scRNA-seq processing workflows and scRNA-seq consensus clustering methods were selected for comparison.

- [MRCC](https://github.com/theislab/multires-consensus-clustering) - This is the primary method evaluated by the pipeline. It is run in two stages. First multiple clusterings are performed using a standard method. Second those clusterings are combined using the newly developed consensus approach. This allows to us to test multiple parameter sets for the combining stage without having to repeat the more computationally intensive clustering stage see [Usage](#usage).
- Random - Random assignment of labels as a negative control. The number of labels is the same as the number of labels in the dataset.
- [SC3](https://bioconductor.org/packages/release/bioc/html/SC3.html) - Single-Cell Consensus Clustering. A consensus clustering method designed for scRNA-seq data. The dataset is clustered multiple times using k-means on different distance metrics and transformations of the data.
- [Scanpy](https://scanpy.readthedocs.io/) - Scanpy is the most used Python toolbox for scRNA-seq analysis. The standard Scanpy workflow makes use of graph-based clustering and is comparable to the Seurat workflow.
- [Seurat](https://satijalab.org/seurat/) - Seurat is the most used R toolbox for scRNA-seq analysis. The standard workflow makes use of graph-based clustering and is comparable to the Scanpy workflow.
- [SIMLR](https://bioconductor.org/packages/release/bioc/html/SIMLR.html) - Consensus clustering based on optimising different distance kernels.

### Metrics

Metrics are divided into two categories: unsupervised metrics which compare clustering assignments to the ground truth labels but do not require them to be matched and supervised metrics which treat the task as a classification problem and require clustering assignments to be matched to the ground truth labels.

#### Unsupervised metrics

- [Adjusted Mutual Information](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_mutual_info_score.html)
- [Adjusted Rand Index](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html)
- [Completeness score](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.completeness_score.html)
- [Element-Centric Clustering Similarity](https://doi.org/10.1038/s41598-019-44892-y) (Implemented in the [ClustAssess package](https://github.com/Core-Bioinformatics/ClustAssess))
- [Fowlkes-Mallows Index](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.fowlkes_mallows_score.html)
- [Homogeneity score](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.homogeneity_score.html)

#### Supervised metrics

- [F1 score](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html)
- [Matthews Correlation Coefficient](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.matthews_corrcoef.html)

## Setup

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)
2. Install [`Conda`](https://conda.io/miniconda.html)
3. Download the pipeline by cloning the repository or as a [ZIP file](https://github.com/theislab/consensus-clustering-evaluation/archive/refs/heads/main.zip)
4. Run the scripts in `datasets/` to create input dataset files
## Usage

The pipeline can be run using:

```sh
nextflow run main.nf
```

By default this will just run a small test dataset.
To run on the full datasets a parameters file needs to be provided (see [Parameters](#parameters)).
For example, to run on all datasets using the provided parameters file use:

```sh
nextflow run main.nf -params-file conf/all-datasets.yml
```

To run the pipeline on a high-performance computing system with a submission queue you need to supply a profile configuration.
An example for the HMGU slurm cluster is provided but you should refer to the [NextFlow docs](https://www.nextflow.io/docs) for how to design a profile for your system.

```sh
nextflow run main.nf -profile hmgu-slurm
```
### Parameters

The parameters file can be used to define both datasets and parameter sets for the MRCC method.

#### Datasets

Datasets are defined using the following YAML:

```yaml
input:
    - name: Dataset1       # Shouldn't have spaces, other unusual symbols
      file: Dataset1.h5ad  # Path to a H5AD file containing the dataset
      labels: Labels       # Name of the `.obs` column containing cell labels
    - name: Dataset2
      file: Dataset2.h5ad
      labels: CellLabels
```

#### MRCC parameters

Parameters sets for the MRCC method are defined using the following YAML.
See the `bin/method-mrcc.py` script for more description of the parameters.

```yaml
mrcc:
    - name: MRCC_N_Le_SR_1_MR_1    # Name of the parameter set
      graph_type: neighbour        # Method for building the multi-resolution graph
      community_type: leiden       # Community detection method
      single_resolution: 1         # Community detection resolution for single-resolution graphs
      multi_resolution: 1          # Community detection resolution for the multi-resolution graph
    - name: MRCC_A_Le_SR_1_MR_1
      graph_type: all
      community_type: leiden
      single_resolution: 1
      multi_resolution: 1
```

## Output

Output of the pipeline will be created in the `results/` directory.
This includes the clustering output from each method, the metric scores and some basic summary plots.
The pipeline trace (runtime etc.) is also available in the `pipeline_trace/` directory.
