  # To Do List
  
- Add documentation
  - indicate the two input options (full metadata for squidbase, only minimal metadata)
  - parameters 
- Should user adapt the config file or are there command line option?
- Run also on plasmodium datasets. Here, the larger datasets appeared to have crashed during the minimap2 step. This is not related to the pipeline nessecarily, but should be investigated


# Log

## Minimal dataset

Now included in github repo

## support for squidbase metadata csv output

Done via additional code in the nextflow process section. Outputting the entire CSV is optional.

## docker & apptainer 

Used the following commands 

- docker build -t wim1994/squidpipe:latest .
- docker push wim1994/squidpipe:latest
- docker pull wim1994/squidpipe:latest
- docker run -it wim1994/squidpipe:latest /bin/bash

Apptainer is enabled via the config file

## new squidpipe yaml 

Made a new yaml, and used the following commands

conda create --name squidpipe_env
conda activate squidpipe_env

- conda install bioconda::kraken2
- conda install bioconda::seqtk
- conda install conda-forge::pigz
- conda install bioconda::seqkit
- conda install bioconda::minimap2
- conda install bioconda::samtools
- pip install h5py
- pip install pod5
- pip install vbz-h5py-plugin
- pip install pandas
- conda install -c conda-forge ncbi-datasets-cli

conda env export -n squidpipe_env > squidpipe_env.yaml

I have rebuilt the docker container after this.

## parameter validation


- conda install nf-core
- conda remove nf-core
- pip install nf-core
- nf-core pipelines schema build

conda create --name nf-core python=3.12 nf-core nextflow
conda activate nf-core

nf-core pipelines schema build