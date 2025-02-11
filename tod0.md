- make a minimal dataset using squidbase's POD5 files
- fix the mapping issue - what nf operator can be use to "remove" the duplicate taxonommic ids?
- tests
- rerun squidbase dataset
- run also on plasmodium


# support for squidbase output CSV

Columns required: 



# docker container support (+ apptainer support?)

pipeline still works with the conda environment, but when used with the docker container it crashes during the awk step. Why?

rebuild without starting an interactive bash shell in the dockerfile

- docker build -t wim1994/squidpipe:latest .

- docker push wim1994/squidpipe:latest


docker pull wim1994/squidpipe:latest
docker run -it wim1994/squidpipe:latest /bin/bash

# nieuwe squidpipe yaml 


conda create --name squidpipe_env
conda activate squidpipe_env


// conda install bioconda::kraken2

// conda install bioconda::seqtk
// conda install conda-forge::pigz

// conda install bioconda::seqkit
// conda install bioconda::minimap2

// conda install bioconda::samtools




pip: 

pip install h5py
pip install pod5
pip install vbz-h5py-plugin
pip install pandas
conda install -c conda-forge ncbi-datasets-cli

h5py==3.10.0
      - iso8601==2.1.0
      - lib-pod5==0.3.15
      - more-itertools==10.4.0
      - pod5==0.3.15
      - polars==0.20.31
      - pyarrow==16.1.0
      - tqdm==4.66.5
      - vbz-h5py-plugin==1.0.1

pod5==0.3.15

import csv
import os
import subprocess
import zipfile
import shutil
import logging
import json
import argparse

import pandas as pd

conda env export -n squidpipe_env > squidpipe_env.yaml

# issue

only sarcov data outputted 