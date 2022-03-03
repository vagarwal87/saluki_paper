# Saluki: Predicting mRNA half-life from mRNA sequence

This repository is intended to accompany our publication, primarily to enhance the reproducibility of our results. For more information please refer to:

Agarwal V, Kelley D. The genetic and biochemical determinants of mRNA degradation rates in mammals. 2022. **_bioRxiv_**.

Each folder is organized as follows:

* Collection and visualization of relationships between human and mouse half-life datasets (as shown in **Fig 1**)
* Analysis of relationships between half-life samples (as shown in **Fig 2**)
* Lasso regression modeling to predict half-life from mRNA sequence features (as shown in **Fig 3**)
* Lasso regression modeling to predict half-life from biochemical features (as shown in **Fig 4**)
* Benchmarking the performance of our deep learning model Saluki, and interpreting what it has learned through insertional motif analysis (as shown in **Fig 5**)
* Benchmarking Saluki's ability to predict the function of non-coding elements and variants therein, as measured by massively parallel reporter assays (as shown in **Fig 6**)

If you find our code or predictions to be helpful for your work, please cite the paper above.

# Dependencies for running entire pipeline:
* Python3 modules: argparse, collections, numpy, h5py, json, glob, logomaker, matplotlib, pandas, seaborn, [basenji](https://github.com/calico/basenji), [tfmodisco](https://github.com/kundajelab/tfmodisco), [deeplift](https://github.com/kundajelab/deeplift)

* R libraries: abind, biomaRt, devtools, gdata, ggplot2, glmnet, imputeR, missMDA, readxl, plyr, LSD, data.table, preprocessCore, psych, Biostrings, rhdf5, tidyr, reshape2, RColorBrewer

* [TensorFlow (>=2.0)](https://www.tensorflow.org/install/)

* TomTom from [The MEME Suite](http://meme-suite.org/doc/download.html?man_type=web)

# Instructions for use

For R code to work properly, please copy the contents of .Rprofile in this folder to your local .Rprofile.

Users are advised to read the code closely and modify commented pieces as appropriate to acquire
desired output for your environment. For example, you will need to download all of the additional
R library and Python module dependencies for the code to work. This being said, if you find crucial
files are missing, making the code unusable, or if you identify a major problem in the code, please
raise a Github issue.

In each Figure's folder, change directories to it and please read the file "runme.sh" first as it provides a general overview of relevant commands that were used sequentially to pre-process the data and generate the figures.

To train the Saluki model, please visit this [link](https://github.com/calico/basenji/tree/master/manuscripts/saluki) within the basenji repo.

**OPTIONAL**: For full functionality and to fix symbolic links, download the associated datapack and save the "datasets" folder in the base Github directory:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6326409.svg)](https://doi.org/10.5281/zenodo.6326409)

The figures will link to this folder accordingly. Some of the files need to be decompressed, and not all files are provided due to minimize the package size. If you need additional files not provided for the purpose of reproduction, please contact Vikram Agarwal (vagar {at} calicolabs {dot} com).
