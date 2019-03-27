# spl_aggregation_sampling_dsgs

This repository contains all the code related with the letter "_Sampling and Reconstruction of Diffused Sparse Graph Signals from Successive Local Aggregations_", written by _Samuel Rey, Fernando J. Iglesias, Cristobal Cabrera and Antonio G. Marques_, and submitted to the _IEEE Signal Processing Letters_.

Some of the experiments contained in this repository use real data which was originally generated in [1], so we want to thank them in advance. We have loaded and processed these data for expressing them with a diffused sparse model when possible.

## Repository organization
This repository is organized as follows:
* **dataset**: processed dataset used in the experiment shown on the figure 1.b of the letter, which uses real-world data. 
* **dataset_additional**: additional useful datasets. It contents the dataset in its original format as obtained from [here](https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets), the complete dataset in MATLAB format after reading the original dataset with the code presented in this repository, and other variants of the processed dataset which have also been tested in the experiment.
* **tools**: includes different scripts and functions which have been used for reading and process the original dataset, and check if it admitted a diffused sparse representation.
* **utils**: folder with the functions used in the main scripts check_node_influence.m and real_data_example.m.
* **check_node_influence.m**: script for running the first experiment of the letter, which is shown in the figure 1.a.
* **real_data_example.m**: script for running the second experiment of the letter, which is shown in the figure 1.b.

## Dependencies
For running the different scripts of this repository it is neccessary to have installed [CVX](http://cvxr.com/cvx/) and [The Graph Signal Processing Toolbox](https://epfl-lts2.github.io/gspbox-html/).

## References
[1] P. D. Dobson and A. J. Doig, “Distinguishing enzyme structures from non-enzymes without alignments,” Journal of molecular biology, vol. 330, no. 4, pp. 771–783, 2003.
