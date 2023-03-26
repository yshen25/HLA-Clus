# HLA-Clus
> **Package for clustering HLA alleles based on structural landscape of peptide binding groove**

## Introduction
The pipeline contains **three** major steps:  
  1. To process the structures as cloud of points.  
  2. To calculate the structure distances (SDs) between HLA alleles, which is demonstrated to highly correlate with peptide binding specificity.  
  3. Two approaches are provided to cluster HLA alleles based on SD: hierarchical clustering and nearest-neighbor clustering.  
    * The **hierarchical clustering** method is used to classify a group of user-defined HLA alleles  
    * The **nearest-neighbor clustering** method is used to classify user-defined alleles according to pre-defined supertypes  

## Installation
### 1.Prerequisites
HLAcc relies on the following public available packages
 1. [Python](https://www.python.org/) = 3.8
 2. [Jupyter Notebook](https://jupyter.org/)
 3. [NumPy](https://numpy.org/)
 4. [Scipy](https://scipy.org/)
 5. [pandas](https://pandas.pydata.org/)
 6. [BioPandas](http://rasbt.github.io/biopandas/)
 7. [PyMOL](https://pymol.org) >= 2.5.2  
 You may need a lisence to use PyMOL package, please refer to [https://pymol.org/2/buy.html](https://pymol.org/2/buy.html)  
 8. [Biopython](https://biopython.org/)
 9. [matplotlib](https://matplotlib.org/)
 10. [seaborn](https://seaborn.pydata.org/)
 11. [scikit-learn](https://scikit-learn.org/)

### 2.Installation
 * via **git**
 
 This method is recomended for most users\
 First, make sure all required packages are installed\
 Then, go to the directory that you would like to put the HLAcc directory
 ```shell
 cd [destination directory]
 ```
 Finally, clone the whole repository
 ```shell
 git clone https://github.com/yshen25/HLAcc
 ```
 * via **pip**
 
 This method is for the users who would like to customize functions provided in HLAcc. In this way, HLAcc will be installed as a standard pip library, and is callable via import
 First, install PyMOL using conda, as it is not supported by pip install
 ```shell
 conda install -c conda-forge -c schrodinger pymol-bundle
 ```
 Additionally, on macOS, PyQt5 need to be installed:

 ```shell
 pip install PyQt5
 ```
 Then, you can install HLAcc via pip, other required packages will also be installed automatically if you haven't done so.
 ```shell
 pip install git+https://github.com/yshen25/HLAcc
 ```
 
## Usage
 0. Structural modeling  
 HLAcc takes 3D HLA structures in **.pdb** format as input. We have modeled 451 HLA I structures using ColabFold and Rosetta FastRelax, which could be found [here].  
 You could also use costomized HLA structures.  
 1. Structure processing  
 Put all structures in one directory, then specify a destination directory for coarse grained structures. Then use **Processing_pipeline** function to accomplish the task:  
 ```python
 from HLAcc.pipeline import Processing_pipeline
 
 Processing_pipeline([Structure_directory], [Destination_directory])
 ```
 2. Hierarchical clustering
 To hierarchically cluster HLA alleles stored in Destination_directory into N clusters, use **HC_pipeline** function as follows:  
 ```python
 from HLAcc.pipeline import HC_pipeline
 
 HC_pipeline([Destination_directory], [N])
 ```
 3. Nearest neighbor clustering
 To cluster HLA alleles stored in Destination_directory based on clusters defined in anchor_dictionary, use **NN_pipeline** function:
 ```python
 from HLAcc.pipeline import NN_pipeline
 
 NN_pipeline([Destination_directory], [anchor_dictionary])
 ```
