<div align="center">

# CoCo-ST: Compare and Contrast Spatial Transcriptomics

**A scalable contrastive learning framework for identifying spatial domains in spatial transcriptomics data**

[![GitHub Stars](https://img.shields.io/github/stars/WuLabMDA/CoCo-ST?style=social)](https://github.com/WuLabMDA/CoCo-ST/stargazers)
[![GitHub Issues](https://img.shields.io/github/issues/WuLabMDA/CoCo-ST)](https://github.com/WuLabMDA/CoCo-ST/issues)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/WuLabMDA/CoCo-ST/blob/master/LICENSE)

</div>

---

## Overview

**CoCo-ST (Compare and Contrast Spatial Transcriptomics)** is a computational algorithm that leverages **contrastive learning** to identify both **high-variance** and **low-variance** spatial structures across spatial transcriptomics datasets.  
Unlike conventional methods that primarily detect dominant spatial domains, CoCoST uncovers subtle yet biologically meaningful spatial niches that are critical for understanding **early tumor evolution**, **precancerous changes**, and **cell–cell interactions**.

---

<p align="center">
  <img src="figures/workflow.png" width="900">
</p>

## Key Features

-  **Contrastive spatial learning** – disentangles unique and shared structures across datasets (e.g., tumor vs. normal).  
-  **Multi-scale domain detection** – supports analysis at subcellular, cellular, and tissue (spots) scales.  
-  **Cross-sample integration** – harmonizes data from multiple samples and platforms (Visium, Visium HD, Xenium).  
-  **Scalable** – efficiently handles millions of spatial spots.  
-  **Biologically interpretable** – links contrastive components and spatial domains to known biological structures.  

---

## Installation

You can install the latest development version of **CoCoST** directly from GitHub:

```r
# Option 1: using devtools
# install.packages("devtools")
devtools::install_github("WuLabMDA/CoCo-ST")
```
or by running:

```r
# Option 2: using remotes
# install.packages("remotes")
remotes::install_github("WuLabMDA/CoCo-ST")

```
### Support

For any question, request or bug report please create a new issue in this repository. 

### Contributions

We welcome contributions and suggestions from the community. If you have any idea, please submit it as an issue, which we will look into and ask for further explannation if necessary. 

### Reporting issues

CoCo-ST is under continuous development. If you encounter an issue, please make sure you install all the required packages necessary to run the codes. If that does not solve your problem, please open a new issue detailing your encountered problem by providing a code and a demo example. We will try to look into your issue and hopefully provide you with a solution. Thanks.

##  Cite CoCo-ST
If you use this work, please cite:
/n Aminu, M., Zhu, B., Vokes, N. et al. CoCo-ST detects global and local biological structures in spatial transcriptomics datasets. Nat Cell Biol (2025). https://doi.org/10.1038/s41556-025-01781-z



