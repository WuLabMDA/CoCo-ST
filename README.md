# PanDA
PanDA is a joint discriminant analysis method aimed at fusing multi-omics datasets through finding a discriminant common latent space. PanDA captures cross-omics interaction and consistency and uses an uncorrelated constraint to ensure that the extracted latent components for each omics (omics-specific components) are not highly correlated. Because the components extracted using PanDA contain valuable discriminant information, we refer to them as discriminant components. These components can be used as inputs to several multi-omics analysis tools to enable efficient, improved downstream analysis. Here and in our paper, we demonstrated the advantages of PanDA over ten integrative multi-omics methods through four distinct downstream analyses: single-cell multi-omics data visualization, patient (or tumor) classification, biomarker identification, and clinical outcome prediction.

<p align="center">
  <img src="figures/workflow.png" width="900">
</p>
