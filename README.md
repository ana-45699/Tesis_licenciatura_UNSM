# Tesis_licenciatura_UNSM

##1. Introduction

Gene fusions are hybrid genes produced by chromosomal recombination or aberrant
alternative splicing. Activated kinase fusion genes have been found in multiple types of
cancer. There are several protein kinases used as tumoral biomarkers due to their high
mutation frequency in cancer patients. Furthermore, these kinase fusion genes have
been found to be intrinsically related to processes such as growth regulation and cell
replacement. The main element of the protein kinase function is the kinase domain which
is normally controlled by a regulatory domain of the protein kinase itself. Acral lentiginous
melanoma (ALM) is a poorly studied subtype of cutaneous melanoma. Currently, the
number of fusion genes found in ALM is relatively low and their association with this type
of cancer is unknown.

##2. Objective

Detect and identify the fusion genes with potential to be tumor markers in tumor tissue of patients with ALM

##3. Analysis or Workflow Steps

* Quality analysis and determine the parameters to run fusion genes predictors
* Detection of the fusion genes by STAR-Fusion and Fuseq
* Exploratory analysis of the fusion genes in all the samples
* Comparative analysis between the results of Fuseq and STAR-Fusion
* Search of protein kinase genes inside the fusion genes,that are in both results (Fuseq and STAR-Fusion)
* Validation with FusionInspector of the selected fusion genes
* Prediction of functional domains of the fusion genes that contain protein kinase genes

##4. Results

* 337 gene fusions were detected in 80 ALM patient samples of the 112 samples analyzed. The detected fusions have not previously been found in databases such as FusionGDB and COSMIC.
* We found gene fusions with a higher frequency (8-15%) in the samples such as RP11-680G10.1,GSE1, SAMD5, and SASH1.
* We found protein kinase fusion genes with a coverage of 16.96% of tumor sampples.

##5. Conclusion

* Gene fusions were detected with a frequency greater than 10%, having a potential to be used as markers in ALM, however, it is necessary to carry out a greater number of studies to test whether they are ideal markers.
* Kinase gene fusions were found at a low frequency, but the study of their function could be relevant to understanding the development of ALM.

##6. Dependencies

  Fuseq v1.1.4
  STAR-fusion v1.10.1
  FusionInspector v2.6.0
  Program language: R v3.6.2
  Libraries:
    tidyverse v1.3.1
    dplyr v1.0.7
    ggplot2 v3.3.5
    chimeraviz v1.12.3

##7. License

[![License: CC BY-NC 4.0](https://licensebuttons.net/l/by-nc/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc/4.0/)

##8. Author

Ana Cecilia Romani Vasquez
