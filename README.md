# NatImmunol-Aug2018
Code and data to reproduce the analysis/figures from this manuscript:
[Human germinal center transcriptional programs are de-synchronized in B cell lymphoma, Nature Immunology](https://www.nature.com/articles/s41590-018-0181-4)

For greater simplicity, quality control was already performed on the input data for these R scripts.
All the data files needed are in the folder `/data` , expect for the two scRNA-seq datasets.

The scRNA-seq data are available on this [GEO link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115795).

You can also download them directly:

-To download `GSM3190075_humanBcells_matrixGeneCell.csv`, please click on this [link](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3190075&format=file&file=GSM3190075%5FhumanBcells%5FmatrixGeneCell%2Ecsv%2Egz)

-To download `GSM3190076_humanFL5_matrixGeneCell.csv`, please click on this [link](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3190076&format=file&file=GSM3190076%5FhumanFL5%5FmatrixGeneCell%2Ecsv%2Egz)


The R scripts are ordered by Figures and by type of data (sc-qPCR and scRNA-seq). Most of the scripts can be run independently.

Please use the [Issues Tracker](https://github.com/MilpiedLab/NatImmunol-Aug2018/issues) for questions/issues.

cervera [at] ciml.univ-mrs.fr 
milpied [at] ciml.univ-mrs.fr
