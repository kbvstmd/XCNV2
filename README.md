1. Introduction
X-CNV model requires two R packages, data.table and xgboost, and Bedtools v2.26.0. If the R packages and bedtools can not be installed automatically, users can install them by yourself. The executable file of bedtools should be placed in ./tools/.

2. Installation
sh Install.sh

3. Predicting CNV pathogenicity using X-CNV model
Usage:
./bin/XCNV CNV.bed

