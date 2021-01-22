1. Introduction
X-CNV is a tool to predict CNV pathogenicity using an XGBoost classifier, which calculates a meta-voting prediction (MVP) score to quantitatively evaluate disease-causing probability. It consists of the most comprehensive CNV data and annotations by integrating various publicly available genetic variant repositories. The features covering the genomics, genome region, variation types, and population genetics properties are taken into account to boost the prediction power. More importantly, a meta-voting prediction (MVP) score is proposed to measure the CNV pathogenic effect quantitatively, which can be used to determine the CNV pathogenicity. 

2. Requirements
The local version X-CNV requires two R packages, data.table and xgboost, and Bedtools v2.26.0. If the R packages and bedtools cannot be installed automatically, users can install them manually. The executable file of bedtools should be placed in ./tools/.

3. Installation
git clone xx
sh Install.sh

4. Usage and example
Usage:
./bin/XCNV prefix.bed
The output filename: prefix.output.csv

Example:
./bin/XCNV ./example_data/1.beds
The results can be seen in the 1.output.csv

5. Input & output
Input file format (The columns are separated by TAB key and the header is not required): 
2	2222999	3000222	gain

Column 1: The chromosome (no “chr”)
Column 2: Start
Column 3: End
Column 4: CNV type (gain or loss)

The output file has 35 columns and is provided as Comma-Separated Values (CSV) format. 

Columns	Description	Category
Chr	Chromosome	Input
Start	Start position	Input
End	End position	Input
CNV type	CNV type (gain or loss)	Input
FATHMM score	FATHMM Dnase score for the CNV region	Coding
LR score	LR score for the CNV region	Coding
LRT score	LRT Dnase score for the CNV region	Coding
MutationAssessor score	MutationAssessor score for the CNV region	Coding
MutationTaster score	MutationTaster score for the CNV region	Coding
Polyphen2-HDIV score	Polyphen2_HDIV score for the CNV region	Coding
Polyphen2-HVAR score	Polyphen2_HVAR score for the CNV region	Coding
RadialSVM score	RadialSVM Dnase score for the CNV region	Coding
SIFT score	SIFT Dnase score for the CNV region	Coding
VEST3 score	VEST3 score for the CNV region	Coding
pLI	Probability of being loss-of-function intolerant	Coding
Episcore	A computational method to predict haploinsufficiency leveraging epigenomic data from a broad range of tissue and cell types by machine learning methods.	Coding
GHIS	An integrative approach to predicting haploinsufficient genes	Coding
CADD score	Average CADD score for the CNV region	Genome-wide
GERP	GERP++_RS Dnase score for the CNV region	Genome-wide
phyloP100way	phyloP100way_vertebrate score for the CNV region	Genome-wide
phyloP46way	phyloP46way_placental score for the CNV region	Genome-wide
SiPhy29way	SiPhy_29way_logOdds score for the CNV region	Genome-wide
cdts-1st	The coverage ratio between  CDTS percentile < 1% and the CNV region	Noncoding
cdts-5th	The coverage ratio between  CDTS percentile < 5% and the CNV region	Noncoding
pELS	The coverage of proximal enhancer-like sequence (pELS) within the CNV region	Noncoding
CTCF-bound	The coverage of CTCF-bound sequence within the CNV region	Noncoding
PLS	The coverage of promoter-like sequence within the CNV region	Noncoding
dELS	The coverage of distal enhancer-like sequence within the CNV region	Noncoding
CTCF-only	The coverage of CTCF-only sequence within the CNV region	Noncoding
DNase-H3K4me3	The coverage of DNase-H3K4me3 sequence within the CNV region	Noncoding
gain-PAF	Population allele frequency for duplication	Universal
Length	CNV length	Universal
loss-PAF	Population allele frequency for deletion	Universal
Type.1	CNV type (gain or loss code as 1 or 0)	Universal
MVP_score	The MVP score	Output



