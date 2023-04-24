BRIGHT
=======================

# Description

`BRIGHT` is a group of methods for using individual-level data, published genotype-trait summary statistics (GWAS), or combined data types from different ethnic populations to improve the recalibration, discrimination, and prediction accuracy on the target minority cohort. We implemented group LASSO, Elastic Net, MCP and SCAD penalties for marker fine mapping.

# Reference
Li, Q., Patrick, M. T., Zhang, H., Khunsriraksakul, C., Stuart, P. E., Gudjonsson, J. E., ... & He, K. (2022). Bregman Divergence-Based Data Integration with Application to Polygenic Risk Score (PRS) Heterogeneity Adjustment. arXiv preprint arXiv:2210.06025 (https://arxiv.org/abs/2210.06025)

# Installation

`BRIGHT` requires the following R packages: `Rcpp`, `Matrix`, `mvtnorm`, `BEDMatrix`, `grpreg`. Install them by: 

```r
install.packages(c("Rcpp", "Matrix", "mvtnorm", "BEDMatrix", "grpreg"), dependencies=TRUE)
```
For Windows and Mac users, it would be easiest to download the following binaries ([Windows](To be determined), [Mac](To be determined)) and install using: 
```r
install.packages("/path/to/downloaded_binary_file", repos=NULL)
```

If you are on Linux or you would like to compile from source, you can download the source codes [BRIGHT_0.0.1.tar.gz](To be determined). Mac users should refer to [this page](https://cran.r-project.org/bin/macosx/tools/) for the various dependencies required. Install then via: 
```r
install.packages("/path/to/downloaded_source.tar.gz", repos=NULL, type="source")
```

If you have `devtools`, you can also type: 
```r
install_github("To be determined")
```
or
```r
install_github("To be determined")
```
for the latest development version. Or you can clone the latest development version here and install yourself using `devtools`. 

# Warning!

Most functions in `BRIGHT` impute missing genotypes in PLINK bfiles with a homozygous A2 genotype, which is the same as using the `--fill-missing-a2` option in PLINK. It is the user's responsibility to filter out individuals and SNPs with too many missing genotypes beforehand. 

# BRIGHTs tutorial
BRIGHTs group of methods utilize a wide variety of summary-level data from different populations to carry out transfer-learning. We accounted for Linkage Disequilibrium (LD) via a reference panel (1000 genome project as default).
The reference panel is assumed to be in PLINK 1 [format](https://www.cog-genomics.org/plink/1.9/input#bed).
Summary statistics are expected to be loaded into memory as a data.frame/data.table. 

Below we discuss the required data and implementation tutorials separately for quantitative traits and binary traits.

## BRIGHTs with quantitative traits
For quantitative traits, BRIGHTs requires the GWAS summary statistics or marginal genotype-trait inner product, $\frac{\boldsymbol X^\top\boldsymbol y}{n}$, from the target minority population, while from the prior majority populations either GWAS summary statistics, marginal genotype-trait inner product, or coefficients estimated from joint models (e.g. PRS or LASSO regression) can be used for model fitting. We note that more than 1 prior majority data can be incorporated in the BRIGHTs model.

First we read the minority summary statistics and majority summary statistics into R, and provide the `ref` names of the reference panel. If `ref` names are provided as "EUR", "AFR", "EAS", "SAS" ,or "AMR", then the default 1000 genome project reference panels will be used; otherwise `ref` needs to be provided as a directory to the plink1 format files (.bim, .bed, .fam). 


```r
library(BRIGHT)
library(data.table)

### Read target minority GWAS summary statistics file or marginal genotype-trait inner product file###

# Read in target GWAS
Tind="GWAS"
Tss <- fread("Target_GWAS.txt")
head(Tss)

# Alternatively read in target marginal genotype-trait inner product
Tind="IProd"
Tss <- fread("Target_IProd.txt")
head(Tss)

### Read prior majority GWAS summary statistics file, marginal genotype-trait inner product, or joint coefficient estimates, more than 1 prior majority data can be read in###

Pind=c("GWAS","IProd","Coef")
Pss1 <- fread("Prior_GWAS1.txt")
head(Pss1)
Pss2 <- fread("Prior_IProd2.txt")
head(Pss2)
Pss3 <- fread("Prior_Coef3.txt")
head(Pss3)
Pss=list("1"=Pss1,"2"=Pss2,"3"=Pss3) # The order of list Pss need to be matched with Pind

### Specify the PLINK file stub of the reference panel or "EUR", "AFR", "EAS", "SAS" ,or "AMR" ###
ref.bfile <- "refpanel"

### Read LD region file, only required if ref.bfile is provided as PLINK1 format ###
LDblocks <- "AFR.hg19" # This will use LD regions as defined in Berisa and Pickrell (2015) for the African population and the hg19 genome build.
# Other alternatives available. Type ?BRIGHTs for more details. 
```
Reference: [Berisa and Pickrell (2015)](https://academic.oup.com/bioinformatics/article/32/2/283/1743626/Approximately-independent-linkage-disequilibrium)

Then, a preprocessing step is required to remove the SNPs that are not in the reference panel from all data, convert target data into marginal SNPs-trait inner product, convert prior data into joint coefficient estimates, and match the effect alleles between the reference panel and data.

```r
dat <- PreprocessS(Tss = Tss, Tind = Tind, Pss = Pss, Pind = Pind, ref.bfile=ref.bfile, LDblocks=LDblocks)
```

Running BRIGHTs using standard pipeline with LASSO penalty on quantitative traits: 
```r
out <- BRIGHTs(data = dat, type.trait="quantitative", penalty="LASSO")
```

## BRIGHTs with binary traits
For binary traits, in addition to the GWAS summary statistics or marginal genotype-trait inner product, $\frac{\boldsymbol X^\top\boldsymbol y}{n}$, BRIGHTs requires an estimate based on logistic LASSO regression, $\boldsymbol{\hat b}$, and the marginal genotype-predicted traits inner product, $\frac{\boldsymbol X^\top expit(\boldsymbol X \boldsymbol{\hat b})}{n}$, from the target minority population. From the prior majority populations coefficients estimated from joint models (e.g. logistic LASSO regression) is required for model fitting. We note that more than 1 prior majority data can be incorporated in the BRIGHTs model. 

First we read the minority summary statistics and majority summary statistics into R, and provide the `ref` names of the reference panel. If `ref` names are provided as "EUR", "AFR", "EAS", "SAS" ,or "AMR", then the default 1000 genome project reference panels will be used; otherwise `ref` needs to be provided as a directory to the plink1 format files (.bim, .bed, .fam). 


```r
library(BRIGHT)
library(data.table)

### Read target minority GWAS summary statistics file or marginal genotype-trait inner product file###

# Read in target GWAS
Tind=c("GWAS","LASSO","IProdPred")
Tss1 <- fread("Target_GWAS.txt")
head(Tss1)

# Alternatively read in target marginal genotype-trait inner product
Tind=c("IProd","LASSO","IProdPred")
Tss <- fread("Target_IProd.txt")
head(Tss)

### Read target minority LASSO estimates file and marginal genotype-predicted outcome inner product file###
bhat <- fread("Target_LASSO.txt")
head(bhat)
rhat <- fread("Target_IProdPred.txt")
head(bhat)

Tss <- list("1" = Tss1, "2" = bhat, "3" = rhat) # The order of list Tss need to be matched with Tind

### Read prior majority GWAS summary statistics file, marginal genotype-trait inner product, or joint coefficient estimates, more than 1 prior majority data can be read in###

Pind=c("Coef", "Coef", "Coef")
Pss1 <- fread("Prior_Coef1.txt")
head(Pss1)
Pss2 <- fread("Prior_Coef2.txt")
head(Pss2)
Pss3 <- fread("Prior_Coef3.txt")
head(Pss3)
Pss=list("1"=Pss1,"2"=Pss2,"3"=Pss3) # The order of list Pss need to be matched with Pind

### Specify the PLINK file stub of the reference panel or "EUR", "AFR", "EAS", "SAS" ,or "AMR" ###
ref.bfile <- "refpanel"

### Read LD region file, only required if ref.bfile is provided as PLINK1 format ###
LDblocks <- "AFR.hg19" # This will use LD regions as defined in Berisa and Pickrell (2015) for the African population and the hg19 genome build.
# Other alternatives available. Type ?BRIGHTs for more details. 
```

Then, a preprocessing step is required to remove the SNPs that are not in the reference panel from all data, convert target data into marginal SNPs-trait inner product, convert prior data into joint coefficient estimates, and match the effect alleles between the reference panel and data.

```r
dat <- PreprocessS(Tss = Tss, Tind = Tind, Pss = Pss, Pind = Pind, ref.bfile=ref.bfile, LDblocks=LDblocks)
```

Running BRIGHTs using standard pipeline with LASSO penalty on quantitative traits: 
```r
out <- BRIGHTs(data = dat, type.trait="binary", penalty="LASSO")
```
This procedure requires additional and quite stringent summary statistics from both target and prior data, in genetics studies its quite common to treat binary outcome as continuous and perform continuous models on the data; therefore, in the case where the above additonal summary statistics are not available, the BRIGHTS with quantitative traits procedure can also be used to analyze the binary data.

# BRIGHTi tutorial
BRIGHTi group of methods utilize individual-level data from target minority populations and a wide variety of summary-level data from prior majority population to carry out transfer-learning.
Summary statistics are expected to be loaded into memory as a data.frame/data.table for the prior majority population. 

Below we discuss the required data and implementation tutorials.

BRIGHTi requires the individual-level genotype and phenotype from the target minority population, while from the prior majority populations either GWAS summary statistics, marginal genotype-trait inner product, or coefficients estimated from joint models (e.g. PRS or LASSO regression) can be used for model fitting. We note that more than 1 prior majority data can be incorporated in the BRIGHTi model.

First we read the minority genotype data from plink1 files and phenotype data from text files and majority summary statistics into R.


```r
library(BRIGHT)
library(data.table)

### Read target minority GWAS summary statistics file or marginal genotype-trait inner product file###

# Read in target individual-level data
Tgeno <- "/path/to/plink"
Tpheno <- fread("Target_phenotype.txt")
head(Tpheno)


### Read prior majority GWAS summary statistics file, marginal genotype-trait inner product, or joint coefficient estimates, more than 1 prior majority data can be read in###

Pind=c("GWAS","IProd","Coef")
Pss1 <- fread("Prior_GWAS1.txt")
head(Pss1)
Pss2 <- fread("Prior_IProd2.txt")
head(Pss2)
Pss3 <- fread("Prior_Coef3.txt")
head(Pss3)
Pss=list("1"=Pss1,"2"=Pss2,"3"=Pss3) # The order of list Pss need to be matched with Pind
```

Then, a preprocessing step is required to remove the SNPs that are not in the target minority genotype files from prior majority data, convert prior data into joint coefficient estimates, and match the effect alleles between the minority genotype and prior data.

```r
dat <- PreprocessI(Tpheno = Tpheno, Tgeno = Tgeno, Pss = Pss, Pind = Pind)
```

Running BRIGHTi using standard pipeline with LASSO penalty on different types of traits including "quantitative", "binary", "count": 
```r
out <- BRIGHTi(data = dat, type.trait="quantitative", penalty="LASSO")
```

# Model validation
## Model validation with individual-level test data
When individual-level test data is available, `BRIGHT` package provide automated validation functions and generates evaluation plots:

```r
# Read in target individual-level data
Testgeno <- "/path/to/test/plink"
Testpheno <- fread("Test_phenotype.txt")
head(Testpheno)

# Perform testing
Val <- Valid.Ind(out, Testpheno, Testgeno)
```

## Model validation with summary-level test data
When summary test data is available, `BRIGHT` package provide automated validation functions for parameter fine-tunning:

```r
# Read in test GWAS
Testind="GWAS"
Testss <- fread("Target_GWAS.txt")
head(Testss)

# Alternatively read in test marginal genotype-trait inner product
Testind="IProd"
Testss <- fread("Target_IProd.txt")
head(Testss)

# Perform testing
Val <- Valid.Sum(Testss, Testind, Testpheno, Testgeno)
```

We note that summary level test data is only supported for quantitative traits

