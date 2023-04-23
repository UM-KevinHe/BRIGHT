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
out <- BRIGHTs(data = dat,type.trait="quantitative",penalty="LASSO")
```

## BRIGHTs with binary traits
For binary traits, in addition to the GWAS summary statistics or marginal genotype-trait inner product, $\frac{\boldsymbol X^\top\boldsymbol y}{n}$, BRIGHTs requires an estimate based on logistic LASSO regression, $\boldsymbol{\hat b}$, and the marginal genotype-predicted traits inner product, $\frac{\boldsymbol X^\top expit(\boldsymbol X \boldsymbol{\hat b})}{n}$, from the target minority population. From the prior majority populations coefficients estimated from joint models (e.g. logistic LASSO regression) is required for model fitting. We note that more than 1 prior majority data can be incorporated in the BRIGHTs model. This procedure requires additional and quite stringent summary statistics from both target and prior data, in genetics studies its quite common to treat binary outcome as continuous and perform continuous models on the data; therefore, in the case where the above additonal summary statistics are not available, the BRIGHTS with quantitative traits procedure can also be used to analyze the binary data.

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
out <- BRIGHTs(data = dat,type.trait="quantitative",penalty="LASSO")
```

#### Parallel processing with the `parallel` package
Note that parallel processing is done by `LDblocks`. 
```r
library(parallel)
cl <- makeCluster(2, type="FORK") # Parallel over 2 nodes
out <- lassosum.pipeline(cor=cor, chr=ss$Chr, pos=ss$Position, 
                         A1=ss$A1, A2=ss$A2, # A2 is not required but advised
                         ref.bfile=ref.bfile, test.bfile=test.bfile, 
                         LDblocks = LDblocks, cluster=cl)
```
#### Including covariates in validation
It is possible to include covariates in validation or splitvalidation (though not in pseudovalidation). Simply pass the covariate matrix as an argument to `validate` or `splitvalidate`. 
```r 
v <- validate(out, covar=covar)
# covar <- rnorm(nrow.bfile(out$test.bfile)) # If you need a dummy for testing
```
Since v0.4.2, the `covar` argument in `validate` and `splitvalidate` can also take a `data.frame` with the first 2 columns headed by FID and IID, and the other columns being covariates (any headers). It can also be a file name for such a data.frame.

#### Apply validated betas to new data 
To apply the best lassosum predictor (indexed by `s` and `lambda`) to a new dataset, first subset the `lassosum.pipeline` object. Then `validate` again: 
```r 
out2 <- subset(out, s=v$best.s, lambda=v$best.lambda)
v2 <- validate(out2, covar=covar, test.bfile="Some_new_bfile")
```

#### lassosum in subsamples  
To use subsamples of either the test dataset or the reference panel, use the `keep.test` and `keep.ref` options respectively in `lassosum.pipeline`. Type `?lassosum.pipeline` to see what format these options can take. To use only a subset of SNPs for inference, remove the unnecessary SNPs from the summary statistics file. 

#### Estimation by chromosome 
In large datasets, it is typical for the PLINK data to be organized by chromosomes. For example, say we have `chr1.bed` and `chr2.bed`. It is much faster to carry out lassosum by chromosome, and then combine them together. An example: 
```r 
out1 <- lassosum.pipeline(..., test.bfile="chr1")
out2 <- lassosum.pipeline(..., test.bfile="chr2")
out <- merge(out1, out2)
v <- validate(out)
```
#### For datasets with large sample size in `ref.bfile`
To guard against taking too long to run, `lassosum.pipeline` will complain if the size of the reference panel is too large. Currently "too large" is taken to mean greater than 20000. `lassosum.pipeline` will suggest that you take a sample of, say, 5000 as the reference panel using the `sample=` option. Note that if you run `lassosum.pipeline` across different chromosomes, you need to keep the sample the same. You can do so by running `set.seed()` before `lassosum.pipeline`. Alternatively, specify the exact sample with the `keep.ref=` option. 

#### For datasets with large sample size in `test.bfile`
In general, it is a good idea to anticipate the test sample you will need in the `lassosum.pipeline` stage rather than specify them at the `validate` stage, using the `keep.test=` option. This is especially the case if you need to specify the `pheno` or the `covar` option in `validate`. The reason is that if a different test sample is used (or if a different `test.bfile` is specified) in `validate`, then the polygenic score is re-calculated from the betas, this can be very time consuming. In a future version, I may allow `lassosum.pipeline` to take `pheno` and `covar` as argument so as to identify the exact test sample to circumvent this issue. 

#### Extracting the lassosum betas 
The exact SNPs used in calculating the PGS are given in the `out$sumstats` data.frame, where `out` is the output from `lassosum.pipeline`. The order of the SNPs is the same as that of `v$best.beta` where `v` is the output from `validate`/`splitvalidate`/`psuedovalidate`. 

### Support
Further documentation for various functions can be obtained by running `help` from `R`, e.g. 
```r
help(lassosum.pipeline)
help(validate)
help(pgs) # A low-level command for calculating PGS from betas 
```
If there are any further questions or problems with running or installing `lassosum`, please do email me at <timmak@yahoo.com>. 
