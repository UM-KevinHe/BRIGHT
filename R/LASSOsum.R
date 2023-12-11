#' Fitting LASSOsum PRS (equivalent to BRIGHTs with eta=0)
#'
#' Fitting LASSOsum PRS. Quality control is performed first on the summary statistics (please see QC section). The hyperparameters are tuned through BRCp (equivalent to AIC).
#'
#' @param ss a data frame of the summary statistics, including "GWAS" or "Corr" corresponding to the GWAS or correlation summary statistics, the input summary type must match with ind below.
#' @param ind a character scaler of the input summary type, it must be either "GWAS" or "Corr" corresponding to the GWAS or correlation summary statistics. It must match with the summary statistics in ss.
#' @param ref a character scaler indicating the LD reference. It takes values as either a directory to user provided Plink 1 binary files (.bim, .bed, .fam) or the 1000 genome project reference provided by the package: "EUR.hg19", "EUR.hg38", "AFR.hg19", "AFR.hg38", "EAS.hg19", "EAS.hg38", "AMR.hg19", "AMR.hg38", "SAS.hg19" or "SAS.hg38". The abbreviations before dot are populations (more detials can be found here: https://useast.ensembl.org/Help/Faq?id=532); The abbreviations after dot are the genetic build (more details can be found: https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19).
#' @param LDblocks a character scaler indicating the LD block structure to be used from Berisa et al for the population. It takes values from "EUR.hg19", "EUR.hg38", "AFR.hg19", "AFR.hg38", "EAS.hg19", "EAS.hg38", "AMR.hg19", "AMR.hg38", "SAS.hg19" or "SAS.hg38". The abbreviations before dot are populations (more detials can be found here: https://useast.ensembl.org/Help/Faq?id=532); The abbreviations after dot are the genetic build (more details can be found: https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19). The genetic build must match with those provided in ref.
#' @param N a integer scaler indicating the sample sizes of the summary statistics input in ss.
#' @return A list of LASSOsum models corresponding to different hyperparameters. The BRcp are also output for hyperparameter fine tuning.

LASSOsum=function(ss, ind, LDblocks, KG, N ,ref, lab = "LASSOsum"){
  if(names(ss)[6] != "Coef"){
    ss = QC(ss, ind, KG, lab = lab)
    #writeLines(strwrap("Coefficients not identified, starting fitting LASSOsum on prior data"))
    data = list()
    data[["Tss"]] = ss
    data[["Tind"]] = ind
    data[["TLDblocks"]] = LDblocks
    data[["Tref"]] = ref
    data[["LD_dir"]] = paste(ref,".bed", sep = "")
    data[["TN"]] = N
    Prior_Lassosum = BRIGHTs(data, m = NA, group = NA, penalty = 1, tau = 0, eta = 0, lambda = NA, nlambda = 50, lambda.min = 0.0001, alpha = 1, gamma = 0, eps = 0.001, max_iter = 1000000, dfmax = 5000, gmax = 5000, user = T)
    if(is.na(N)){
      writeLines(strwrap("sample size of prior data must be provided for the AIC calculation, otherwise please directly input Coefficients from joint models"))
    }else{
      AIC = 2 * (Prior_Lassosum[[as.character("0")]]$dev2*N + colSums(as.matrix(Prior_Lassosum[[as.character(0)]]$Beta) != 0))
    }
    rst = list()
    rst[["AIC"]] = AIC
    rst[["model"]] = Prior_Lassosum
    return(rst)
  }else{
    stop('"Coef" found but "IProd" or "GWAS required for LASSOsum fitting, aborting"')
  }
}
