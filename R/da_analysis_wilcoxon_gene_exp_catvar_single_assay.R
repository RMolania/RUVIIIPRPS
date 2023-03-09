#' is used to compute the Wilcoxon Rank Sum and Signed Rank Test on the gene expression
#' of a single assay given a categorical of two groups
#'
#'
#' @param expr.data is the gene expression matrix genes by samples
#' @param apply.log Indicates whether to apply a log-transformation to the data
#' @param cat_var is a categorical cat_var such as sample types or batches
#' @param n.cores is the number of cpus used for mclapply parallelization
#'
#' @return dataframe containing the genes, the p-values before and after BH correction
#' @importFrom stats wilcox.test p.adjust
#' @importFrom parallel mclapply
#' @import ggplot2
#' @export


da_analysis_wilcoxon_gene_exp_catvar_single_assay <- function(
  expr.data,
  apply.log=FALSE,
  cat_var,
  n.cores
){
  if(apply.log==FALSE){
    expr.data <- expr.data
  }
  else{
    expr.data <- log2(expr.data + 1)
  }
  pval <- mclapply(
    row.names(expr.data),
    function(x) wilcox.test(expr.data[x ,] ~ cat_var)[[3]], mc.cores = n.cores)
    results <- data.frame(
        genes = row.names(expr.data),
        pvalue = unlist(pval),
        ad.pvalue = p.adjust(p = unlist(pval), method = 'BH')
  )
  return(results)
}
