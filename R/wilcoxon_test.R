#' is used to compute the Wilcoxon Rank Sum and Signed Rank Test to data on a categorical variable
#'
#'
#' @param expr.data is the gene expression matrix genes by samples
#' @param apply.log Indicates whether to apply a log-transformation to the data
#' @param variable is a categorical variable such as sample types or batches
#' @param n.cores is the number of cpus used for mclapply parallelization
#'
#' @return dataframe containing the genes, the p-values before and after BH correction
#' @importFrom stats wilcox.test p.adjust
#' @importFrom parallel mclapply
#' @import ggplot2
#' @export


wilcoxon_test <- function(
  expr.data,
  apply.log=FALSE,
  variable,
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
    function(x) wilcox.test(expr.data[x ,] ~ variable)[[3]], mc.cores = n.cores)
    results <- data.frame(
        genes = row.names(expr.data),
        pvalue = unlist(pval),
        ad.pvalue = p.adjust(p = unlist(pval), method = 'BH')
  )
  return(results)
}
