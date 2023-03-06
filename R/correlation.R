
#' is used to compute the correlation (Spearman or Pearson) between the data and a continuous variable (i.e. library size)
#'
#' @param expr.data is the gene expression matrix genes by samples
#' @param apply.log Indicates whether to apply a log-transformation to the data
#' @param variable is a categorical variable such as library size or purity
#' @param method to select which method to use either Spearman or Pearson correlation
#' @param n.cores is the number of cpus used for mclapply parallelization
#'
#' @return dataframe containing the genes, the rho, the pvalues before and after BH correction
#' @importFrom stats cor.test p.adjust
#' @importFrom parallel mclapply
#' @import ggplot2
#' @export

correlation<- function(
        expr.data,
        apply.log,
        variable,
        method,
        n.cores
){
    if(apply.log==FALSE){
        expr.data <- expr.data
    }
    else{
        expr.data <- log2(expr.data + 1)
    }
    rho <- mclapply(
        1:nrow(expr.data),
        function(x){
            round(cor.test(
                x = expr.data[x, ],
                y = variable,
                method = method)[[4]], 6)},
        mc.cores = n.cores
    )
    pval <- mclapply(
        1:nrow(expr.data),
        function(x){
            cor.test(
                x = expr.data[x, ],
                y = variable,
                method = method)[[3]]},
        mc.cores = n.cores)

    results <- data.frame(
        genes = row.names(expr.data),
        rho = unlist(rho),
        pvalue = unlist(pval),
        adj.pvalue = p.adjust(unlist(pval), 'BH')
    )
    return(results)
}
