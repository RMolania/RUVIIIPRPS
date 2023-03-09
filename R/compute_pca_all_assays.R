#' is used to compute PCA of the data on all assays
#'
#'
#' @param sce Dataset that will be used to compute the PCA
#' @param apply.log Indicates whether to apply a log-transformation to the data
#'
#'
#' @return pca.all PCA components computed on the provided data
#' @importFrom SummarizedExperiment assays
#' @importFrom BiocSingular bsparam
#' @import ggplot2
#' @export
#'
compute_pca_all_assays=function(
    sce,
    #assay_names,
    apply.log=FALSE
){
    normalizations=names(
        SummarizedExperiment::assays(sce)
    )
    pca.all <- lapply(
        normalizations,
        function(x){
          compute_pca_single_assay <- function(y) {
                dat=as.matrix(SummarizedExperiment::assay(sce, x))
                if (apply.log==FALSE)
                    dat <- dat
                else
                    dat <- log2(dat + 1)
                svd <- base::svd(scale(
                    x = t(dat),
                    center = TRUE,
                    scale = FALSE
                ))
                percent <- svd$d ^ 2 / sum(svd$d ^ 2) * 100
                percent <-
                    sapply(seq_along(percent),
                           function(i) {
                               round(percent[i], 1)
                           })
                pca=list(
                    sing.val = svd,
                    variation = percent)
                return(pca)
          }
          return(compute_pca_single_assay)
        })
    names(pca.all) <- normalizations
    return(pca.all)
}
