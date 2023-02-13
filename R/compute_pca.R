#' is used to compute PCA of the data
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
compute_pca=function(
    sce,
    apply.log=FALSE
){
    normalizations=names(
        SummarizedExperiment::assays(sce)
    )
    pca.all <- lapply(
        normalizations,
        function(x){
            pca(sce= as.matrix(SummarizedExperiment::assay(sce, x)),apply.log = apply.log)
        })
    names(pca.all) <- normalizations
    return(pca.all)
}


#=================== PCA =================
# Principal component analysis using singular value decomposition (SVD)
## sce: Dataset that will be used to compute the PCA
## is.log: Indicates whether to apply a log-transformation to the data
pca <- function(sce, apply.log=FALSE) {
    if (apply.log==FALSE)
        sce <- sce
    else
        sce <- log2(sce + 1)
    svd <- base::svd(scale(
        x = t(sce),
        center = TRUE,
        scale = FALSE
    ))
    percent <- svd$d ^ 2 / sum(svd$d ^ 2) * 100
    percent <-
        sapply(seq_along(percent),
               function(i) {
                   round(percent[i], 1)
               })
    return(list(
        sing.val = svd,
        variation = percent))
}

