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
            single_pca(sce= as.matrix(SummarizedExperiment::assay(sce, x)),apply.log = apply.log)
        })
    names(pca.all) <- normalizations
    return(pca.all)
}
