#' is used to compute PCA of the gene expression (assay) of a SummarizedExperiment class object.
#'
#' @param se A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay_names Optional string or list of strings for selection of the names
#' of the assays of the SummarizedExperiment class object to compute the PCA.
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#'
#' @return pca.all List of the PCA components computed
#' @importFrom SummarizedExperiment assays
#' @importFrom BiocSingular bsparam
#' @import ggplot2
#' @export
#'

compute_pca=function(
    se,
    assay_names=NULL,
    apply.log=FALSE
){
    if (!is.null(assay_names)){
        normalizations=assay_names
    }else{
        normalizations=names(
        SummarizedExperiment::assays(se))
    }
    pca.all <- lapply(
        normalizations,
        function(x){
          compute_pca_single_assay <- function(y) {
                dat=as.matrix(SummarizedExperiment::assay(y, x))
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
          return(compute_pca_single_assay(se))
        })
    names(pca.all) <- normalizations
    return(pca.all)
}
