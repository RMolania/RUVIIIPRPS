#' is used to compute PCA of the gene expression (assay) of a SummarizedExperiment class object.
#'
#' @param se A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay_names Optional string or list of strings for selection of the names
#' of the assays of the SummarizedExperiment class object to compute the PCA.
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#' @param scale Either a logical value or a numeric-alike vector of length equal
#' to the number of columns of the gene expression (assay) of a SummarizedExperiment class object.
#' It is a generic function to scale the columns of a numeric matrix, default is set to 'FALSE'.
#' @param center Either a logical value or a numeric-alike vector of length equal
#' to the number of columns of the gene expression (assay) of a SummarizedExperiment class object.
#' It is a generic function to center the columns of a numeric matrix, default is set to 'FALSE'.
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
    apply.log=FALSE,
    scale=FALSE,
    center=TRUE
){
    ### Check se and assay names
    if (!class(se)[1] == 'SummarizedExperiment') {
        stop('Please provide a summarized experiment object.\n')
    } else if((!is.null(assay_names))&& (any(assay_names %in% names(assays(se)))=='FALSE')){
        stop('The selected assay is/are not in the assay names of the SummarizedExperiment class object.\n')
    }
    ## Assays
    if (!is.null(assay_names)){
        normalization=as.factor(unlist(assay_names))
    }else{
        normalizations=as.factor(names(assays(se)))
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
                    center = center,
                    scale = scale
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
