#' is used to compute PCA of a single assay
#'
#'
#' @param sce Dataset that will be used to compute the PCA
#' @param apply.log Indicates whether to apply a log-transformation to the data
#'
#'
#' @return pca List containing the svd and the total variation explained by the PCA components computed on the single data
#' @export

compute_pca_single_assay <- function(sce, apply.log=FALSE) {
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
    pca=list(
        sing.val = svd,
        variation = percent)
    return(pca)
}

