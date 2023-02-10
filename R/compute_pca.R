#' @rdname compute_pca
#' @title compute_pca
#'
#' is used to compute PCA of the data
#'
#'
#' @param sce Dataset that will be used to compute the PCA
#'
#'
#'
#' @return pca.all PCA components computed on the provided data
#' @export
#'
compute_pca=function(
        ##compute_pca(skcm.se)
    sce,
    apply.log=FALSE
){
    normalizations=names(
        SummarizedExperiment::assays(sce)
    )
    pca.all <- lapply(
        normalizations,
        function(x){
            .pca(
                data = as.matrix(
                    SummarizedExperiment::assay(sce, x)
                ),
                is.log = !apply.log)
        })
    names(pca.all) <- normalizations
    return(pca.all)
}
