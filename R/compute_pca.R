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


    #=================== PCA =================
    # Principal component analysis using singular value decomposition (SVD)
    ## data: gene expression matrix genes by samples
    ## is.log:  logical factor indicating whether the data is log transformed or not
    .pca <- function(data, is.log) {
        if (is.log)
            data <- data
        else
            data <- log2(data + 1)
        svd <- base::svd(scale(
            x = t(data),
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

    # Fast principal component analysis
    ## data: gene expression matrix genes by samples
    ## nPcs: integer scalar specifying the number of singular values to return
    ## is.log:  logical factor indicating whether the data is log transformed or not
    fast.pca <- function(data, nPcs, is.log) {
        if(is.log){
            data <- data
        }else{
            data <- log2(data + 1)
        }
        svd <- BiocSingular::runSVD(
            x = t(data),
            k = nPcs,
            BSPARAM = BiocSingular::bsparam(),
            center = TRUE,
            scale = FALSE
        )
        percent <- svd$d^2/sum(svd$d^2)*100
        percent <-
            sapply(
                seq_along(percent),
                function(i) {round(percent[i], 1)})
        return(list(
            sing.val = svd,
            variation = percent))
    }
}
