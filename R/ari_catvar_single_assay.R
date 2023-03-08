#' is used to compute the adjusted rank index (ARI) from the first PC of a single assay
#' given a categorical variable
#'
#'
#' @param select_pca PCs of the dataset that will be used
#' @param cat_var is a categorical variable such as sample types or batches
#' @param nPCs is the number of PCs used to measure the distance
#'
#' @return a single value, the ari computed
#' @importFrom stats dist
#' @importFrom mclust mclustBIC Mclust adjustedRandIndex
#' @export


ari_catvar_single_assay <- function(
        select_pca,
        cat_var,
        nPCs=3
){
    BIC <- mclustBIC(data = select_pca)
    mod <- Mclust(data = select_pca, x = BIC)
    ari=adjustedRandIndex(
        mod$classification,
        cat_var)
    return(ari)
}
