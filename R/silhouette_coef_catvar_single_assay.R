#' is used to compute the Silhouette coefficient from the PC of a single assay
#' given a categorical variable
#'
#'
#' @param pca PCs of the dataset that will be used
#' @param cat_var is a categorical variable such as sample types or batches
#' @param nPCs is the number of PCs used to measure the distance
#'
#' @return a single value, the silhouette coefficient computed
#' @importFrom stats dist
#' @importFrom cluster silhouette
#' @export


silhouette_coef_catvar_single_assay <- function(
        pca,
        cat_var,
        nPCs=3
){
    d.matrix <- as.matrix(dist(pca[, seq_len(nPCs)]))
    avg=summary(silhouette(
        as.numeric(as.factor(cat_var)),
        d.matrix))$avg.width
    return(avg)
}
