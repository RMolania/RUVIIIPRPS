#' is used to compute the Silhouette coefficient given a categorical variable
#'
#'
#' @param pca PCs of the dataset that will be used
#' @param variable is a categorical variable such as sample types or batches
#' @param nPCs is the number of PCs used to measure the distance
#'
#' @return a single value, the silhouette coefficient computed
#' @importFrom stats dist
#' @importFrom cluster silhouette
#' @export


silhouette_coeff <- function(
        pca,
        variable,
        nPCs
){
    d.matrix <- as.matrix(dist(pca[, seq_len(nPCs)]))
    avg=summary(silhouette(
        as.numeric(as.factor(variable)),
        d.matrix))$avg.width
    return(avg)
}
