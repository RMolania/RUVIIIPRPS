
#' is used to compute the linear regression of a continuous variable and the first cumulative PCs
#' of the data on a single assay
#'
#'
#' @param pca PCs of the dataset that will be used in the plot
#' @param regression_var The continous variable that will be computed to the PCA of the data
#' (i.e. library size)
#' @param nb_pca_comp The number of components of the PCA used to compute the regression
#'
#' @return vector containing the computed regression
#' @importFrom stats lm
#' @export


regression_pc_contvar_single_assay<-function(
        pca,
        regression_var,
        nb_pca_comp=10
){
    pcs <- pca$sing.val$u
    rSquared <- sapply(
        1:nb_pca_comp,
        function(y) {
            lm.ls <- summary(lm(
                regression_var ~ pcs[, 1:y])
            )$r.squared
        })
    return(rSquared)
}
