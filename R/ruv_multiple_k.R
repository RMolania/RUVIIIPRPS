
#' is used to normalise the gene expression (assay) of a SummarizedExperiment class object
#' using RUVIII-PRPS method for multiple k values in one go.
#'
#'
#' @param Y A m by n matrix of the Raw gene expression matrix where m is the number of samples and n is the
#' number of features of a SummarizedExperiment variable to be normalised.
#' @param M Replicate matrix.
#' @param ctl Logical vector of length n of the negative control genes.
#' @param k A vector containing the range of k - the number of unwanted factors - to be tested.
#' @param eta A matrix with n columns, gene-wise (as opposed to sample-wise) covariates. By default is set to NULL.
#' @param include.intercept Add an intercept term to eta if it does not include one already. By default is set to TRUE.
#' @param average Average replicates after adjustment. By default is set to FALSE.
#' @param fullalpha Can be included to speed up execution. By default is set to NULL.
#' @param return.info If FALSE, only the adjusted data matrix is returned. If TRUE, additional information
#' is returned. By default is set to FALSE.
#' @param inputcheck Perform a basic sanity check on the inputs, and issue a warning if there is a problem.
#' By default is set to TRUE.
#'
#' @return list List containing the corrected gene expression matrix and the M, alpha and W for each k values tested.

#' @importFrom ruv RUV1 residop
#' @export

ruv_multiple_k <- function(
        Y,
        M,
        ctl,
        k = NULL,
        eta = NULL,
        include.intercept = TRUE,
        average = FALSE,
        fullalpha = NULL,
        return.info = FALSE,
        inputcheck = TRUE
){

    ruv.adj.kother <- lapply(
        k,
        function(x){
            ruv.adj.at.other <- ruv_III_prps(
                Y,
                M,
                ctl,
                k = k,
                eta = eta,
                include.intercept = include.intercept,
                average = average,
                fullalpha = fullalpha,
                return.info = return.info,
                inputcheck = inputcheck
            )
        })
    names(ruv.adj.kother) <- paste0('RUV_K', k)
    return(ruv.adj.kother)
}

