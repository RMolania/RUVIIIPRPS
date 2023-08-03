
#' is used to normalise the gene expression (assay) of a SummarizedExperiment class object
#' using RUVIII-PRPS method for multiple k values in one go.
#'
#'
#' @param Y A m by n matrix of the Raw gene expression matrix where m is the number of samples and n is the
#' number of features of a SummarizedExperiment variable to be normalised.
#' @param M Replicate matrix.
#' @param ctl Logical vector of length n of the negative control genes.
#' @param k A single value or a vector of values containing a single k or a range of k - the number of unwanted factors - to be tested.
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

normalise <- function(
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

    ## Run RUVIII_PRPS for the first k value provided (maybe the only one)
    ruv.adj<- ruv_III_prps(
                Y,
                M,
                ctl,
                k = k[1],
                eta = eta,
                include.intercept = include.intercept,
                average = average,
                fullalpha = fullalpha,
                return.info = return.info,
                inputcheck = inputcheck)
    ruv.adj=list(ruv.adj)
    names(ruv.adj) <- paste0('RUV_K', k[1])
    return(ruv.adj)

    ## if there are multiple k values provided
    if (length(k)>1) {
    ruv.adj.others.k <- lapply(
        k[2:length(k)],
        function(x){
            ruv.adj.k<- ruv_III_prps(
                Y,
                M,
                ctl,
                k = x,
                eta = eta,
                include.intercept = include.intercept,
                average = average,
                fullalpha = fullalpha,
                return.info = return.info,
                inputcheck = inputcheck
            )
        })
    names(ruv.adj.others.k) <- paste0('RUV_K', k)
    ruv.adj.allk=append(ruv.adj.others.k, ruv.adj, after = 0)
    ruv.adj=ruv.adj.allk
    }
    ## Return a list containing the adjusted dataset(s) for single k or multiple k values
    return(ruv.adj)
}

