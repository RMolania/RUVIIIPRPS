
#' is used to normalise the gene expression (assay) of a SummarizedExperiment class object
#' using RUVIII-PRPS method for multiple k values in one go.
#'
#'
#' @param Y A m by n matrix of the Raw gene expression where m is the number of samples and n is the
#' number of features of a SummarizedExperiment variable to be normalised.
#' @param Yrep A matrix of the gene expression of the obtained PRPS
#' @param M Replicate matrix.
#' @param ctl Logical vector of length n of the negative control genes.
#' @param k A vector containing the range of k - the number of unwanted factors - to be tested.
#' @param eta A matrix with n columns, gene-wise (as opposed to sample-wise) covariates. By default is set to NULL.
#'
#' @return list List containing the corrected gene expression matrix and the M, alpha and W for each k values tested.

#' @importFrom ruv RUV1 residop
#' @export

sf.RUV_III_MK <- function(Y, Yrep, M, ctl, k = c(1:20), eta = NULL){

    ruv.adj.kmax <- ruv_III(
        Y = Y,
        Yrep = Yrep,
        M = M,
        ctl = ctl,
        k = max(k),
        eta = eta,
        return.info = TRUE
    )
    ruv.adj.kother <- lapply(
        k[k!=max(k)],
        function(x){
            ruv.adj.at.other <- ruv_III(
                Y = Y,
                Yrep = Yrep,
                M = M,
                eigVec = ruv.adj.kmax$eigenvector,
                ctl = ctl,
                k = x,
                eta = eta,
                return.info = TRUE
            )
        })
    names(ruv.adj.kother) <- paste0('RUV_K', k[k!=max(k)])
    ruv.adj.kother[[paste0('RUV_K',  max(k))]] <- ruv.adj.kmax
    return(ruv.adj.kother)
}

