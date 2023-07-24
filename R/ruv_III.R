
#' is used to normalise the gene expression (assay) of a SummarizedExperiment class object
#' using RUVIII-PRPS method.
#'
#'
#' @param Y A m by n matrix of the Raw gene expression matrix where m is the number of samples and n is the
#' number of features of a SummarizedExperiment variable to be normalised.
#' @param Yrep A matrix of the gene expression of the obtained PRPS
#' @param M Replicate matrix.
#' @param ctl Logical vector of length n of the negative control genes.
#' @param k The number of unwanted factors to use.
#' @param eta A matrix with n columns, gene-wise (as opposed to sample-wise) covariates. By default is set to NULL.
#' @param Ynord A matrix with n columns, gene-wise (as opposed to sample-wise) covariates. By default is set to NULL.
#' @param eigVec A matrix with n columns, gene-wise (as opposed to sample-wise) covariates. By default is set to NULL.
#' @param include.intercept Add an intercept term to eta if it does not include one already. By default is set to TRUE.
#' @param average Average replicates after adjustment. By default is set to FALSE.
#' @param return.info If FALSE, only the adjusted data matrix is returned. If TRUE, additional information
#' is returned. By default is set to FALSE.
#' @param inputcheck Perform a basic sanity check on the inputs, and issue a warning if there is a problem.
#' By default is set to TRUE.
#'
#' @return list List containing the corrected gene expression matrix and the M, alpha and W.

#' @importFrom ruv RUV1 residop
#' @importFrom Rfast colmeans transpose
#' @importFrom BiocSingular bsparam
#' @importFrom  Rcpp sourceCpp
#' @export

ruv_III <- function (
        Y,
        Yrep,
        M,
        ctl,
        k = NULL,
        eta = NULL,
        Ynord = NULL,
        eigVec = NULL,
        include.intercept = TRUE,
        average = FALSE,
        return.info = FALSE,
        inputcheck = FALSE)
{

    ### Rcpp functions
    sourceCpp('fastCppFunctions.cpp',
              rebuild = T,
              verbose = T
    )

    m <- nrow(Y)
    n <- ncol(Y)


    tological <-function (ctl, n)
    {
        ctl2 = rep(FALSE, n)
        ctl2[ctl] = TRUE
        return(ctl2)
    }

    ctl <- tological(ctl, n)

    #message('check the inputs finished')
    if (inputcheck) {
        if (m > n)
            warning(
                "m is greater than n!  This is not a problem itself, but may
              indicate that you need to transpose your data matrix.
              Please ensure that rows correspond to observations
              (e.g. RNA-seq assay) and columns correspond to features (e.g. genes).")
        if (sum(is.na(Y)) > 0)
            warning("Y contains missing values.  This is not supported.")
        if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) >
            0)
            warning("Y contains infinities.  This is not supported.")
    }
    ############# RUV-I
    Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
    Yrep = RUV1(Yrep, eta, ctl, include.intercept = include.intercept)
    ############# Standardize
    mu <- colmeans(Y)
    mu_mat <- matrixMult(rep(1, m), matTranspose(mu))
    Y_stand <- matSubtraction(Y, mu_mat)

    if (ncol(M) >= m)
        newY = Y
    else if (is.null(k)) {
        ycyctinv = solve(matrixMult(Y[, ctl], matTranspose(Y[, ctl])))
        newY = matrixMult(
            (M %*% solve(matTranspose(M) %*% ycyctinv %*% M) %*% (matrixMult(matTranspose(M) ,  ycyctinv))),
            Y)
        fullalpha = NULL
    }
    else if (k == 0) {
        newY = Y
        fullalpha = NULL
    }
    else {
        if (is.null(Ynord) & is.null(eigVec)) {
            ## residual operators
            Y0 = fastResidop(Yrep, M)
            ## eigen vectors
            eigenvector <- BiocSingular::runSVD(
                x = matrixMult(Y0, matTranspose(Y0)),
                k = k,
                BSPARAM = bsparam(),
                center = TRUE,
                scale = FALSE
            )$u
            ## fullalpha
            fullalpha = matrixMult(
                matTranspose(eigenvector[, 1:k, drop = FALSE]),
                Yrep
            )
        }
        if (!is.null(Ynord) & is.null(eigVec)) {
            ## eigen vectors
            eigenvector <- BiocSingular::runSVD(
                x = matrixMult(Ynord, Rfast::transpose(Ynord)),
                k = k,
                BSPARAM = bsparam(),
                center = TRUE,
                scale = FALSE
            )$u
            ## fullalpha
            fullalpha = matrixMult(transpose(eigenvector[, 1:k, drop = FALSE]),  Yrep)
        }
        if (is.null(Ynord) & !is.null(eigVec)) {
            fullalpha = matrixMult(transpose(eigVec[, 1:k, drop = FALSE]),  Yrep)
        }
        ## alpha
        alpha = fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
        ac = alpha[, ctl, drop = FALSE]
        ## W
        W <- matrixMult(matrixMult(Y_stand[ , ctl], transpose(ac)),
                        Matrix::solve(matrixMult(ac, transpose(ac))))
        ## data adjustment
        newY = matSubtraction(Y, matrixMult(W, alpha))
    }
    if (average)
        newY = ((1 / apply(M, 2, sum)) * t(M)) %*% newY
    if (!return.info)
        return(newY)
    if (is.null(Ynord) & is.null(eigVec))
        return(list(
            newY = newY,
            eigenvector = eigenvector,
            W = W,
            Ynord = Y0
        ))
    if (is.null(eigVec))
        return(list(
            newY = newY,
            eigenvector = eigenvector,
            W = W
        ))
    else
        return(list(
            newY = newY,
            W = W))
}
