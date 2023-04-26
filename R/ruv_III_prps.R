
#' is used to normalise the gene expression (assay) of a SummarizedExperiment class object
#' using RUVIII-PRPS method.
#'
#'
#' @param Y A m by n matrix of the Raw gene expression matrix where m is the number of samples and n is the
#' number of features of a SummarizedExperiment variable to be normalised.
#' @param M Replicate matrix.
#' @param ctl Logical vector of length n of the negative control genes.
#' @param k The number of unwanted factors to use.
#' @param eta A matrix with n columns, gene-wise (as opposed to sample-wise) covariates. By default is set to NULL.
#' @param include.intercept Add an intercept term to eta if it does not include one already. By default is set to TRUE.
#' @param average Average replicates after adjustment. By default is set to FALSE.
#' @param fullalpha Can be included to speed up execution. By default is set to NULL.
#' @param return.info If FALSE, only the adjusted data matrix is returned. If TRUE, additional information
#' is returned. By default is set to FALSE.
#' @param inputcheck Perform a basic sanity check on the inputs, and issue a warning if there is a problem.
#' By default is set to TRUE.
#'
#' @return list List containing the corrected gene expression matrix and the M, alpha and W.

#' @importFrom ruv RUV1 residop
#' @export


ruv_III_prps<-function(
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
  if (is.data.frame(Y) ) {
        Y <- data.matrix(Y)
    }
    m <- nrow(Y)
    n <- ncol(Y)
    M <- ruv::replicate.matrix(M)

    ctl2 <- rep(FALSE, n)
    ctl2[ctl] <- TRUE
    ctl=ctl2

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
    Y <- RUV1(Y, eta, ctl, include.intercept = include.intercept)
    mu <- colMeans(Y)
    mu_mat <- rep(1, m) %*% t(mu)
    Y_stand <- Y - mu_mat
    if (ncol(M) >= m)
        newY <- Y
    else if (is.null(k)) {
        ycyctinv <- solve(Y[, ctl] %*% t(Y[, ctl]))
        newY <- (M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*% ycyctinv)) %*% Y
        fullalpha <- NULL
    } else if (k == 0) {
        newY <- Y
        fullalpha <- NULL
    } else {
        if (is.null(fullalpha) ) {
            Y0 <- residop(Y, M)
            fullalpha <- t(svd(Y0 %*% t(Y0))$u[, 1:min(m - ncol(M), sum(ctl)), drop = FALSE]) %*% Y
        }
        alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
        ac <- alpha[, ctl, drop = FALSE]
        W <- Y_stand[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
        newY <- Y - W %*% alpha
    }
    if (average)
        newY <- ((1/apply(M, 2, sum)) * t(M)) %*% newY
    if (!return.info) {
        return(newY)
    } else {
        return(list(newY = newY, M = M, fullalpha = fullalpha,  W =  W))
    }
}

