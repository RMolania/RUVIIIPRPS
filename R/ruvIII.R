#' is used to normalise the gene expression (assay) of a SummarizedExperiment class object
#' using the RUVIII-PRPS method.
#'
#' @param se.obj A SummarizedExperiment object that will be used to computer fastRUV-III
#' @param assay.name String for the selection of the name of the assay data
#' of the SummarizedExperiment class object
#' @param replicate.data TO BE DEFINED
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
#' @param assess.se.obj Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na TO BE DEFINED.
#' @param save.se.obj Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the normalised gene expression
#' into a new assay (RUVIII_K).

#' @importFrom DelayedArray t colMeans
#' @importFrom BiocSingular bsparam
#' @importFrom Matrix solve
#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @importFrom ruv replicate.matrix RUV1
#' @importFrom BiocSingular runSVD
#' @export


ruvIII<-function(
        se.obj,
        assay.name,
        replicate.data,
        Y,
        M,
        ctl,
        k = NULL,
        eta = NULL,
        include.intercept = TRUE,
        average = FALSE,
        fullalpha = NULL,
        return.info = FALSE,
        inputcheck = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'measurements',
        save.se.obj = TRUE,
        verbose = TRUE
){

    printColoredMessage(message = '------------The ruvIII function starts:',
                        color = 'white',
                        verbose = verbose)

    if(k == 0 || is.null(k)){
        stop('k cannot be 0. This means no adjustment will be made.')
    } else if(min(table(colnames(replicate.data))) == 1){
        stop('There are only replicated samples of a single sample in the replicate.data')
    }

    # check the SummarizedExperiment ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose
        )
    }

    ## Functions
    tological <- function(ctl, n) {
        ctl2 <- rep(FALSE, n)
        ctl2[ctl] <- TRUE
        return(ctl2)
    }
    fast.residop <- function(A, B){
        tBB = DelayedArray::t(B) %*% B
        tBB_inv = Matrix::solve(tBB)
        BtBB_inv = B %*% tBB_inv
        tBA = DelayedArray::t(B) %*% A
        BtBB_inv_tBA = BtBB_inv %*% tBA
        return(A - BtBB_inv_tBA)
    }

    # data preparation ####
    if (is.data.frame(Y) ) {
        Y <- data.matrix(Y)
    } else if (is.data.frame(replicate.data)){
        replicate.data = data.matrix(replicate.data)
    }
    m <- nrow(Y)
    m1 = m
    n <- ncol(Y)
    M <- replicate.matrix(M)
    ctl <- tological(ctl, n)

    if (inputcheck) {
        if (m > n)
            warning(
                "Please ensure that rows correspond to observations
              (e.g. RNA-seq assay) and columns correspond to features (e.g. genes).")
        if (sum(is.na(Y)) > 0)
            warning("Y contains missing values. This is not supported.")
        if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) >
            0)
            warning("Y contains infinity values. This is not supported.")
    }
    # run RUV1 ####
    Y <- RUV1(Y, eta, ctl, include.intercept = include.intercept)
    replicate.data = RUV1(replicate.data, eta, ctl, include.intercept = include.intercept)
    # data standardization ####
    mu <- colMeans(Y)
    mu.mat <- rep(1, m) %*% t(mu)
    Y.stand <- Y - mu.mat
    BSPARAM=bsparam()

    # RUVIII normalization ####
    printColoredMessage(
        message = '### Performing the RUVIII normalization',
        color = 'magenta',
        verbose = verbose)
    if (ncol(M) >= m){
        newY <- Y
    } else if (is.null(k)) {
        ycyctinv <- solve(Y[, ctl] %*% t(Y[, ctl]))
        newY <- (M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*% ycyctinv)) %*% Y
        fullalpha <- NULL
    } else if (k == 0) {
        newY <- Y
        fullalpha <- NULL
    } else {
        if (is.null(fullalpha) ) {
            Y0 = fast.residop(replicate.data, M)
            k.eigVec = min(m - ncol(M), sum(ctl))
            eigVec = runSVD(
                x = Y0,
                k = k,
                BSPARAM = BSPARAM,
                center = FALSE,
                scale = FALSE
            )$u
            fullalpha = t(eigVec[, seq_len(k.eigVec), drop = FALSE]) %*% replicate.data
        }
        alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
        ac <- alpha[, ctl, drop = FALSE]
        W <- Y.stand[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
        newY <- Y - W %*% alpha
    }
    # average replicates ####
    if (average)
        newY <- ((1/apply(M, 2, sum)) * t(M)) %*% newY

    # Return data sets ####
    if (save.se.obj) {
        ### Saving the norm data into a new assay
        new.assay.name <- paste0('RUVIII_K:', k, '_Data:', assay.name)
        if(!new.assay.name %in% (names(se.obj@assays@data)) ){
            se.obj@assays@data[[new.assay.name]] <- t(newY)
        }
        printColoredMessage(message= paste0(
            'The normalized data ', new.assay.name, ' is saved into a new assay of the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose)

        ### Saving the W and fullalpha as well
        if(return.info) {

            ## Check if metadata RUVIII already exist
            if(length(se.obj@metadata)==0 ) {
                se.obj@metadata[['RUVIII']] <- list()
            }
            ## Check if metadata RUVIII already exist
            if(!'RUVIII' %in% names(se.obj@metadata) ) {
                se.obj@metadata[['RUVIII']] <- list()
            }
            ## Check if metadata RUVIII already exist for this assay
            if(!new.assay.name %in% names(se.obj@metadata[['RUVIII']]) ) {
                se.obj@metadata[['RUVIII']][[new.assay.name]] <- list()
            }
            ## Check if metadata RUVIII already exist for this assay for W
            if(!'W'%in% names(se.obj@metadata[['RUVIII']][[new.assay.name]]) ) {
                se.obj@metadata[['RUVIII']][[new.assay.name]][['W']]  <- W
            }
            ## Check if metadata RUVIII already exist for this assay for fullalpha
            if(!'fullalpha'%in% names(se.obj@metadata[['RUVIII']][[new.assay.name]]) ) {
                se.obj@metadata[['RUVIII']][[new.assay.name]] [['fullalpha']] <- fullalpha
            }
        }
        return(se.obj)
    } else if(!return.info & !save.se.obj){
        return(Y=t(newY))
    } else if(return.info & !save.se.obj){
        return(list(
            newY = newY,
            fullalpha = fullalpha,
            W = W))
    }
    printColoredMessage(message = '------------The ruvIII function finished:',
                        color = 'white',
                        verbose = verbose)

}

