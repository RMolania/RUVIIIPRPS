#' is used to run the RUVIII-PRPS method for a given assay and given k.
#'
#' @param se.obj A SummarizedExperiment object that will be used to computer fastRUV-III
#' @param assay.name String for the selection of the name of the assay data
#' of the SummarizedExperiment class object
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param replicate.data TO BE DEFINED
#' @param ctl Logical vector of length n of the negative control genes.
#' @param BSPARAM TO BE DEFINED.
#' @param k The number of unwanted factors to use.
#' @param eta A matrix with n columns, gene-wise (as opposed to sample-wise) covariates. By default is set to NULL.
#' @param include.intercept Logical. Add an intercept term to eta if it does not include one already. By default is set to TRUE.
#' @param apply.average.rep Logical. Average replicates after adjustment. By default is set to FALSE.
#' @param fullalpha Can be included to speed up execution. By default is set to NULL.
#' @param return.info Logical. If FALSE, only the adjusted data matrix is returned. If TRUE, additional information
#' is returned. By default is set to FALSE.
#' @param inputcheck Logical. Perform a basic sanity check on the inputs, and issue a warning if there is a problem.
#' By default is set to TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na TO BE DEFINED.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the normalised gene expression
#' into a new assay (RUVIII_K).

#' @importFrom DelayedArray t colMeans
#' @importFrom BiocSingular bsparam
#' @importFrom Matrix solve
#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @importFrom ruv replicate.matrix RUV1
#' @importFrom BiocSingular runSVD bsparam
#' @export


ruvIII<-function(
        se.obj,
        assay.name,
        apply.log=TRUE,
        pseudo.count = 1,
        replicate.data,
        ctl,
        BSPARAM=NULL,
        k = NULL,
        eta = NULL,
        include.intercept = TRUE,
        apply.average.rep = FALSE,
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

    if( is.null(k)){
        stop('k cannot be 0. This means no adjustment will be made.')
    } else if(min(table(rownames(replicate.data))) == 1){
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
        tBB = t(B) %*% B
        tBB_inv = solve(tBB)
        BtBB_inv = B %*% tBB_inv
        tBA = t(B) %*% A
        BtBB_inv_tBA = BtBB_inv %*% tBA
        return(A - BtBB_inv_tBA)
    }

    ## Get the expression data
    # data transformation ####
    if(apply.log){
        expr.data <- log2(assay(se.obj, assay.name) + pseudo.count)
    }else{
        expr.data <- assay(se.obj, assay.name)
    }
    Y=t(expr.data)

    # data preparation ####
    if (is.data.frame(Y) ) {
        Y <- data.matrix(Y)
    } else if (is.data.frame(replicate.data)){
        replicate.data = data.matrix(replicate.data)
    }
    m <- nrow(replicate.data)
    #m <- nrow(Y)
    m1=nrow(Y)
    #m1 = m
    n <- ncol(Y)
    M=row.names(replicate.data)
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
    mu.mat <- rep(1, m1) %*% t(mu)
    Y.stand <- Y - mu.mat
    #BSPARAM=bsparam()

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
            #k.eigVec = min(m - ncol(M), sum(ctl))
            if (is.null(BSPARAM)){
                BSPARAM=bsparam()
            }
            eigVec = runSVD(
                x = Y0,
                k = k,
                BSPARAM = BSPARAM,
                center = FALSE,
                scale = FALSE
            )$u
            # eigVec = eigen(Y0 %*% t(Y0), symmetric = TRUE)$vectors
            # fullalpha = t(eigVec[, seq_len(k.eigVec), drop = FALSE]) %*% replicate.data
            fullalpha = t(eigVec) %*% replicate.data
        }
        alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
        ac <- alpha[, ctl, drop = FALSE]
        W <- Y.stand[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
        newY <- Y - W %*% alpha
    }
    # average replicates ####
    if (apply.average.rep)
        newY <- ((1/apply(M, 2, sum)) * t(M)) %*% newY

    # Return data sets ####
    if (save.se.obj) {
        ### Saving the norm data into a new assay
        new.assay.name <- paste0('RUV_K', k, '_on_', assay.name)
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

