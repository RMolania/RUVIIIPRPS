#' RUVIII-PRPS normalization.

#' @author Ramyar Molania

#' @description
#' This function applies the RUV-III-PRPS on RNA-seq data.

#' @details
#' Additional details...
#'

#' @references
#' Molania R., ..., Speed, T. P., A new normalization for Nanostring nCounter gene expression data, Nucleic Acids Research,
#' 2019.
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023


#' @param se.obj A SummarizedExperiment object that will be used to computer fastRUV-III
#' @param assay.name String for the selection of the name of the assay data of the SummarizedExperiment class object
#' @param prps.set.names TTTTT
#' @param prps.group TTTTT
#' @param ncg.set.names TTTT
#' @param ncg.group TTTT
#' @param k The number of unwanted factors to use. Can be 0, in which case no adjustment is made. Can also be NULL (the
#' default value), in which case the maximum possible value of k is used; note that in this case no singular value
#' decomposition is necessary and execution is faster.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param alpha TTTT
#' @param return.info If FALSE, only the adjusted data matrix is returned. If TRUE, additional information is returned
#' (see below).
#' @param eta Gene-wise (as opposed to sample-wise) covariates. These covariates are adjusted for by RUV-1 before any
#' further analysis proceeds. Can be either (1) a matrix with n columns, (2) a matrix with n rows, (3) a dataframe with
#' n rows, (4) a vector or factor of length n, or (5) simply 1, for an intercept term.
#' @param include.intercept When eta is specified (not NULL) but does not already include an intercept term, this will
#' automatically include one.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the RUV-III-PRPS normalized data as assay(s) in the
#' SummarizedExperiment object or to output the result as list. By default it is set to 'TRUE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return SummarizedExperiment A SummarizedExperiment object containing the normalised gene expression
#' into a new assay (RUVIII_K).

#' @importFrom DelayedArray t colMeans
#' @importFrom BiocSingular bsparam
#' @importFrom Matrix solve
#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @importFrom ruv replicate.matrix RUV1 residop
#' @export

RUVIII.PRPS <- function(
        se.obj,
        assay.name,
        prps.group = 'supervised',
        prps.set.names = NULL,
        ncg.group = 'supervised',
        ncg.set.names = NULL,
        k = 1,
        apply.log = 'assay',
        pseudo.count = 1,
        alpha = NULL,
        eta = NULL,
        include.intercept = TRUE,
        return.info = FALSE,
        assess.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The RUVIII.PRPS function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if (length(assay.name) > 1) {
        stop('The "assay.name" should contain only one assay name.')
    }
    if(is.logical(apply.log)){
        stop('The "apply.log" must be one of the "assay", "prps", "both", "none" ')
    }

    # check the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose)
    }

    # prps set data ####
    ## use all the prps sets in the SummarizedExperiment object ####
    if (is.null(prps.set.names)) {
        if (is.null(prps.group) | length(prps.group) == 0) {
            stop('The "prps.group" should be specified ("supervised", "unSupervised" or "both").')
        } else if (!prps.group %in% c('supervised', 'unSupervised', 'both')) {
            stop('The "prps.group" should be one of the "supervised", "unSupervised" or "both".')
        }
        if (prps.group == 'supervised') {
            prps.data <- se.obj@metadata$PRPS$supervised
            printColoredMessage(
                message = paste0(length(prps.data), ' supervised PRPS set(s) are found in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
        } else if (prps.group == 'unSupervised') {
            prps.data <- se.obj@metadata$PRPS$unSupervised
            printColoredMessage(
                message = paste0(length(prps.data), ' unSupervised PRPS set(s) are found in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
        } else if (prps.group == 'both') {
            all.prps <- intersect(names(se.obj@metadata$PRPS), c('supervised', 'unSupervised'))
            if (length(all.prps) != 2) {
                stop('The "supervised" or "unSupervised" PRPSsets are not found in the SummarizedExperiment object.')
            }
            prps.data <- c(se.obj@metadata$PRPS$supervised, se.obj@metadata$PRPS$unSupervised)
        }
    }
    ## use a specific set of prps in the SummarizedExperiment object ####
    if (inherits(prps.set.names, what = 'character')) {
        if (is.null(prps.group)) {
            stop('The "prps.group" should be specified ("supervised" or "unSupervised").')
        } else if (!prps.group %in% c('supervised', 'unsupervised')) {
            stop('The "prps.group" should be one of the "supervised" or "unSupervised".')
        }
        if (!intersect(names(se.obj@metadata$PRPS), prps.group) == prps.group) {
            stop('The "prps.group" is not found in the SummarizedExperiment object.')
        }
        if (!prps.set.names %in% names(se.obj@metadata$PRPS[[prps.group]])) {
            stop('The "prps.set.names" is not found in the SummarizedExperiment object.')
        }
        prps.data <- list(prps = se.obj@metadata$PRPS[[prps.group]][[prps.set.names]])
        names(prps.data) <- prps.set.names
    }
    ## use a list of PRPS sets in the SummarizedExperiment object ####
    if (inherits(prps.set.names, what = 'list')) {
        if (is.null(prps.group)) {
            stop('The "prps.group" should be specified ("supervised" or "unSupervised").')
        } else if (!prps.group %in% c('supervised', 'unsupervised')) {
            stop('The "prps.group" should be one of the "supervised" or "unSupervised".')
        }
        if (!intersect(names(se.obj@metadata$PRPS), prps.group) == prps.group) {
            stop('The "prps.group" is not found in the SummarizedExperiment object.')
        }
        m.out <- lapply(
            prps.set.names,
            function(x) {
                if (!x %in% names(se.obj@metadata$PRPS[[prps.group]]))
                    stop(paste0('The prps ste', x, 'is not found in the SummarizedExperiment object.'))
                })
        prps.data <- lapply(
            prps.set.names,
            function(x) se.obj@metadata$PRPS[[prps.group]][[x]])
        names(prps.data) <- prps.set.names
    }

    ## check prps data ####
    m.out <- lapply(
        1:length(prps.data),
        function(x) {
            if (nrow(se.obj) != nrow(prps.data[[x]])) {
                stop('The number of genes between the SummarizedExperiment object and the prps data is not the same.')
            }
            if (!all.equal(row.names(se.obj), row.names(prps.data[[x]]))) {
                stop('The order of genes are not the same between the SummarizedExperiment object and the PRPS data.')
            }
            rep.samples <- findRepeatingPatterns(
                vec = colnames(prps.data[[x]]),
                n.repeat = 2)
            if (length(rep.samples) == 0) {
                stop('The names of the columns in the prps.data are all unique.')
            } else if (length(rep.samples) != length(unique(colnames(prps.data[[x]])))) {
                stop('Some names of the columns in the prps.data are all unique.')
            }
            printColoredMessage(
                message = paste0('The "', names(prps.data)[x], '" contains ', length(rep.samples),' PRPS sets.'),
                color = 'blue',
                verbose = verbose)
        })
    prps.data <- do.call(cbind, prps.data)

    # ncg sets ####
    ## use all NCG sets in the SummarizedExperiment object ####
    if (is.null(ncg.set.names)) {
        if (is.null(ncg.group) | length(ncg.group) == 0) {
            stop('The "ncg.group" must be specified.')
        } else if (!ncg.group %in% c('supervised', 'unSupervised', 'both')) {
            stop('The "ncg.group" must be one of the:"supervised", "unSupervised" or "both".')
        }
        if (ncg.group == 'supervised') {
            ncg.list <- se.obj@metadata$NCG$supervised
        } else if (ncg.group == 'unsupervised') {
            ncg.list <- se.obj@metadata$NCG$unSupervised
        } else if (ncg.group  == 'both') {
            all.ncg <- intersect(names(se.obj@metadata$NCG), c('supervised', 'unSupervised'))
            if (length(all.ncg) != 2) {
                stop('The "Supervised" or "UnSupervised" NCGs are not found in the SummarizedExperiment object.')
            }
            ncg.list <- c(se.obj@metadata$NCG$supervised, se.obj@metadata$NCG$unSupervised)
        }
    }
    ## use a specific set of NCG in the SummarizedExperiment object ####
    if (inherits(ncg.set.names, what = 'character')) {
        if (is.null(ncg.group)) {
            stop('The "ncg.group" should be specified ("supervised" or "unSupervised").')
        } else if (!ncg.group %in% c('supervised', 'unsupervised')) {
            stop('The "prps.group" should be one of the "supervised" or "unSupervised".')
        }
        if (!intersect(names(se.obj@metadata$NCG), ncg.group) == ncg.group) {
            stop('The "ncg.group" is not found in the SummarizedExperiment object.')
        }
        if (!ncg.set.names %in% names(se.obj@metadata$NCG[[ncg.group]])) {
            stop('The "ncg.set.names" is not found in the SummarizedExperiment object.')
        }
        ncg.list <- list(ncg = se.obj@metadata$NCG[[ncg.group]][[ncg.set.names]])
        names(ncg.list) <- ncg.set.names
    }

    ## use a list of PRPS sets in the SummarizedExperiment object ####
    if (inherits(ncg.set.names, what = 'list')) {
        if (is.null(ncg.group)) {
            stop('The "ncg.group" should be specified ("supervised" or "unSupervised").')
        } else if (!ncg.group %in% c('supervised', 'unsupervised')) {
            stop('The "ncg.group" should be one of the "supervised" or "unSupervised".')
        }
        if (!intersect(names(se.obj@metadata$NCG), ncg.group) == ncg.group) {
            stop('The "ncg.group" is not found in the SummarizedExperiment object.')
        }
        m.out <- lapply(
            ncg.set.names,
            function(x) {
                if (!x %in% names(se.obj@metadata$NCG[[ncg.group]]))
                    stop(paste0('The ncg ste', x, 'is not found in the SummarizedExperiment object.'))
            })
        ncg.list <- lapply(
            ncg.set.names,
            function(x) se.obj@metadata$NCG[[ncg.group]][[x]])
        names(ncg.list) <- ncg.set.names
    }

    # check K ####
    if (min(k) <= 0) stop('k cannot be 0 or negative values.')

    # check max k for each normalization ####
    all.ncg.prps <- expand.grid('prps.data', names(ncg.list))
    colnames(all.ncg.prps) <- c('PRPS', 'NCG')
    max.k.values <- unlist(lapply(
        c(1:nrow(all.ncg.prps)),
        function(x) {
            M <- ruv::replicate.matrix(colnames(prps.data))
            min(ncol(M), sum(ncg.list[[all.ncg.prps$NCG[x]]]))
        }))
    all.ncg.prps$max.k <- max.k.values
    m.out <- lapply(
        1:nrow(all.ncg.prps),
        function(x)
            print( paste0(
                    'The maximum k for the curret PRPS data',
                    ' and the NCG set:',
                    all.ncg.prps$NCG[x],
                    ' is ',
                    all.ncg.prps$max.k[x], '.')
            ))

    # possible runs of RUV-III ####
    printColoredMessage(message = '-- Possible runs of the RUV-III-PRPS method:',
                        color = 'magenta',
                        verbose = verbose)
    if(is.null(k)){
        final.k <- lapply(
            c(1:nrow(all.ncg.prps)),
            function(x){
                printColoredMessage(
                    message = paste0(
                        'The RUV-III-PRPS with k = ', all.ncg.prps$max.k[x],
                        ' will be applied on the current PRPS data ', 'and the NCG: ', all.ncg.prps$NCG[x], '.'),
                    color = 'blue',
                    verbose = verbose)
                all.ncg.prps$max.k[x]
            })
    } else if (length(k) == 1 ){
        final.k <- lapply(
            c(1:nrow(all.ncg.prps)),
            function(x) {
                if ( k > all.ncg.prps$max.k[x]) {
                    printColoredMessage(
                        message = paste0('The RUV-III-PRPS with k = ',  all.ncg.prps$max.k[x],
                            ' will be applied on the current PRPS data ',
                            'and the NCG: ', all.ncg.prps$NCG[x], '.'),
                        color = 'blue',
                        verbose = verbose)
                    all.ncg.prps$max.k[x]
                } else {
                    printColoredMessage(
                        message = paste0('The RUV-III-PRPS with k = ', k,
                            ' will be applied on the PRPS data:', all.ncg.prps$PRPS[x], ' and the NCG set: ',
                            all.ncg.prps$NCG[x],'.'),
                        color = 'blue',
                        verbose = verbose)
                    k
                }
            })
    } else if (length(k) > 1){
        max.k <- max(k)
        final.k <- lapply(
            c(1:nrow(all.ncg.prps)),
            function(x) {
                if (all.ncg.prps$max.k[x] <= max.k) {
                    poss.k <- k[k %in% 1:all.ncg.prps$max.k[x]]
                    if (length(poss.k) == 1) {
                        printColoredMessage(
                            message = paste0( 'The RUV-III-PRPS with k = ', poss.k,
                                ' will be applied on the current PRPS data ', ' and the NCG: ', all.ncg.prps$NCG[x], '.'),
                            color = 'blue',
                            verbose = verbose)
                        poss.k
                    } else if (length(poss.k) > 1) {
                        printColoredMessage(
                            message = paste0('The RUV-III-PRPS for individual k = ', paste0(poss.k, collapse = ','),
                                ' will be applied on the current PRPS data: ', ' and the NCG: ', all.ncg.prps$NCG[x], '.'),
                            color = 'blue',
                            verbose = verbose)
                        poss.k
                    }
                } else{
                    printColoredMessage(
                        message = paste0('The RUV-III-PRPS for individual k = ', paste0(k, collapse = ','),
                            ' will be applied on the current PRPS data: ', ' and the NCG set: ', all.ncg.prps$NCG[x], '.'),
                        color = 'blue',
                        verbose = verbose)
                    k
                }
            })

    }

    all.ncg.prps$k <- final.k
    printColoredMessage(
        message = paste0('In totall ', sum(sapply(final.k, length)), ' runs of RUV-III-PRPS will applied.'),
        color = 'blue',
        verbose = verbose)

    # data transformation ####
    if (is.null(pseudo.count)) pseudo.count = 0
    printColoredMessage(message = '-- Data transformation:',
                        color = 'magenta',
                        verbose = verbose)
    printColoredMessage(message = 'Note, make sure both the assay and PRPS data are log transformed.',
                        color = 'red',
                        verbose = verbose)
    if (apply.log == 'both') {
        printColoredMessage(
            message = paste('Apply log2 + ', pseudo.count, ' (pseudo.count) on both assay and prps data.'),
            color = 'blue',
            verbose = verbose
        )
        Y <- log2(assay(se.obj, assay.name) + pseudo.count)
        prps.data <- log2(prps.data + pseudo.count)
    } else if (apply.log == 'assay') {
        printColoredMessage(
            message = paste0('Apply log2 + ', pseudo.count, ' (pseudo.count) on only assay.'),
            color = 'blue',
            verbose = verbose
        )
        Y <- log2(assay(se.obj, assay.name) + pseudo.count)
    } else if (apply.log == 'prps') {
        printColoredMessage(
            message = paste('Apply log2 + ', pseudo.count, ' (pseudo.count) on only assay.'),
            color = 'blue',
            verbose = verbose
        )
        prps.data <- log2(prps.data + pseudo.count)
    } else if (apply.log == 'none') {
        printColoredMessage(message = 'It seems both assay and PRPS are already log transformed.',
                            color = 'blue',
                            verbose = verbose)
        Y <- assay(x = se.obj, i = assay.name)
    }

    # RUVIII normalization ####
    printColoredMessage(message = '-- Apply the RUV-III-PRPS method:',
                        color = 'magenta',
                        verbose = verbose)
    all.ruv <- lapply(
        1:nrow(all.ncg.prps),
        function(x) {
            printColoredMessage(
                message = paste0(
                    'Run RUV-III-PRPS for individual k: ',
                    all.ncg.prps$k[x],
                    ' for the current PRPS data',
                    ' and the NCG set: ',
                    all.ncg.prps$NCG[x],
                    '.'),
                color = 'blue',
                verbose = verbose
            )
            if (length(all.ncg.prps$k[[x]]) > 1) {
                printColoredMessage(message = 'The alpha will be calculated for the maximum value of k.',
                                    color = 'blue',
                                    verbose = verbose)
                # ncg ####
                ncg.set <- ncg.list[[all.ncg.prps$NCG[x]]]
                ncg.logi <- rep(FALSE, nrow(se.obj))
                ncg.logi[ncg.set] <- TRUE
                ncg.set <- ncg.logi
                # prps dara ####
                prps.data <- prps.data
                k.vals <- all.ncg.prps$k[[x]]
                # apply RUV1 on Y ####
                Y <- t(cbind(Y, prps.data))
                if (is.data.frame(Y))
                    Y <- data.matrix(Y)
                if (sum(is.na(Y)) > 0) {
                    stop("The assay or PRPS data contains missing values. This is not supported.")
                } else if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 0) {
                    stop("The assay or PRPS data contains infinity values. This is not supported.")
                }
                printColoredMessage(message = '-Apply RUV1 on the both assay and prps data:',
                                    color = 'blue',
                                    verbose = verbose)
                Y <- ruv::RUV1(
                    Y = Y,
                    eta = eta,
                    ctl = ncg.set,
                    include.intercept = include.intercept
                )
                # data standardization ####
                printColoredMessage(message = '-Standardize the data:',
                                    color = 'blue',
                                    verbose = verbose)
                Y.stand <- scale(Y, center = TRUE, scale = FALSE)
                printColoredMessage(message = '-Create M martix:',
                                    color = 'blue',
                                    verbose = verbose)
                M <- replicate.matrix(row.names(Y))
                printColoredMessage(message = '-Sanity check on the M matrix:',
                                    color = 'blue',
                                    verbose = verbose)
                sum(rowSums(M) != 1)
                sum(colSums(M) == 1)
                printColoredMessage(message = '-Obtain residuals from the PRPS sets:',
                                    color = 'blue',
                                    verbose = verbose)
                Y0 <- fastResidop2(Y, M)
                printColoredMessage(message = '-Apply svd on the residuals to obtain alpha:',
                                    color = 'blue',
                                    verbose = verbose)
                left.sing.value <- BiocSingular::runSVD(
                        x = Y0,
                        k = max(k.vals),
                        BSPARAM = bsparam(),
                        center = FALSE,
                        scale = FALSE
                    )$u
                alpha <- t(left.sing.value[, 1:max(k.vals), drop = FALSE]) %*% Y
                ac <- alpha[, ncg.set, drop = FALSE]
                printColoredMessage(message = '-Obtain W :',
                                    color = 'blue',
                                    verbose = verbose)
                W <- Y.stand[, ncg.set] %*% t(ac) %*% solve(ac %*% t(ac))
                printColoredMessage(message = '-Obtain the normalized data for the maximum value of k:',
                                    color = 'blue',
                                    verbose = verbose)
                newY.max <- Y - W %*% alpha
                newY.max <- t(newY.max[1:ncol(se.obj) ,])
                newY.max <- list(newY = newY.max, W = W[1:ncol(se.obj) , , drop = FALSE])
                ## other k
                other.k <- k.vals[!k.vals %in% max(k.vals)]
                printColoredMessage(message = '-Obtain the normalized data for the values of k:',
                                    color = 'blue',
                                    verbose = verbose)
                ruv.other.k <- lapply(
                    other.k,
                    function(y) {
                        alpha <- t(left.sing.value[, 1:y, drop = FALSE]) %*% Y
                        ac <- alpha[1:y, ncg.set, drop = FALSE]
                        W <- Y.stand[, ncg.set] %*% t(ac) %*% solve(ac %*% t(ac))
                        newY <- Y - W %*% alpha
                        newY <- t(newY[1:ncol(se.obj) , ])
                        return(list(newY = newY, W = W[1:ncol(se.obj) , , drop = FALSE]))
                    })
                ruv.other.k[[length(k.vals)]] <- newY.max
                names(ruv.other.k) <- paste0('RUVIIIRPPS_K.', k.vals)
                ruv.other.k
            } else {
                ncg.set <- ncg.list[[all.ncg.prps$NCG[x]]]
                ncg.logi <- rep(FALSE, nrow(se.obj))
                ncg.logi[ncg.set] <- TRUE
                ncg.set <- ncg.logi
                prps.data <- prps.data
                k.vals <- all.ncg.prps$k[[x]]
                # apply RUV1 on Y ####
                Y <- t(cbind(Y, prps.data))
                if (is.data.frame(Y))
                    Y <- data.matrix(Y)
                if (sum(is.na(Y)) > 0) {
                    stop("The assay or PRPS data contains missing values. This is not supported.")
                } else if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 0) {
                    stop("The assay or PRPS data contains infinity values. This is not supported.")
                }
                printColoredMessage(message = 'Apply RUV1 on the both assay and prps data:',
                                    color = 'blue',
                                    verbose = verbose)
                Y <- RUV1(
                    Y = Y,
                    eta = eta,
                    ctl = ncg.set,
                    include.intercept = include.intercept
                )
                # data standardization ####
                printColoredMessage(message = '-Standardize the data:',
                                    color = 'blue',
                                    verbose = verbose)
                Y.stand <- scale(Y, center = TRUE, scale = FALSE)
                printColoredMessage(message = '-Create M martix:',
                                    color = 'blue',
                                    verbose = verbose)
                M <- replicate.matrix(row.names(Y))
                printColoredMessage(message = '-Sanity check on the M matrix:',
                                    color = 'blue',
                                    verbose = verbose)
                sum(rowSums(M) != 1)
                sum(colSums(M) == 1)
                printColoredMessage(message = '-Obtain residuals from the PRPS sets:',
                                    color = 'blue',
                                    verbose = verbose)
                Y0 <- fastResidop2(Y, M)
                printColoredMessage(message = '-Apply svd on the residuals to Obtain alpha:',
                                    color = 'blue',
                                    verbose = verbose)
                left.sing.value <- BiocSingular::runSVD(
                        x = Y0,
                        k = max(k.vals),
                        BSPARAM = bsparam(),
                        center = FALSE,
                        scale = FALSE
                    )$u
                alpha <- t(left.sing.value[, 1:k.vals, drop = FALSE]) %*% Y
                ac <- alpha[, ncg.set, drop = FALSE]
                printColoredMessage(message = '-Obtain W :',
                                    color = 'blue',
                                    verbose = verbose)
                W <- Y.stand[, ncg.set] %*% t(ac) %*% solve(ac %*% t(ac))
                newY <- Y - W %*% alpha
                newY <- t(newY[1:ncol(se.obj) ,])
                newY <- list(newY = newY, W = W[1:ncol(se.obj) , , drop = FALSE])
            }
        })
    names(all.ruv) <- paste0(all.ncg.prps$PRPS, '_', all.ncg.prps$NCG)

    # save data sets ####
    if (save.se.obj) {
        ## Saving the norm data into a new assay ####
        for (x in 1:length(all.ruv)) {
            for (y in 1:length(all.ruv[[x]])) {
                new.assay.name <- names(all.ruv[[x]][y])
                se.obj@assays@data[[new.assay.name]] <- all.ruv[[x]][[y]]$newY
            }
        }
        printColoredMessage(
            message = paste0(
                'The normalized data ',
                new.assay.name,
                ' is saved into a new assay of the SummarizedExperiment object.'
            ),
            color = 'blue',
            verbose = verbose
        )

        # saving the W and alpha as well ####
        if (return.info) {
            ## Check if metadata RUVIII already exist
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['RUVIII']] <- list()
            }
            ## Check if metadata RUVIII already exist
            if (!'RUVIII' %in% names(se.obj@metadata)) {
                se.obj@metadata[['RUVIII']] <- list()
            }
            for (x in 1:length(all.ruv)) {
                for (y in 1:length(all.ruv[[x]])) {
                    new.assay.name <- names(all.ruv[[x]][y])
                    ## Check if metadata RUVIII already exist for this assay
                    if (!new.assay.name %in% names(se.obj@metadata[['RUVIII']])) {
                        se.obj@metadata[['RUVIII']][[new.assay.name]] <- list()
                    }
                    ## Check if metadata RUVIII already exist for this assay for W
                    if (!'W' %in% names(se.obj@metadata[['RUVIII']][[new.assay.name]])) {
                        se.obj@metadata[['RUVIII']][[new.assay.name]][['W']]  <- list()
                    }
                    se.obj@metadata[['RUVIII']][[new.assay.name]][['W']]  <- all.ruv[[x]][[y]]$W
                }
            }
        }
        printColoredMessage(message = '------------The RUVIII.PRPS function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else if (!return.info & !save.se.obj) {
        printColoredMessage(message = '------------The RUVIII.PRPS function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.ruv)
    } else if (return.info & !save.se.obj) {
        printColoredMessage(message = '------------The RUVIII.PRPS function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.ruv)
    }
}
