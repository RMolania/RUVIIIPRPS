#' is used to run the RUVIII-PRPS method for a given assay and given k.
#'
#' @param se.obj A SummarizedExperiment object that will be used to computer fastRUV-III
#' @param assay.name String for the selection of the name of the assay data of the SummarizedExperiment class object
#' @param prps.data TTTTT
#' @param prps.group TTTTT
#' @param ncg.sets TTTT
#' @param ncg.approach TTTT
#' @param k TTTT
#' @param apply.log TTTT
#' @param pseudo.count TTTT
#' @param alpha TTTT
#' @param return.info TTTT
#' @param eta TTTT
#' @param include.intercept TTTT
#' @param assess.se.obj TTTT
#' @param remove.na TTTT
#' @param save.se.obj TTTT
#' @param verbose TTTT

#' @return SummarizedExperiment A SummarizedExperiment object containing the normalised gene expression
#' into a new assay (RUVIII_K).

#' @author Ramyar Molania

#' @importFrom DelayedArray t colMeans
#' @importFrom BiocSingular bsparam
#' @importFrom Matrix solve
#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @importFrom ruv replicate.matrix RUV1 residop
#' @export

RUVIII.PRPS <- function(
        se.obj,
        assay.name,
        prps.data = NULL,
        prps.group = 'supervised',
        ncg.sets = NULL,
        ncg.approach = 'supervised',
        k = 1,
        apply.log = 'both',
        pseudo.count = 1,
        alpha = NULL,
        eta = NULL,
        include.intercept = TRUE,
        return.info = FALSE,
        assess.se.obj = TRUE,
        remove.na = 'measurements',
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
    # check the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose)
    }
    # prps data ####
    ## testing all prps sets in the SummarizedExperiment object ####
    if (is.null(prps.data)) {
        if (is.null(prps.group) | length(prps.group) == 0) {
            stop('The "prps.group" should be specified ("supervised", "unsupervised" or "both").')
        } else if (!prps.group %in% c('supervised', 'unsupervised', 'both')) {
            stop('The "prps.group" should be one of the:"supervised", "unsupervised" or "both".')
        }
        if (prps.group == 'supervised') {
            prps.list <- se.obj@metadata$PRPS$supervised
        } else if (prps.group == 'unsupervised') {
            prps.list <- se.obj@metadata$PRPS$unSupervised
        } else if (prps.group == 'both') {
            all.prps <- intersect(names(se.obj@metadata$PRPS), c('supervised', 'unSupervised'))
            if (length(all.prps) != 2) {
                stop('The "Supervised" or "UnSupervised" PRPSsets are not found in the SummarizedExperiment object.')
            }
            prps.list <- lapply(
                c('Supervised', 'unsupervised'),
                function(x) se.obj@metadata$PRPS[[x]])
        }
    }
    ## using specific sets of prps in the SummarizedExperiment object ####
    if (inherits(prps.data, what = 'character')) {
        if (is.null(prps.group)) {
            stop('The "prps.group" should be specified ("supervised", "unsupervised" or "both").')
        } else if (!prps.group %in% c('supervised', 'unsupervised', 'both')) {
            stop('The "prps.group" should be one of the:"supervised", "unsupervised" or "both".')
        }
        prps.data <- lapply(
            1: length(se.obj@metadata$PRPS),
            function(x){
                if(prps.group == 'supervised'){
                    se.obj@metadata$PRPS[[prps.group]][[x]]
                } else if(prps.group == 'unsupervised'){
                    se.obj@metadata$PRPS[[prps.group]][[x]]
                } else if(prps.group == 'both'){
                    se.obj@metadata$PRPS[[prps.group]][[x]]
                    se.obj@metadata$PRPS[[prps.group]][[x]]
                }
            })
        if (se.obj@metadata$PRPS[[prps.group]][[]]) {
            stop('The prps data cannot be found in the SummarizedExperiment object.')
        }
    }
    ## using a list of PRPS data provided by users ####
    if (inherits(prps.data, what = 'list')) {
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
                    message = paste0(
                        'The ',
                        names(prps.data)[x],
                        ' contains ',
                        length(rep.samples),
                        ' PRPS sets.'),
                    color = 'blue',
                    verbose = verbose)
            })
    }

    # NCG sets ####
    ## testing all NCG sets in the SummarizedExperiment object ####
    if (is.null(ncg)) {
        if (is.null(ncg.group) | length(ncg.group) == 0) {
            stop('The "ncg.group" should be specified.')
        } else if (!prps.group %in% c('supervised', 'unsupervised', 'both')) {
            stop(
                'The "ncg.group" should be one of the:"supervised", "unsupervised" or "both".'
            )
        }
        if (ncg.group == 'supervised') {
            ncg.list <- se.obj@metadata$NCG$Supervised
        } else if (ncg.group == 'unsupervised') {
            ncg.list <- se.obj@metadata$NCG$UnSupervised
        } else if (ncg.group  == 'both') {
            all.ncg <-
                intersect(names(se.obj@metadata$NCG),
                          c('Supervised', 'UnSupervised'))
            if (length(all.ncg) != 2) {
                stop(
                    'The "Supervised" or "UnSupervised" NCGs are not found in the SummarizedExperiment object.'
                )
            }
            ncg.list <- lapply(c('Supervised', 'unsupervised'),
                               function(x)
                                   se.obj@metadata$NCG[[x]])
        }
    }
    ## using sets specific sets of NCG in the SummarizedExperiment object ####
    if (inherits(ncg, what = 'character')) {
        if (se.obj@metadata$PRPS[[ncg.approach]][[]]) {
            stop('The prps data cannot be found in the SummarizedExperiment object.')
        }
    }
    ## using a list of NCGs provided by users ####
    if (inherits(prps.data, what = 'list')) {
        ncg.list <- ncg
    }

    # check K ####
    if (min(k) <= 0)
        stop('k cannot be 0 or negative values.')

    # check max k for each normalization ####
    all.ncg.prps <- expand.grid(names(prps.list), names(ncg.list))
    colnames(all.ncg.prps) <- c('PRPS', 'NCG')
    max.k.values <- unlist(lapply(c(1:nrow(all.ncg.prps)),
                                  function(x) {
                                      M <-
                                          ruv::replicate.matrix(colnames(prps.list[[all.ncg.prps$PRPS[x]]]))
                                      min(ncol(M), sum(ncg.list[[all.ncg.prps$NCG[x]]]))
                                  }))
    all.ncg.prps$max.k = max.k.values

    # possible runs of RUV-III ####
    printColoredMessage(message = '-- Possible runs of RUV-III :',
                        color = 'magenta',
                        verbose = verbose)
    if (length(k) == 1 | is.null(k)) {
        m.out <- lapply(c(1:nrow(all.ncg.prps)),
                        function(x) {
                            if (is.null(k) | k > all.ncg.prps$max.k[x]) {
                                printColoredMessage(
                                    message = paste0(
                                        'The RUV-III-PRPS with k = ',
                                        all.ncg.prps$max.k[x],
                                        ' will be applied on the PRPS data: ',
                                        all.ncg.prps$PRPS[x],
                                        'and the NCG: ',
                                        all.ncg.prps$NCG[x],
                                        '.'
                                    ),
                                    color = 'blue',
                                    verbose = verbose
                                )
                                all.ncg.prps$max.k[x]

                            } else {
                                printColoredMessage(
                                    message = paste0(
                                        'The RUV-III-PRPS with k = ',
                                        k,
                                        ' will be applied on the PRPS data:',
                                        all.ncg.prps$PRPS[x],
                                        ' and the NCG set: ',
                                        all.ncg.prps$NCG[x],
                                        '.'
                                    ),
                                    color = 'blue',
                                    verbose = verbose
                                )
                                k
                            }
                        })
    } else {
        max.k <- max(k)
        m.out <- lapply(c(1:nrow(all.ncg.prps)),
                        function(x) {
                            if (all.ncg.prps$max.k[x] <= max.k) {
                                poss.k <- k[k %in% 1:all.ncg.prps$max.k[x]]
                                if (length(poss.k) == 1) {
                                    printColoredMessage(
                                        message = paste0(
                                            'The RUV-III-PRPS with k = ',
                                            poss.k,
                                            ' will be applied on the PRPS data: ',
                                            all.ncg.prps$PRPS[x],
                                            ' and the NCG: ',
                                            all.ncg.prps$NCG[x],
                                            '.'
                                        ),
                                        color = 'blue',
                                        verbose = verbose
                                    )
                                    poss.k

                                } else if (length(poss.k) > 1) {
                                    printColoredMessage(
                                        message = paste0(
                                            'The RUV-III-PRPS for individual k = ',
                                            paste0(poss.k, collapse = ','),
                                            ' will be applied on the PRPS data: ',
                                            all.ncg.prps$PRPS[x],
                                            ' and the NCG: ',
                                            all.ncg.prps$NCG[x],
                                            '.'
                                        ),
                                        color = 'blue',
                                        verbose = verbose
                                    )
                                    poss.k

                                }
                            } else{
                                stop('fff')
                            }
                        })
    }
    all.ncg.prps$k <- m.out
    printColoredMessage(
        message = paste0(
            'In totall ',
            sum(sapply(m.out, length)),
            ' runs of RUV-III-PRPS will appliep.'
        ),
        color = 'blue',
        verbose = verbose
    )

    # data transformation ####
    printColoredMessage(message = '-- Data transformation:',
                        color = 'magenta',
                        verbose = verbose)
    printColoredMessage(message = 'Note, make sure both the assay and PRPS data are log transformed.',
                        color = 'red',
                        verbose = verbose)
    if (apply.log == 'both') {
        printColoredMessage(
            message = paste(
                'Apply log2 + ',
                pseudo.count,
                ' (pseudo.count) on both assay and prps data.'
            ),
            color = 'blue',
            verbose = verbose
        )
        Y <- log2(assay(se.obj, assay.name) + pseudo.count)
        prps.list <-
            lapply(prps.list, function(x)
                log2(x + pseudo.count))
    } else if (apply.log == 'assay') {
        printColoredMessage(
            message = paste(
                'Apply log2 + ',
                pseudo.count,
                ' (pseudo.count) on only assay.'
            ),
            color = 'blue',
            verbose = verbose
        )
        Y <- log2(assay(se.obj, assay.name) + pseudo.count)
    } else if (apply.log == 'prps') {
        printColoredMessage(
            message = paste('Apply log2 + ', prps, ' (pseudo.count) on only assay.'),
            color = 'blue',
            verbose = verbose
        )
        prps.list <-
            lapply(prps.list, function(x)
                log2(x + pseudo.count))
    } else if (apply.log == 'none') {
        printColoredMessage(message = 'It seems both assay and PRPS are already log transformed.',
                            color = 'blue',
                            verbose = verbose)
        Y <- assay(se.obj, assay.name)
    }

    # RUVIII normalization ####
    printColoredMessage(message = '-- Apply the RUV-III-PRPS method:',
                        color = 'magenta',
                        verbose = verbose)
    all.ruv <- lapply(1:nrow(all.ncg.prps),
                      function(x) {
                          printColoredMessage(
                              message = paste0(
                                  'Run RUV-III-PRPS for individual k: ',
                                  all.ncg.prps$k[x],
                                  ' for the ',
                                  all.ncg.prps$PRPS[x],
                                  ' and the ',
                                  all.ncg.prps$NCG[x],
                                  '.'
                              ),
                              color = 'blue',
                              verbose = verbose
                          )
                          printColoredMessage(message = 'The alpha will be calculated for the maximum value of k.',
                                              color = 'blue',
                                              verbose = verbose)
                          if (length(all.ncg.prps$k[[x]]) > 1) {
                              ncg.set <- ncg.list[[all.ncg.prps$NCG[x]]]
                              ncg.logi <- rep(FALSE, nrow(se.obj))
                              ncg.logi[ncg.set] <- TRUE
                              ncg.set <- ncg.logi
                              prps.data <- prps.list[[all.ncg.prps$PRPS[x]]]
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
                              printColoredMessage(
                                  message = '-Apply RUV1 on the both assay and prps data:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              Y <- RUV1(
                                  Y = Y,
                                  eta = eta,
                                  ctl = ncg.set,
                                  include.intercept = include.intercept
                              )
                              # data standardization ####
                              printColoredMessage(
                                  message = '-Standardize the data:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              Y.stand <- scale(Y, center = TRUE, scale = FALSE)
                              printColoredMessage(
                                  message = '-Create M martix:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              M <- replicate.matrix(row.names(Y))
                              printColoredMessage(
                                  message = '-Sanity check on the M matrix:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              sum(rowSums(M) != 1)
                              sum(colSums(M) == 1)
                              printColoredMessage(
                                  message = '-Obtain residuals from the PRPS sets:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              Y0 <- fastResidop(Y, M)
                              printColoredMessage(
                                  message = '-Apply svd on the residuals to obtain alpha:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              left.sing.value <- BiocSingular::runSVD(
                                  x = Y0,
                                  k = max(k.vals),
                                  BSPARAM = bsparam(),
                                  center = FALSE,
                                  scale = FALSE
                              )$u
                              alpha <-
                                  t(left.sing.value[, 1:max(k.vals), drop = FALSE]) %*% Y
                              ac <- alpha[, ncg.set, drop = FALSE]
                              printColoredMessage(
                                  message = '-Obtain W :',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              W <-
                                  Y.stand[, ncg.set] %*% t(ac) %*% solve(ac %*% t(ac))
                              printColoredMessage(
                                  message = '-Obtain the normalized data for the maximum value of k:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              newY.max <- Y - W %*% alpha
                              newY.max <- t(newY.max[1:ncol(se.obj) ,])
                              newY.max <-
                                  list(newY = newY.max, W = W[1:ncol(se.obj) , , drop = FALSE])
                              ## other k
                              other.k <- k.vals[!k.vals %in% max(k.vals)]
                              printColoredMessage(
                                  message = '-Obtain the normalized data for the values of k:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              ruv.other.k <- lapply(other.k,
                                                    function(y) {
                                                        alpha <- t(left.sing.value[, 1:y, drop = FALSE]) %*% Y
                                                        ac <- alpha[1:y, ncg.set, drop = FALSE]
                                                        W <-
                                                            Y.stand[, ncg.set] %*% t(ac) %*% solve(ac %*% t(ac))
                                                        newY <- Y - W %*% alpha
                                                        newY <- t(newY[1:ncol(se.obj) ,])
                                                        return(list(newY = newY, W = W[1:ncol(se.obj) , , drop = FALSE]))
                                                    })
                              ruv.other.k[[length(k.vals)]] <- newY.max
                              names(ruv.other.k) <- paste0('RUV_',
                                                           all.ncg.prps$PRPS[x],
                                                           '_',
                                                           all.ncg.prps$NCG[x],
                                                           '_K:',
                                                           k.vals)
                              ruv.other.k
                          } else {
                              ncg.set <- ncg.list[[all.ncg.prps$NCG[x]]]
                              ncg.logi <- rep(FALSE, nrow(se.obj))
                              ncg.logi[ncg.set] <- TRUE
                              ncg.set <- ncg.logi
                              prps.data <- prps.list[[all.ncg.prps$PRPS[x]]]
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
                              printColoredMessage(
                                  message = 'Apply RUV1 on the both assay and prps data:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              Y <- RUV1(
                                  Y = Y,
                                  eta = eta,
                                  ctl = ncg.set,
                                  include.intercept = include.intercept
                              )
                              # data standardization ####
                              printColoredMessage(
                                  message = '-Standardize the data:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              Y.stand <- scale(Y, center = TRUE, scale = FALSE)
                              printColoredMessage(
                                  message = '-Create M martix:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              M <- replicate.matrix(row.names(Y))
                              printColoredMessage(
                                  message = '-Sanity check on the M matrix:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              sum(rowSums(M) != 1)
                              sum(colSums(M) == 1)
                              printColoredMessage(
                                  message = '-Obtain residuals from the PRPS sets:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              Y0 <- fastResidop(Y, M)
                              printColoredMessage(
                                  message = '-Apply svd on the residuals to Obtain alpha:',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              left.sing.value <- BiocSingular::runSVD(
                                  x = Y0,
                                  k = max(k.vals),
                                  BSPARAM = bsparam(),
                                  center = FALSE,
                                  scale = FALSE
                              )$u
                              alpha <-
                                  t(left.sing.value[, 1:max(k.vals), drop = FALSE]) %*% Y
                              ac <- alpha[, ncg.set, drop = FALSE]
                              printColoredMessage(
                                  message = '-Obtain W :',
                                  color = 'blue',
                                  verbose = verbose
                              )
                              W <-
                                  Y.stand[, ncg.set] %*% t(ac) %*% solve(ac %*% t(ac))
                              newY <- Y - W %*% alpha
                              newY <- t(newY[1:ncol(se.obj) ,])
                          }
                      })
    names(all.ruv) <- paste0(all.ncg.prps$PRPS, '_', all.ncg.prps$NCG)

    # save data sets ####
    if (save.se.obj) {
        ## Saving the norm data into a new assay ####
        for (x in 1:length(all.ruv)) {
            for (y in 1:length(all.ruv[[x]])) {
                new.assay.name <- names(all.ruv[[x]][y])
                se.obj@assays@data[[new.assay.name]] <-
                    all.ruv[[x]][[y]]$newY
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
                    se.obj@metadata[['RUVIII']][[new.assay.name]][['W']]  <-
                        all.ruv[[x]][[y]]$W
                }
            }
        }
        printColoredMessage(message = '------------The RUVIII.PRPS function starts:',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else if (!return.info & !save.se.obj) {
        printColoredMessage(message = '------------The RUVIII.PRPS function starts:',
                            color = 'white',
                            verbose = verbose)
        return(all.ruv)
    } else if (return.info & !save.se.obj) {
        printColoredMessage(message = '------------The RUVIII.PRPS function starts:',
                            color = 'white',
                            verbose = verbose)
        return(all.ruv)
    }
}
