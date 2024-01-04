#' is used to assess the association between variables in a SummarizedExperiment object.
#'
#'
#' @description
#' The function assesses the association between all selected biological and unwanted variation variables separately. If two categorical
#' variables are highly correlated, the function keeps one that has the highest number of factors. For two continuous variables,
#' the one with higher variance will be kept.
#'
#'
#' @param se.obj A SummarizedExperiment object.
#' @param bio.variables String of the label of (a) categorical variable(s) that specifies major biological groups
#' such as samples types from colData(se).
#' @param uv.variables String or vector of strings of the label of continuous or categorical variable(s)
#' such as samples types, batch or library size from colData(se) that will be used to define PRPS.
#' @param cat.cor.coef Vector of two numerical values. Indicates the cut-off of the correlation coefficient between each
#' pair of categorical variables. The first one is between each pair of 'uv.variables' and the second one is between each
#' pair of 'bio.variables'. The correlation is computed by the function ContCoef from the DescTools package. If the
#' correlation of a pair of variable is higher than the cut-off, then only the variable that has the highest number of
#' factor will be kept and the other one will be excluded from the remaining analysis. By default they are both set to 0.7.
#' @param cont.cor.coef Vector of two numerical values. Indicates the cut-off of the Spearman correlation coefficient
#' between each pair of continuous variables. The first one is between each pair of 'uv.variables' and the second one is
#' between each pair of 'bio.variables'. If the correlation of a pair of variable is higher than the cut-off, then only
#' the variable that has the highest variance will be kept and the other one will be excluded from the remaining analysis.
#' By default they are both set to 0.7.
#' @param assess.se.obj Logical. Whether to assess the SummarizedExperiment object or not.
#' @param remove.na String. Indicates whether to remove NA or missing values from either the 'sample.annotation', 'both'
#' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be excluded. If 'sample.annotation'
#' is selected, the samples that contains NA or missing values for any 'bio.variables' and 'uv.variables' will be excluded.
#' By default, it is set to
#' 'sample.annotation'.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed
#' during the execution of the functions, by default it is set to TRUE.
#'
#'
#'
#' @details
#' For each pair of categorical variables from 'uv.variables' and each pair of categorical variables from 'bio.variables',
#' the correlation is computed using the function ContCoef from the DescTools R package. For each pair of continuous variables
#' from 'uv.variables' and each pair of continuous variables from 'bio.variables', the correlation is computed using the
#' Spearman correlation. The user defines a minimum cut-off of the correlation coefficient between each pair of categorical
#' variables in the 'cat.cor.coef' for the unwanted variables, followed by the minimum cut-off of the correlation coefficient
#' for the biological variables. If the correlation of a pair of those variables is higher than the minimum cut-off, only
#' the variable that has the highest number of factor will be kept and the other one will be excluded from the remaining analysis.
#' The user defines a minimum cut-off of the correlation coefficient between each pair of continuous variables in the 'cont.cor.coef'
#' for the unwanted variables, followed by the minimum cut-off of the correlation coefficient for the biological variables.
#' If the correlation of a pair of variable is higher than the minimum cut-off, only the variable that has the highest variance
#' will be kept and the other one will be excluded from the remaining analysis.
#'

#' @return a SummarizedExperiment object and the selected 'uv.variables' and 'bio.variables'.

#' @author Ramyar Molania

#' @importFrom SummarizedExperiment assay colData
#' @importFrom DescTools ContCoef
#' @export

variablesCorrelation <- function(
        se.obj,
        bio.variables,
        uv.variables,
        cat.cor.coef = c(0.7, 0.7),
        cont.cor.coef = c(0.7, 0.7),
        assess.se.obj = TRUE,
        remove.na = 'sample.annotation',
        verbose = TRUE) {
    printColoredMessage(message = '------------The variablesCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if (is.null(bio.variables) & is.null(uv.variables)) {
        stop('Both "bio.variables" and "uv.variables" cannot be empty.')
    }
    if (max(cat.cor.coef) > 1) {
        stop('The maximum value of "cat.cor.coef" cannot be more than 1.')
    } else if (max(cont.cor.coef) > 1) {
        stop('The maximum value of "cont.cor.coef" cannot be more than 1.')
    } else if (length(cont.cor.coef) == 1 | length(cat.cor.coef) == 1) {
        stop('Please provide two correlations in "cat.cor.coef" and "cont.cor.coef" for both unwanted and biological variables.')
    } else if (length(intersect(bio.variables, uv.variables)) != 0) {
        if (length(intersect(bio.variables, uv.variables)) == 1) {
            a <- 'a common variable'
        } else a <- 'common variables'
        stop(paste0(
                'The "bio.variables" and "uv.variables" have ',
                paste0(intersect(bio.variables, uv.variables), collapse = ' &'),
                a,
                '. Each variable is either biological and defined in the "bio.variables" ',
                'or unwanted variation and defined in the "uv.variables".'))
    }
    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = NULL,
            variables = c(bio.variables, uv.variables),
            remove.na = remove.na,
            verbose = verbose)
    }
    # unwanted variables ####
    printColoredMessage(message = '--Assess the unwanted variation variables:',
                        color = 'magenta',
                        verbose = verbose)
    if (length(uv.variables) > 0) {
        ## find classes of unwanted variation variables ####
        uv.var.class <- unlist(lapply(
            uv.variables,
            function(x) class(colData(se.obj)[[x]])))
        categorical.uv <- uv.variables[uv.var.class %in% c('factor', 'character')]
        continuous.uv <- uv.variables[uv.var.class %in% c('numeric', 'integer')]
        if (length(categorical.uv) > 0) {
            if (length(categorical.uv) == 1) {
                a <- 'is a categorical variable'
            } else a <- 'are categorical variables'
            printColoredMessage(
                message = paste0(
                    'The ',
                    paste0(categorical.uv, collapse = ' & '),
                    a,
                    ' of unwanted variation.'),
                color = 'blue',
                verbose = verbose)
        }
        if (length(continuous.uv) > 0) {
            if (length(categorical.uv) == 1) {
                a <- 'is a continuous variable'
            } else a <- 'are continuous variables'
            printColoredMessage(
                message = paste0(
                    'The ',
                    paste0(categorical.uv, collapse = ' & '),
                    a,
                    ' of unwanted variation.'),
                color = 'blue',
                verbose = verbose)
        }
        # variation in unwanted variables ####
        printColoredMessage(
            message = '-Checke the levels and variance of categorical and continuous unwanted variation variables, respectively:',
            color = 'magenta',
            verbose = verbose)
        check.uv.vars <- lapply(
            uv.variables,
            function(x) {
                class.type <- class(colData(se.obj)[[x]])
                if (class.type %in% c('factor', 'character')) {
                    if (length(unique(colData(se.obj)[[x]])) == 1) {
                        stop(paste0(
                                'Each categorical unwanted variables "uv.variables" must contain at least two groups.',
                                'However the variable ',
                                x,
                                ' contains only one group:',
                                unique(colData(se.obj)[[x]]),
                                'Please remove this variable from "uv.variables" argument and re-run the function.'))
                    } else
                        printColoredMessage(message = paste0('The ', x, ' contains ', length(unique(colData(se.obj)[[x]])), ' groups.'),
                            color = 'blue',
                            verbose = verbose)
                }
                if (class.type %in% c('numeric', 'integer')) {
                    if (var(colData(se.obj)[[x]]) == 0) {
                        stop(paste0('Each continuous unwanted variables "uv.variables" must contain some variation.',
                                'However the variable ',
                                x,
                                ' contains only one group:',
                                unique(colData(se.obj)[[x]]),
                                'Please remove this variable from "uv.variables" argument and re-run the function.'
                            )
                        )
                    } else
                        printColoredMessage(
                            message = paste0( 'The variance of ', x,' is ', round(var(colData(se.obj)[[x]]), digits = 4), '.'),
                            color = 'blue',
                            verbose = verbose)
                }
            })
        # assess the correlation between categorical UV ####
        printColoredMessage(message = '--Assess the correlation between categoricals variables:',
                            color = 'magenta',
                            verbose = verbose)
        if (length(categorical.uv) > 1) {
            all.pairs <- combn(categorical.uv , 2)
            remove.cat.uv.variable <- lapply(
                1:ncol(all.pairs),
                function(x) {
                    cat.cor <- ContCoef(x = se.obj[[all.pairs[, x][1]]], y = se.obj[[all.pairs[, x][2]]])
                    cat.cor <- round(x = cat.cor, digits = 3)
                    if (cat.cor > cat.cor.coef[1]) {
                        printColoredMessage(
                            paste0(
                                'The variables ',
                                all.pairs[, x][1],
                                ' and ',
                                all.pairs[, x][2],
                                ' are highly associated (corr.coef: ~',
                                cat.cor,
                                '). The one with the higher number of factors will be selected.'),
                            color = 'blue',
                            verbose = verbose
                        )
                        variable.freq <- c(
                            length(unique(se.obj[[all.pairs[, x][1]]])) ,
                            length(unique(se.obj[[all.pairs[, x][2]]]))
                            )
                        if (diff(variable.freq) == 0) {
                            remove.variables <- all.pairs[1, x]
                        } else remove.variables <- all.pairs[, x][which(variable.freq != max(variable.freq))]
                    } else {
                        printColoredMessage(
                            paste0(
                                'The variables ',
                                all.pairs[, x][1],
                                ' and ',
                                all.pairs[, x][2],
                                ' are not highly associated (corr.coef:',
                                cat.cor,
                                ').'),
                            color = 'blue',
                            verbose = verbose )
                        remove.variables <- NULL
                    }
                    return(remove.variables)
                })
            remove.cat.uv.variable <- unique(unlist(remove.cat.uv.variable))
            categorical.uv <- categorical.uv[!categorical.uv %in% remove.cat.uv.variable]
            if (length(categorical.uv) == 1) {
                a <- 'is selected as source'
            } else a <- 'are selected as sources'
            printColoredMessage(
                message = paste0(
                        'The ',
                        paste0(categorical.uv, collapse = ' & '),
                        a,
                        ' of categorical unwanted variation variables.'),
                color = 'blue',
                verbose = verbose)
        } else {
            printColoredMessage(message = 'There is only one source of categorical variable.')
        }
        # assess the correlation between continuous UV ####
        if (length(continuous.uv) > 1) {
            printColoredMessage(message = '--Assess the correlation between categoricals variables:',
                                color = 'magenta',
                                verbose = verbose)
            all.pairs <- combn(continuous.uv , 2)
            remove.cont.uv.variable <- lapply(
                1:ncol(all.pairs),
                function(x) {
                    corr.coef <- suppressWarnings(cor.test(
                        x = se.obj[[all.pairs[1 , x]]],
                        y = se.obj[[all.pairs[2 , x]]],
                        method = 'spearman'
                    ))[[4]]
                    corr.coef <- round(x = abs(corr.coef), digits = 3)
                    if (corr.coef > cont.cor.coef[1]) {
                        printColoredMessage(
                            message =
                                paste0(
                                    'The variables ',
                                    all.pairs[, x][1],
                                    ' and ',
                                    all.pairs[, x][2],
                                    ' are highly correlated (spearman.corr.coef:',
                                    corr.coef,
                                    '). The one with the higest variance will be selected.'),
                            color = 'blue',
                            verbose = verbose)
                        variable.freq <- c(var(se.obj[[all.pairs[, x][1]]]), var(se.obj[[all.pairs[, x][2]]]))
                        if (diff(variable.freq) == 0) {
                            remove.variables <- all.pairs[1, x]
                        } else remove.variables <- all.pairs[, x][which(variable.freq != max(variable.freq))]
                    } else {
                        printColoredMessage(
                            message = gsub(
                                '"',
                                '',
                                paste0(
                                    'The variables ',
                                    all.pairs[, x][1],
                                    ' and ',
                                    all.pairs[, x][2],
                                    ' are not highly correlated (spearman.corr.coef:',
                                    corr.coef,
                                    ').')),
                            color = 'blue',
                            verbose = verbose
                        )
                        remove.variables <- NULL
                    }
                    return(remove.variables)
                })
            remove.cont.uv.variable <- unique(unlist(remove.cont.uv.variable))
            continuous.uv <- continuous.uv[!continuous.uv %in% remove.cont.uv.variable]
            if (length(continuous.uv) == 1) {
                a <- 'is selected as continuous source'
            } else a <- 'are selected as continuous sources'
            printColoredMessage(
                message =
                    paste0(
                        'The ',
                        paste0(continuous.uv, collapse = ' & '),
                        a,
                        ' of unwanted variation.'),
                color = 'blue',
                verbose = verbose
            )
        } else
            printColoredMessage(
                message = 'There is only one source of continuous variable.',
                color = 'blue',
                verbose = verbose)
    } else
        printColoredMessage(
            message = 'Any unwanted variables are provided.',
            color = 'blue',
            verbose = verbose)

    # biological variables ####
    printColoredMessage(message = '--Assess the biological variables:',
                        color = 'magenta',
                        verbose = verbose)
    if (length(bio.variables) > 0) {
        ## classes of biological variables and variation
        bio.var.class <- unlist(lapply(
            bio.variables,
            function(x) class(colData(se.obj)[[x]])))
        categorical.bio <- bio.variables[bio.var.class %in% c('factor', 'character')]
        continuous.bio <- bio.variables[bio.var.class %in% c('numeric', 'integer')]
        if (length(categorical.bio) > 0) {
            if (length(categorical.bio) == 1) {
                a <- 'is a categorical variable'
            } else a <- 'are categorical variables'
            printColoredMessage(
                message = paste0(
                    'The ',
                    paste0(categorical.bio, collapse = ' & '),
                    a,
                    ' of biological variation.'),
                color = 'blue',
                verbose = verbose)
        }
        if (length(continuous.bio) > 0) {
            if (length(continuous.bio) == 1) {
                a <- 'is a continuous source '
            } else a <- 'are continuous sources '
            printColoredMessage(
                message = paste0(
                    'The ',
                    paste0(continuous.bio, collapse = ' & '),
                    a,
                    ' biological variation.'),
                color = 'blue',
                verbose = verbose)
        }
        ## variation in biological variables ####
        printColoredMessage(message = '--Check the levels and variance of categorical and continuous biological variables, respectively:',
                            color = 'magenta',
                            verbose = verbose)
        check.bio.vars <- lapply(
            bio.variables,
            function(x) {
                class.type <- class(colData(se.obj)[[x]])
                if (class.type %in% c('factor', 'character')) {
                    if (length(unique(colData(se.obj)[[x]])) == 1) {
                        stop(paste0(
                                'biological variables must contain at least two groups.',
                                'The biological variable ',
                                x,
                                ' contains only one group:',
                                unique(colData(se.obj)[[x]]),
                                'Please remove this variable from "bio.variables" argument and re-run the function.'))
                    } else {
                        printColoredMessage(
                            message = paste0('The ',
                                             x,
                                             ' contains ',
                                             length(unique(
                                                 colData(se.obj)[[x]]
                                             )),
                                             ' groups.'),
                            color = 'blue',
                            verbose = verbose)
                    }
                }
                if (class.type %in% c('numeric', 'integer')) {
                    if (var(colData(se.obj)[[x]]) == 0) {
                        stop( paste0(
                                'The variance of biological variable ',
                                x,
                                ' is equal to 0.',
                                'This variable must contain some variation.',
                                'Please remove this variable from "bio.variables" argument and re-run the function.'))
                    } else
                        printColoredMessage(
                            message = paste0(
                                'The variance of ',
                                x,
                                ' is ',
                                round(var(colData(se.obj)[[x]]), digits = 4),
                                '.'),
                            color = 'blue',
                            verbose = verbose)
                }
            })
        # assess correlation between categorical sources of biological variables ####
        if (length(categorical.bio) > 1) {
            printColoredMessage(message = '--Assess the correlation between categoricals variables:',
                                color = 'magenta',
                                verbose = verbose)
            all.pairs <- combn(categorical.bio , 2)
            remove.cat.bio.variable <- lapply(
                1:ncol(all.pairs),
                function(x) {
                    cat.cor <- ContCoef(x = se.obj[[all.pairs[, x][1]]], y = se.obj[[all.pairs[, x][2]]])
                    cat.cor <- round(x = cat.cor, digits = 3)
                    if (cat.cor > cat.cor.coef[2]) {
                        printColoredMessage(
                            paste0(
                                'The variables ',
                                all.pairs[, x][1],
                                ' and ',
                                all.pairs[, x][2],
                                ' are highly associated (corr.coef: ~',
                                cat.cor,
                                '). The one with the higher number of factors will be selected.'),
                            color = 'blue',
                            verbose = verbose)
                        variable.freq <- c(
                            length(unique(se.obj[[all.pairs[, x][1]]])),
                            length(unique(se.obj[[all.pairs[, x][2]]]))
                            )
                        if (diff(variable.freq) == 0) {
                            remove.variables <- all.pairs[1, x]
                        } else remove.variables <- all.pairs[, x][which(variable.freq != max(variable.freq))]
                    } else {
                        printColoredMessage(
                            paste0(
                                'The variables ',
                                all.pairs[, x][1],
                                ' and ',
                                all.pairs[, x][2],
                                ' are not highly associated (corr.coef:',
                                cat.cor,
                                ').'),
                            color = 'blue',
                            verbose = verbose
                        )
                        remove.variables <- NULL
                    }
                    return(remove.variables)
                })
            remove.cat.bio.variable <- unique(unlist(remove.cat.bio.variable))
            categorical.bio <- categorical.bio[!categorical.bio %in% remove.cat.bio.variable]
            if (length(categorical.bio) == 1) {
                a <- 'is selected as a source'
            } else a <- 'are selected as a sources'
            printColoredMessage(
                message =
                    paste0(
                        'The ',
                        paste0(categorical.bio, collapse = ' & '),
                        a,
                        ' of categorical biological variation.'),
                color = 'blue',
                verbose = verbose)
        }

        # assess the correlation between continuous sources of biological variation ####
        if (length(continuous.bio) > 1) {
            printColoredMessage(message = '--Assess the correlation between continuous variables:',
                                color = 'magenta',
                                verbose = verbose)
            all.pairs <- combn(continuous.bio , 2)
            remove.cont.bio.variable <- lapply(
                1:ncol(all.pairs),
                function(x) {
                    corr.coef <-
                        suppressWarnings(cor.test(
                            x = se.obj[[all.pairs[1 , x]]],
                            y = se.obj[[all.pairs[2 , x]]],
                            method = 'spearman'
                        ))[[4]]
                    corr.coef <- round(x = abs(corr.coef), digits = 3)
                    if (corr.coef > cont.cor.coef[2]) {
                        printColoredMessage(
                            message =
                                paste0(
                                    'The variables ',
                                    all.pairs[, x][1],
                                    ' and ',
                                    all.pairs[, x][2],
                                    ' are highly correlated (spearman.corr.coef:',
                                    corr.coef,
                                    '). The one with the highest variance will be selected.'),
                            color = 'blue',
                            verbose = verbose
                        )
                        variable.freq <- c(var(se.obj[[all.pairs[, x][1]]]), var(se.obj[[all.pairs[, x][2]]]))
                        if (diff(variable.freq) == 0) {
                            remove.variables <- all.pairs[1, x]
                        } else remove.variables <- all.pairs[, x][which(variable.freq != max(variable.freq))]
                    } else{
                        printColoredMessage(
                            message = gsub(
                                '"',
                                '',
                                paste0(
                                    'The variables ',
                                    all.pairs[, x][1],
                                    ' and ',
                                    all.pairs[, x][2],
                                    ' are not highly correlated (spearman.corr.coef:',
                                    corr.coef,
                                    ').')),
                            color = 'blue',
                            verbose = verbose
                        )
                        remove.variables <- NULL
                    }
                    return(remove.variables)
                })
            remove.cont.bio.variable <- unique(unlist(remove.cont.bio.variable))
            continuous.bio <- continuous.bio[!continuous.bio %in% remove.cont.bio.variable]
            if (length(continuous.bio) == 1) {
                a <- 'as continuous source'
            } else a <- 'as continuous source'
            printColoredMessage(
                message =
                    paste0(
                        'The ',
                        paste0(continuous.bio, collapse = ' & '),
                        a,
                        'of biological variation.'
                    ),
                color = 'blue',
                verbose = verbose)

        }
    }
    # save the results ####
    if (length(uv.variables) != 0 & length(bio.variables) != 0) {
        all.res <- list(
            se.obj = se.obj,
            uv.variables = c(continuous.uv, categorical.uv),
            bio.variables = c(continuous.bio, categorical.bio)
        )
    } else if (length(uv.variables) == 0 &
               length(bio.variables) != 0) {
        all.res <-  list(se.obj = se.obj,
                         bio.variables = c(continuous.bio, categorical.bio))
    } else if (length(uv.variables) != 0 &
               length(bio.variables) == 0) {
        all.res <- list(se.obj = se.obj,
                        uv.variables = c(continuous.uv, categorical.uv))
    }

    printColoredMessage(message = '------------The variablesCorrelation function finished.',
                        color = 'white',
                        verbose = verbose)
    return(all.res)
}
