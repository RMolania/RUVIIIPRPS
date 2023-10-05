#' is used to assess the association between variables of a SummarizedExperiment object using Spearman correlation.
#'
#' @param se.obj A summarized experiment object.
#' @param assay.name String for the selection of the name of the assay of the SummarizedExperiment class object used to define PRPS.
#' @param bio.variables String of the label of (a) categorical variable(s) that specifies major biological groups
#' such as samples types from colData(se).
#' @param uv.variables String or vector of strings of the label of continuous or categorical variable(s)
#' such as samples types, batch or library size from colData(se) that will be used to define PRPS.
#' @param cont.coef TO BE DEFINED.
#' @param spearman.coef Vector of two numeric values. Each value is the cut-off to consider there is an association
#' between the variables based on the Spearman correlation coefficient. The first one is used for the 'uv.variables' and the second for 'bio.variables'.
#' @param assess.se.obj Logical. Whether to assess the SummarizedExperiment object or not.
#' @param remove.na To remove NA or missing values from either the assays or sample annotation or both.
#' @param verbose Logical. Whether to show the messages of the functions or not.

#' @return a SummarizedExperiment object and the selected 'uv.variables' and 'bio.variables'.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom DescTools ContCoef
#' @importFrom stats complete.cases
#' @export

variablesCorrelation <- function(
#checkVariables <- function(
        se.obj,
        assay.name = NULL,
        bio.variables,
        uv.variables,
        cont.coef = c(0.7, 0.7),
        spearman.coef = c(0.7, 0.7),
        assess.se.obj = TRUE,
        remove.na = 'sample.annotation',
        verbose = TRUE) {
    ## checking arguments
    printColoredMessage(message = '------------The checkVariables function starts:',
                          color = 'white',
                          verbose = verbose)
    printColoredMessage(message = '### Checking the arguments inputs of the function:',
                          color = 'magenta',
                          verbose = verbose)
    if (is.null(bio.variables) & is.null(uv.variables)) {
        stop('Both bio.variables and uv.variables are empty, please provide at least one of them')
    } else if (max(cont.coef) > 1) {
        stop('The cont.coef argument cannot be more than 1.')
    } else if (max(spearman.coef) > 1) {
        stop('The spearman.coef argument cannot be more than 1.')
    } else if (length(cont.coef) == 1 | length(cont.coef) == 1) {
        stop( 'Please provide correlation coefs for both unwanted and biological variatiion.')
    } else if (length(intersect(bio.variables, uv.variables)) != 0) {
        stop(if (length(intersect(bio.variables, uv.variables)) == 1) {
            paste0(
                'The bio.variables and uv.variables have ',
                paste0(
                    intersect(bio.variables, uv.variables),
                    collapse = ' &'
                ),
                ' as a common variable.',
                ' A variable is either biological and defined in the "bio.variables" or unwanted variation and defined in the "uv.variables".'
            )
        } else{
            paste0(
                'The bio.variables and uv.variables have ',
                paste0(intersect(bio.variables, uv.variables),collapse = ' &'
                ),
                ' as common variable/s.',
                ' A variable is either biological and defined in the "bio.variables" or unwanted variation and defined in the "uv.variables"'
            )
        })
    }

    ### Checking summarized experiment objects
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(bio.variables, uv.variables),
            remove.na = remove.na,
            verbose = verbose
        )
    }
    ### Checking unwanted variation variables
    if (length(uv.variables) > 0) {
        ## classes of uv variables and variation
        printColoredMessage(
            message = '### Checking the unwanted variation variables "uv.variables":',
            color = 'magenta',
            verbose = verbose
        )
        uv.var.class <- unlist(lapply(uv.variables, function(x)
            class(colData(se.obj)[[x]])))
        categorical.uv <-
            uv.variables[uv.var.class %in% c('factor', 'character')]
        continuous.uv <-
            uv.variables[uv.var.class %in% c('numeric', 'integer')]
        if (length(categorical.uv) > 0) {
            if (length(categorical.uv) == 1) {
                printColoredMessage(
                    message = paste0(
                        'The ',
                        paste0(categorical.uv, collapse = ' & '),
                        ' is a categorical variable of unwanted variation.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            } else {
                printColoredMessage(
                    message = paste0(
                        'The ',
                        paste0(categorical.uv, collapse = ' & '),
                        ' are categorical variables of unwanted variation.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            }
        } else {
            printColoredMessage(
                message = 'No categorical variable of unwanted variation is provided.',
                color = 'red',
                verbose = verbose
            )
        }
        if (length(continuous.uv) > 0) {
            if (length(continuous.uv) == 1) {
                printColoredMessage(
                    message = paste0(
                        paste0('The ',
                               continuous.uv, collapse = ' & '),
                        ' is a continuous variable of unwanted variation.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            } else {
                printColoredMessage(
                    message = paste0(
                        paste0('The ',
                               continuous.uv, collapse = ' & '),
                        ' are continuous variables of unwanted variation.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            }

        } else {
            printColoredMessage(
                message = 'No continuous variable of unwanted variation were provided.',
                color = 'red',
                verbose = verbose
            )
        }
        ###
        uv.variables <- c(continuous.uv, categorical.uv)
        printColoredMessage(
            message = '### Checking the levels and variance of categorical and continuous unwanted variables, respectively:',
            color = 'magenta',
            verbose = verbose
        )
        check.uv.vars <- lapply(uv.variables,
                                function(x) {
                                    class.type <- class(colData(se.obj)[[x]])
                                    if (class.type %in% c('factor', 'character')) {
                                        if (length(unique(colData(se.obj)[[x]])) == 1) {
                                            stop(
                                                paste0(
                                                    'The Unwanted variable "uv.variables" must contain at least two groups.',
                                                    'However the variable ',
                                                    x,
                                                    ' contains only one group:',
                                                    unique(colData(se.obj)[[x]]),
                                                    'Please remove this variable from "uv.variables" argument and re-run the function.'
                                                )
                                            )
                                        } else {
                                            printColoredMessage(
                                                message = paste0(
                                                    'The ',
                                                    x,
                                                    ' contains ',
                                                    length(unique(colData(
                                                        se.obj
                                                    )[[x]])),
                                                    ' groups.'
                                                ),
                                                color = 'blue',
                                                verbose = verbose
                                            )
                                        }
                                    }
                                    if (class.type %in% c('numeric', 'integer')) {
                                        if (var(colData(se.obj)[[x]]) == 0) {
                                            stop(
                                                paste0(
                                                    'The variance of the unwanted variable ',
                                                    x,
                                                    ' is equal to 0.',
                                                    'However, this variable must contain some variation.',
                                                    'Please remove this variable from "uv.variables" argument and re-run the function.'
                                                )
                                            )
                                        } else{
                                            printColoredMessage(
                                                message = paste0(
                                                    'The variance of ',
                                                    x,
                                                    ' is ',
                                                    round(var(colData(
                                                        se.obj
                                                    )[[x]]), digits = 4),
                                                    '.'
                                                ),
                                                color = 'blue',
                                                verbose = verbose
                                            )
                                        }
                                    }
                                })
        ### Check correlation between categorical sources of unwanted variation
        if (length(categorical.uv) > 1) {
            printColoredMessage(
                message = '### As several categorical variables of unwanted variation were provided, their respective associations
                will be assessed.'
                ,
                color = 'magenta',
                verbose = verbose
            )
            printColoredMessage(
                message = 'Applying chi-square test to assess the association between all pairs of categorical unwanted variables',
                color = 'green',
                verbose = verbose
            )
            all.pairs <- combn(categorical.uv , 2)
            remove.cat.uv.variable <- lapply(1:ncol(all.pairs),
                                             function(x) {
                                                 cont.coef <- ContCoef(x = se.obj[[all.pairs[, x][1]]],
                                                                                   y = se.obj[[all.pairs[, x][2]]])
                                                 cont.coef <-
                                                     round(x = cont.coef, digits = 3)
                                                 if (cont.coef > cont.coef[1]) {
                                                     printColoredMessage(
                                                         paste0(
                                                             'The variables ',
                                                             all.pairs[, x][1],
                                                             ' and ',
                                                             all.pairs[, x][2],
                                                             ' are highly associated (corr.coef: ~',
                                                             cont.coef,
                                                             '). The one with higher groups number for PRPS will be selected.'
                                                         ),
                                                         color = 'blue',
                                                         verbose = verbose
                                                     )
                                                     variable.freq <-
                                                         c(length(unique(se.obj[[all.pairs[, x][1]]])) ,
                                                           length(unique(se.obj[[all.pairs[, x][2]]])))
                                                     if (diff(variable.freq) == 0) {
                                                         remove.variables <- all.pairs[1, x]
                                                     } else{
                                                         remove.variables <-
                                                             all.pairs[, x][which(variable.freq != max(variable.freq))]
                                                     }
                                                 } else{
                                                     printColoredMessage(
                                                         paste0(
                                                             'The variables ',
                                                             all.pairs[, x][1],
                                                             ' and ',
                                                             all.pairs[, x][2],
                                                             ' are not highly associated (corr.coef:',
                                                             cont.coef,
                                                             ').'
                                                         ),
                                                         color = 'blue',
                                                         verbose = verbose
                                                     )
                                                     remove.variables <-
                                                         NULL
                                                 }
                                                 return(remove.variables)
                                             })
            remove.cat.uv.variable <-
                unique(unlist(remove.cat.uv.variable))
            categorical.uv <-
                categorical.uv[!categorical.uv %in% remove.cat.uv.variable]
            if (length(categorical.uv) == 1) {
                printColoredMessage(
                    message =
                        paste0(
                            'Finally, the variable ',
                            paste0(categorical.uv, collapse = ' & '),
                            ' is selected as source of categorical unwanted variation for creating PRPS.'
                        ),
                    color = 'blue',
                    verbose = verbose
                )
            } else {
                printColoredMessage(
                    message =
                        paste0(
                            'Finally, the variables ',
                            paste0(categorical.uv, collapse = ' & '),
                            ' are selected as sources of categorical unwanted variation for creating PRPS.'
                        ),
                    color = 'blue',
                    verbose = verbose
                )
            }
        } else if (length(categorical.uv) == 1) {
            printColoredMessage(
                message =
                    paste0(
                        'Finally, the variable ',
                        paste0(categorical.uv, collapse = ' & '),
                        ' is selected as source of categorical unwanted variation for creating PRPS.'
                    ),
                color = 'blue',
                verbose = verbose
            )

        }
        ### Check correlation between continuous sources of unwanted variation
        if (length(continuous.uv) > 1) {
            printColoredMessage(
                message = '### As several continuous sources of unwanted variation were provided,
                their respective associations based on correlation will be assessed.',
                color = 'magenta',
                verbose = verbose
            )
            printColoredMessage(
                message =
                    'Applying Spearman correlation test between all pairs of continuous sources of unwanted variation.',
                color = 'white',
                verbose = verbose
            )
            all.pairs <- combn(continuous.uv , 2)
            remove.cont.uv.variable <- lapply(1:ncol(all.pairs),
                                              function(x) {
                                                  corr.coef <-
                                                      suppressWarnings(cor.test(
                                                          x = se.obj[[all.pairs[1 , x]]],
                                                          y = se.obj[[all.pairs[2 , x]]],
                                                          method = 'spearman'
                                                      ))[[4]]
                                                  corr.coef <-
                                                      round(x = abs(corr.coef), digits = 3)
                                                  if (corr.coef > spearman.coef[1]) {
                                                      printColoredMessage(
                                                          message =
                                                              paste0(
                                                                  'The variables ',
                                                                  all.pairs[, x][1],
                                                                  ' and ',
                                                                  all.pairs[, x][2],
                                                                  ' are highly correlated (spearman.corr.coef:',
                                                                  corr.coef,
                                                                  '). The one with the higest variance will be selected for creating PRPS.'
                                                              ),
                                                          color = 'blue',
                                                          verbose = verbose
                                                      )
                                                      variable.freq <-
                                                          c(var(se.obj[[all.pairs[, x][1]]]), var(se.obj[[all.pairs[, x][2]]]))
                                                      if (diff(variable.freq) == 0) {
                                                          remove.variables <- all.pairs[1, x]
                                                      } else{
                                                          remove.variables <-
                                                              all.pairs[, x][which(variable.freq != max(variable.freq))]
                                                      }
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
                                                                  '). PRPS will be created for individual ones.'
                                                              )
                                                          ),
                                                          color = 'blue',
                                                          verbose = verbose
                                                      )
                                                      remove.variables <-
                                                          NULL
                                                  }
                                                  return(remove.variables)
                                              })
            remove.cont.uv.variable <-
                unique(unlist(remove.cont.uv.variable))
            continuous.uv <-
                continuous.uv[!continuous.uv %in% remove.cont.uv.variable]
            if (length(continuous.uv) == 1) {
                printColoredMessage(
                    message =
                        paste0(
                            'Finally, the variable ',
                            paste0(continuous.uv, collapse = ' & '),
                            ' is selected as continuous source of unwanted variation for creating PRPS.'
                        ),
                    color = 'blue',
                    verbose = verbose
                )

            } else {
                printColoredMessage(
                    message =
                        paste0(
                            'Finally, the variables ',
                            paste0(continuous.uv, collapse = ' & '),
                            ' are selected as continuous sources of unwanted variation for creating PRPS.'
                        ),
                    color = 'blue',
                    verbose = verbose
                )
            }
        } else if (length(continuous.uv) == 1) {
            printColoredMessage(
                message =
                    paste0(
                        'Finally, the variable ',
                        paste0(continuous.uv, collapse = ' & '),
                        ' is selected as continuous source of unwanted variation for creating PRPS.'
                    ),
                color = 'blue',
                verbose = verbose
            )
        }
    }
    ### Checking biological variables
    if (length(bio.variables) > 0) {
        ## classes of uv variables and variation
        printColoredMessage(
            message = '### Checking the biological variables:',
            color = 'magenta',
            verbose = verbose
        )
        bio.var.class <- unlist(lapply(bio.variables, function(x)
            class(colData(se.obj)[[x]])))
        categorical.bio <-
            bio.variables[bio.var.class %in% c('factor', 'character')]
        continuous.bio <-
            bio.variables[bio.var.class %in% c('numeric', 'integer')]
        if (length(categorical.bio) > 0) {
            if (length(categorical.bio) == 1) {
                printColoredMessage(
                    message = paste0(
                        'The ',
                        paste0(categorical.bio, collapse = ' & '),
                        ' is a categorical variable of biological variation.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            } else {
                printColoredMessage(
                    message = paste0(
                        'The ',
                        paste0(categorical.bio, collapse = ' & '),
                        ' are categorical variables for the biological variation.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            }
        } else{
            printColoredMessage(
                message = 'No categorical variable of biological variation were provided.',
                color = 'red',
                verbose = verbose
            )
        }
        if (length(continuous.bio) > 0) {
            if (length(continuous.bio) == 1) {
                printColoredMessage(
                    message = paste0(
                        paste0('The ', continuous.bio, collapse = ' & '),
                        ' is a continuous source of biological variation.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            } else {
                printColoredMessage(
                    message = paste0(
                        paste0('The ', continuous.bio, collapse = ' & '),
                        ' are continuous sources of biological variation.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            }

        } else {
            printColoredMessage(
                message = 'No continuous variable of biological variation were provided.',
                color = 'red',
                verbose = verbose
            )
        }
        ###
        bio.variables <- c(continuous.bio, categorical.bio)
        printColoredMessage(
            message = '### Checking the levels and variance of categorical and continuous biological variables, respectively:',
            color = 'magenta',
            verbose = verbose
        )
        check.bio.vars <- lapply(bio.variables,
                                 function(x) {
                                     class.type <- class(colData(se.obj)[[x]])
                                     if (class.type %in% c('factor', 'character')) {
                                         if (length(unique(colData(se.obj)[[x]])) == 1) {
                                             stop(
                                                 paste0(
                                                     'biological variables must contain at least two groups.',
                                                     'The biological variable ',
                                                     x,
                                                     ' contains only one group:',
                                                     unique(colData(se.obj)[[x]]),
                                                     'Please remove this variable from "bio.variables" argument and re-run the function.'
                                                 )
                                             )
                                         } else {
                                             printColoredMessage(
                                                 message = paste0(
                                                     'The ',
                                                     x,
                                                     ' contains ',
                                                     length(unique(colData(
                                                         se.obj
                                                     )[[x]])),
                                                     ' groups.'
                                                 ),
                                                 color = 'blue',
                                                 verbose = verbose
                                             )
                                         }
                                     }
                                     if (class.type %in% c('numeric', 'integer')) {
                                         if (var(colData(se.obj)[[x]]) == 0) {
                                             stop(
                                                 paste0(
                                                     'The variance of biological variable ',
                                                     x,
                                                     ' is equal to 0.',
                                                     'This variable must contain some variation.',
                                                     'Please remove this variable from "bio.variables" argument and re-run the function.'
                                                 )
                                             )
                                         } else{
                                             printColoredMessage(
                                                 message = paste0(
                                                     'The variance of ',
                                                     x,
                                                     ' is ',
                                                     round(var(colData(
                                                         se.obj
                                                     )[[x]]), digits = 4),
                                                     '.'
                                                 ),
                                                 color = 'blue',
                                                 verbose = verbose
                                             )
                                         }
                                     }
                                 })
        ### Check correlation between categorical sources of unwanted variation
        if (length(categorical.bio) > 1) {
            printColoredMessage(
                message = '### Several categorical sources of biological variation were provided. The respective association between each pair of
                categorical sources of unwanted variation will be assessed.'
                ,
                color = 'magenta',
                verbose = verbose
            )

            printColoredMessage(
                message = 'Applying chi-square test between all pairs of categorical sources of biological variation.',
                color = 'green',
                verbose = verbose
            )
            all.pairs <- combn(categorical.bio , 2)
            remove.cat.bio.variable <- lapply(1:ncol(all.pairs),
                                              function(x) {
                                                  cont.coef <- ContCoef(x = se.obj[[all.pairs[, x][1]]],
                                                                                    y = se.obj[[all.pairs[, x][2]]])
                                                  cont.coef <-
                                                      round(x = cont.coef, digits = 3)
                                                  if (cont.coef > cont.coef[2]) {
                                                      printColoredMessage(
                                                          paste0(
                                                              'The variables ',
                                                              all.pairs[, x][1],
                                                              ' and ',
                                                              all.pairs[, x][2],
                                                              ' are highly associated (corr.coef: ~',
                                                              cont.coef,
                                                              '). The one with the higher groups number will selected for PRPS.'
                                                          ),
                                                          color = 'blue',
                                                          verbose = verbose
                                                      )
                                                      variable.freq <-
                                                          c(length(unique(se.obj[[all.pairs[, x][1]]])) ,
                                                            length(unique(se.obj[[all.pairs[, x][2]]])))
                                                      if (diff(variable.freq) == 0) {
                                                          remove.variables <- all.pairs[1, x]
                                                      } else{
                                                          remove.variables <-
                                                              all.pairs[, x][which(variable.freq != max(variable.freq))]
                                                      }
                                                  } else{
                                                      printColoredMessage(
                                                          paste0(
                                                              'The variables ',
                                                              all.pairs[, x][1],
                                                              ' and ',
                                                              all.pairs[, x][2],
                                                              ' are not highly associated (corr.coef:',
                                                              cont.coef,
                                                              ').'
                                                          ),
                                                          color = 'blue',
                                                          verbose = verbose
                                                      )
                                                      remove.variables <-
                                                          NULL
                                                  }
                                                  return(remove.variables)
                                              })
            remove.cat.bio.variable <-
                unique(unlist(remove.cat.bio.variable))
            categorical.bio <-
                categorical.bio[!categorical.bio %in% remove.cat.bio.variable]
            if (length(categorical.bio) == 1) {
                printColoredMessage(
                    message =
                        paste0(
                            'Finally, the variable ',
                            paste0(categorical.bio, collapse = ' & '),
                            ' is selected as a source of categorical biological variation for creating PRPS.'
                        ),
                    color = 'blue',
                    verbose = verbose
                )
            } else {
                printColoredMessage(
                    message =
                        paste0(
                            'Finally, the variables ',
                            paste0(categorical.bio, collapse = ' & '),
                            ' are selected as sources of categorical biological variation for creating PRPS.'
                        ),
                    color = 'blue',
                    verbose = verbose
                )
            }
        } else if (length(categorical.bio) == 1) {
            printColoredMessage(
                message =
                    paste0(
                        'Finally, the variable ',
                        paste0(categorical.bio, collapse = ' & '),
                        ' is selected as a source of categorical biological variation for creating PRPS.'
                    ),
                color = 'blue',
                verbose = verbose
            )

        }
        ### Check correlation between continuous sources of unwanted variation
        if (length(continuous.bio) > 1) {
            printColoredMessage(
                message = '### Several continuous variables of biological variation were provided.
                Checking the correlation between each pair of continuous variables of unwanted variation.',
                color = 'magenta',
                verbose = verbose
            )
            printColoredMessage(
                message =
                    'Applying Spearman correlation test between all possible pairs of continuous variables of biological variation.',
                color = 'white',
                verbose = verbose
            )
            all.pairs <- combn(continuous.bio , 2)
            remove.cont.bio.variable <- lapply(1:ncol(all.pairs),
                                               function(x) {
                                                   corr.coef <-
                                                       suppressWarnings(cor.test(
                                                           x = se.obj[[all.pairs[1 , x]]],
                                                           y = se.obj[[all.pairs[2 , x]]],
                                                           method = 'spearman'
                                                       ))[[4]]
                                                   corr.coef <-
                                                       round(x = abs(corr.coef), digits = 3)
                                                   if (corr.coef > spearman.coef[2]) {
                                                       printColoredMessage(
                                                           message =
                                                               paste0(
                                                                   'The variables ',
                                                                   all.pairs[, x][1],
                                                                   ' and ',
                                                                   all.pairs[, x][2],
                                                                   ' are highly correlated (spearman.corr.coef:',
                                                                   corr.coef,
                                                                   '). Then, one of them with the higest variance will be selected for creating PRPS.'
                                                               ),
                                                           color = 'blue',
                                                           verbose = verbose
                                                       )
                                                       variable.freq <-
                                                           c(var(se.obj[[all.pairs[, x][1]]]), var(se.obj[[all.pairs[, x][2]]]))
                                                       if (diff(variable.freq) == 0) {
                                                           remove.variables <- all.pairs[1, x]
                                                       } else{
                                                           remove.variables <-
                                                               all.pairs[, x][which(variable.freq != max(variable.freq))]
                                                       }
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
                                                                   '). PRPS will be created for individual ones.'
                                                               )
                                                           ),
                                                           color = 'blue',
                                                           verbose = verbose
                                                       )
                                                       remove.variables <-
                                                           NULL
                                                   }
                                                   return(remove.variables)
                                               })
            remove.cont.bio.variable <-
                unique(unlist(remove.cont.bio.variable))
            continuous.bio <-
                continuous.bio[!continuous.bio %in% remove.cont.bio.variable]

            if (length(continuous.bio) == 1) {
                printColoredMessage(
                    message =
                        paste0(
                            'Finally, the variable',
                            paste0(continuous.bio, collapse = ' & '),
                            ' is selected as continuous source of biological variation for creating PRPS.'
                        ),
                    color = 'blue',
                    verbose = verbose
                )

            } else {
                printColoredMessage(
                    message =
                        paste0(
                            'Finally, the variables ',
                            paste0(continuous.bio, collapse = ' & '),
                            ' are selected as continuous source of biological variation for creating PRPS.'
                        ),
                    color = 'blue',
                    verbose = verbose
                )
            }
        } else if (length(continuous.bio) == 1) {
            printColoredMessage(
                message =
                    paste0(
                        'Finally, the variable ',
                        paste0(continuous.bio, collapse = ' & '),
                        ' is selected as continuous source of biological variation for creating PRPS.'
                    ),
                color = 'blue',
                verbose = verbose
            )
        }
    }
    ### Checking biological variables
    printColoredMessage(message = '------------The checkVariables function finished.',
                          color = 'white',
                          verbose = verbose)
    #### output
    if (length(uv.variables) != 0 & length(bio.variables) != 0) {
        return(list(
            se.obj = se.obj,
            uv.variables = c(continuous.uv, categorical.uv),
            bio.variables = c(continuous.bio, categorical.bio)
        ))
    } else if (length(uv.variables) == 0 &
               length(bio.variables) != 0) {
        return(list(
            se.obj = se.obj,
            bio.variables = c(continuous.bio, categorical.bio)
        ))
    } else if (length(uv.variables) != 0 &
               length(bio.variables) == 0) {
        return(list(
            se.obj = se.obj,
            uv.variables = c(continuous.uv, categorical.uv)
        ))
    }
}
