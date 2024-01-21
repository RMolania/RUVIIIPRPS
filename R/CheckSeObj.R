#' is used to check the assay names, variables and missing/NA values in a SummarizedExperiment object.

#' @author Ramyar Molania

#' @description
#' This functions assesses the structure of a SummarizedExperiment object and removes any missing/NA values from both
#' assay(s) and sample annotation. If there are missing/NA in only of the assays, the corresponding rows in other assays
#' will be remove as well. Please note that, the current RUV-III-PRPS method do not support missing/NA values in the assay(s).

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. Indicates the name(s) of the assay(s) in the SummarizedExperiment object. By default it is
#' set to 'all', then all the assay(s) will be assessed.
#' @param variables Symbol. Indicates the name(s) of the column(s) in the sample annotation of the SummarizedExperiment
#' object. By default it is set to 'all', then all the columns will be checked. We recommend to specify those column(s)
#' that are of interest to your analysis.
#' @param remove.na Symbol. Indicates whether to remove missing/NA values from either the 'assays', 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains missing/NA values will be excluded. If 'sample.annotation'
#' is selected, the samples that contains NA or missing values for each 'variables' will be excluded. By default, it is
#' set to 'both'.
#' @param verbose Logical. If TRUE shows the process messages.

#' @return A SummarizedExperiment object.

#' @importFrom SummarizedExperiment assays assay colData
#' @importFrom stats complete.cases
#' @import ggplot2
#' @export

checkSeObj <- function(
        se.obj,
        assay.names = 'all',
        variables = 'all',
        remove.na = 'both',
        verbose = TRUE) {
    printColoredMessage(message = '------------The checkSeObj function starts:',
                        color = 'white',
                        verbose = verbose)
    # check the class and structure of the SummarizedExperiment object ####
    printColoredMessage(message = '-- Check the class and structure of the SummarizedExperiment object:',
                        color = 'magenta',
                        verbose = verbose)
    if (!class(se.obj)[1] == 'SummarizedExperiment') {
        stop('The "se.obj" provided is not a class of SummarizedExperiment object.')
    } else if (class(se.obj)[1] == 'SummarizedExperiment' & ncol(colData(se.obj)) == 0) {
        stop ('The sample annotation (colData(se.obj)) is not found in the SummarizedExperiment object.')
    } else {
        printColoredMessage(
            message = paste0(
                'The SummarizedExperiment object contains: ',
                ncol(se.obj),
                ' (samples or assays) and ',
                nrow(se.obj),
                ' (genes or measurements).'),
            color = 'blue',
            verbose = verbose)
    }
    # check function inputs ####
    if (!is.null(remove.na) & length(remove.na) > 1) {
        stop('The "remove.na" must be one of "assays", "sample.annotation", "both", or "none".')
    }
    if (!is.null(remove.na) &
        !remove.na %in% c("assays", "sample.annotation" , "both", "none")) {
        stop('The remove.na must be one of "assays", "sample.annotation", "both" or "none".')
    }
    if (remove.na == 'both') {
        if (is.null(assay.names) | is.null(variables))
            stop('The "assay.names" or "variables" cannot be empty when the "remove.na" is "both".')
    }
    if (remove.na == 'assays') {
        if (is.null(assay.names))
            stop('The "assay.names" cannot be empty when the "remove.na" is "assays".')
    }
    if (remove.na == 'sample.annotation') {
        if (is.null(variables))
            stop('The "variables" cannot be empty when the "remove.na" is "sample.annotation".')
    }
    # check the assays ####
    printColoredMessage(message = '-- Check the assay names of the SummarizedExperiment object',
                        color = 'magenta',
                        verbose = verbose)
    ## grammars ####
    assay.name <- 'assay'
    if (!is.null(assay.names)) {
        if (assay.names[1] == 'all') {
            if (length(names(assays(se.obj))) == 1) {
                assay.name = 'assay'
                assay.ps = 'is'
            } else{
                assay.name = 'assays'
                assay.ps = 'are'
            }
        } else if (length(assay.names) == 1) {
            assay.name = 'assay'
            assay.ps = 'is'
        } else if (length(assay.names) > 1) {
            assay.name = 'assays'
            assay.ps = 'are'
        }
    }
    variable.name <- 'variable'
    if (!is.null(variables)) {
        if (variables[1] == 'all') {
            if (length(colnames(colData(se.obj))) == 1) {
                variable.name = 'variable'
                variable.ps = 'is'
            } else{
                variable.name = 'variables'
                variable.ps = 'are'
            }
        } else if (length(variables) == 1) {
            variable.name = 'variable'
            variable.ps = 'is'
        } else if (length(variables) > 1) {
            variable.name = 'variables'
            variable.ps = 'are'
        }
    }
    if (is.null(assay.names)) {
        printColoredMessage(message = 'Please note: the names and missing/NA of all assays in the SummarizedExperiment object will not be checked.',
                            color = 'red',
                            verbose = verbose)
    } else if (assay.names[1] == 'all') {
        printColoredMessage(
            message = paste0(
                'The ',
                paste0(names(assays(se.obj)), collapse = ' & '),
                ' assays are found in the SummarizedExperiment object.'
            ),
            color = 'blue',
            verbose = verbose
        )
    } else if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)) {
        printColoredMessage(
            message = paste0(
                paste0(names(assays(se.obj)), collapse = ' & '),
                ' ',
                assay.ps,
                ' the current ',
                assay.name,
                ' of the SummarizedExperiment object.'
            ),
            color = 'blue',
            verbose = verbose
        )
        stop(if (length(assay.names[!assay.names %in% names(assays(se.obj))]) == 1) {
            paste0(
                'The assay name: ',
                paste0(assay.names[!assay.names %in% names(assays(se.obj))] , collapse = ' & '),
                ' is not found in the SummarizedExperiment object.'
            )
        } else {
            paste0(
                'The assay names: ',
                paste0(assay.names[!assay.names %in% names(assays(se.obj))] , collapse = ' & '),
                ' are not found in the SummarizedExperiment object.'
            )
        })

    } else if (sum(assay.names %in% names(assays(se.obj))) == length(assay.names)) {
        printColoredMessage(
            message = paste0(
                'The ',
                assay.name,
                ' ',
                paste0(assay.names, collapse = '&'),
                ' ',
                assay.ps,
                ' found in the SummarizedExperiment object.'
            ),
            color = 'blue',
            verbose = verbose
        )
    }
    # checking the variables ####
    printColoredMessage(message = '-- Check the variables in the SummarizedExperiment object:',
                        color = 'magenta',
                        verbose = verbose)
    if (is.null(variables)) {
        printColoredMessage(message = 'Please note: missing/NA of all variables in the sample annotation of the SummarizedExperiment object will not be checked.',
                            color = 'red',
                            verbose = verbose)
    } else if (variables[1] == 'all') {
        printColoredMessage(
            message = paste0(
                ncol(colData(se.obj)),
                ' ',
                variable.name,
                ' ',
                variable.ps,
                ' in the SummarizedExperiment object.'
            ),
            color = 'blue',
            verbose = verbose
        )
    } else if (sum(variables %in% colnames(colData(se.obj))) == length(variables)) {
        printColoredMessage(
            message = paste0(
                'The ',
                variable.name,
                ': ',
                paste0(variables, collapse = ' & '),
                ' ',
                variable.ps,
                ' found in the SummarizedExperiment object.'
            ),
            color = 'blue',
            verbose = verbose
        )
    } else {
        printColoredMessage(
            message = paste0(
                'The ',
                paste0(variables, collapse = ' & '),
                ': ',
                variable.name,
                ' ',
                variable.ps,
                ' provided as variables.'
            ),
            color = 'blue',
            verbose = verbose
        )
        if (!sum(variables %in% colnames(colData(se.obj))) == length(variables)) {
            stop(
                paste0(
                    'The ',
                    variable.name,
                    ': ',
                    paste0(variables[!variables %in% colnames(colData(se.obj))] , collapse = ' & '),
                    variable.ps,
                    'not found in the SummarizedExperiment object.'
                )
            )
        }
    }
    # remove missing/NA values ####
    printColoredMessage(message = '-- Remove missing/NA values:',
                        color = 'magenta',
                        verbose = verbose)
    if (!is.null(variables)) {
        if (variables[1] == 'all')
            variables <- colnames(colData(se.obj))
    }
    if (!is.null(assay.names)) {
        if (assay.names[1] == 'all')
            assay.names <- names(assays(se.obj))
    }
    if (remove.na != 'none') {
        if (remove.na == 'both') {
            printColoredMessage(
                message = 'Check the variables in the sample annotation:',
                color = 'blue',
                verbose = verbose
            )
            na.sample.annot <- apply(colData(se.obj)[, variables, drop = FALSE],
                                     2,
                                     function(x)
                                         sum(is.na(x)))
            na.sample.annot <- names(which(na.sample.annot > 0))
            if (length(na.sample.annot) > 0) {
                printColoredMessage(
                    message = paste0(
                        'The ',
                        paste0(na.sample.annot, collapse = ' &'),
                        ' ',
                        variable.name,
                        ' contains missing/NA values.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
                ncol.a <- ncol(se.obj)
                keep.samples <-
                    complete.cases(as.data.frame(colData(se.obj)[, variables, drop = FALSE]))
                se.obj <- se.obj[, keep.samples]
                ncol.b <- ncol(se.obj)
                if (c(ncol.a - ncol.b) == 1) {
                    ps = 'sample is excluded'
                } else
                    ps = 'samples are excluded'
                printColoredMessage(
                    message = paste0(
                        ncol.a - ncol.b,
                        ' ',
                        ps,
                        ' as a result of the presence of missing/NA values in the ',
                        variable.name,
                        '.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            } else{
                printColoredMessage(
                    message = paste0(
                        'Any missing/NA are found in the ',
                        paste0(variables, collapse = '&')
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            }
            printColoredMessage(
                message = 'Check the assays:',
                color = 'blue',
                verbose = verbose
            )
            na.measurements <- sapply(assay.names,
                                      function(x)
                                          rowSums(is.na(assay(
                                              se.obj, x
                                          ))) > 0)
            if (sum(na.measurements) > 0) {
                nrow.a <- nrow(se.obj)
                se.obj <- se.obj[!na.measurements, ]
                nrow.b <- nrow(se.obj)
                if (c(nrow.a - nrow.b) == 1) {
                    ps <- 'measurement is excluded'
                } else
                    ps <- 'measurements are excluded'
                printColoredMessage(
                    message = paste0(
                        nrow.a - nrow.b,
                        ' ',
                        ps,
                        '  from the ',
                        assay.name,
                        ' as a result of the presence of missing/NA values.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            } else{
                printColoredMessage(
                    message = paste0(
                        'Any missing/NA are found in the ',
                        paste0(assay.names, collapse = '&'),
                        assay.name,
                        '.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            }
        } else if (remove.na == 'sample.annotation') {
            printColoredMessage(
                message = 'Check the variables in sample annotation:',
                color = 'blue',
                verbose = verbose
            )
            na.sample.annot <- apply(colData(se.obj)[, variables, drop = FALSE],
                                     2,
                                     function(x)
                                         sum(is.na(x)))
            na.data <- names(which(na.sample.annot > 0))
            if (length(na.data) > 0) {
                printColoredMessage(
                    message = paste0(
                        paste0(na.data, collapse = ' &'),
                        ' contains missing/NA values.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
                ncol.a <- ncol(se.obj)
                keep.samples <-
                    complete.cases(as.data.frame(colData(se.obj)[, variables, drop = FALSE]))
                se.obj <- se.obj[, keep.samples]
                ncol.b <- ncol(se.obj)
                if (c(ncol.a - ncol.b) == 1) {
                    ps <- 'sample is excluded'
                } else
                    ps <- 'samples are excluded'
                printColoredMessage(
                    message = paste0(
                        ncol.a - ncol.b,
                        ' ',
                        ps,
                        ' as a result of the presence of missing/NA values in the ',
                        variable.name,
                        '.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            } else{
                printColoredMessage(
                    message = paste0(
                        'Any missing/NA are found in ',
                        paste0(variables, collapse = '&'),
                        ':'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            }
            printColoredMessage(
                message = 'Note, the current RUV-III-PRPS method do not support NA in the measurements.',
                color = 'red',
                verbose = verbose
            )
        } else if (remove.na == 'assays') {
            printColoredMessage(
                message = 'Check the assays:',
                color = 'blue',
                verbose = verbose
            )
            na.measurements <- sapply(assay.names,
                                      function(x)
                                          rowSums(is.na(assay(
                                              se.obj, x
                                          ))) > 0)
            na.measurements <- rowSums(na.measurements) > 0
            if (sum(na.measurements) > 0) {
                nrow.a <- nrow(se.obj)
                se.obj <- se.obj[!na.measurements, ]
                if (c(nrow.a - nrow.b) == 1) {
                    ps <- 'measurement is excluded'
                } else
                    ps <- 'measurements are excluded'
                printColoredMessage(
                    message = paste0(
                        nrow.a - nrow.b,
                        ps,
                        '  from the ',
                        assay.name,
                        ' as a result of the presence of missing/NA values.'
                    ),
                    color = 'red',
                    verbose = verbose
                )
            } else{
                printColoredMessage(
                    message = paste0(
                        'Any missing/NA values are found in the ',
                        assay.name,
                        '.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            }
        }
        else{
            printColoredMessage(
                message = 'Missing/NA values in both assays and variables are not checked.',
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = 'Note, the current RUV-III-PRPS method do not support missing/NA in the measurements.',
                color = 'red',
                verbose = verbose
            )
        }
    }
    printColoredMessage(message = '------------The checkSeObj function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}
