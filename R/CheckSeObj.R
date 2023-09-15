#' is used to check summarized experiment object.
#'
#'
#' @param se.obj A summarized experiment object.
#' @param assay.names A name or names of the assays in summarized experiment object.
#' @param variables A name or names of the variables in the sample annotation of the summarized experiment object.
#' @param remove.na 'both', 'measurements', 'sample.annotation', 'none' TO BE DEFINED WHAT EACH DOES
#' @param verbose whether to show the messages or not.
#' @importFrom SummarizedExperiment assays assay colData
#' @importFrom stats complete.cases
#' @import ggplot2
#' @export

checkSeObj <- function(se.obj,
                       assay.names = 'All',
                       variables = 'All',
                       remove.na = 'both',
                       verbose = TRUE)
{
    ### Assess Summarized experiment object class
    printColoredMessage(message = '------------The checkSeObj function starts:',
                        color = 'white',
                        verbose = verbose)
    ### Checking the summarized experiment object
    printColoredMessage(message = '### Checking the SummarizedExperiment object:',
                        color = 'magenta',
                        verbose = verbose)
    if (!class(se.obj)[1] == 'SummarizedExperiment') {
        stop('The se.obj provided is not a SummarizedExperiment object.')
    } else if (class(se.obj)[1] == 'SummarizedExperiment' &
               ncol(colData(se.obj)) == 0) {
        stop('Sample annotation (colData(se.obj)) is not found in the SummarizedExperiment object.')
    } else {
        printColoredMessage(
            message = paste0(
                'The class of the se.obj is a SummarizedExperiment with: ',
                ncol(se.obj),
                ' (samples or assays) and ',
                nrow(se.obj),
                ' (genes or measurements).'
            ),
            color = 'blue',
            verbose = verbose
        )
    }
    ### grammars
    assay.name <- 'assay'
    if (!is.null(assay.names)) {
        if (assay.names[1] == 'All') {
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
        if (variables[1] == 'All') {
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
    ### Checking the assay name
    printColoredMessage(
        message = paste0(
            '### Checking the ',
            assay.name ,
            ' name in the SummarizedExperiment object:'
        ),
        color = 'magenta',
        verbose = verbose
    )
    if (is.null(assay.names)) {
        printColoredMessage(
            message = 'The assay.names = NULL so the names and NA/missing values of the assays of the SummarizedExperiment object will not be checked.',
            color = 'red',
            verbose = verbose)
    } else if (assay.names[1] == 'All') {
        printColoredMessage(
            message = paste0(
                paste0(names(assays(se.obj)), collapse = ' & '),
                ' ',
                assay.ps,
                ' the current ',
                assay.name,
                ' of the SummarizedExperiment object. NA or missing values will be removed from all the current ',
                assay.name,
                '.'
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
        } else{
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
    ### Checking the variables.
    printColoredMessage(
        message = paste0(
            '### Checking the ',
            variable.name,
            ' in the SummarizedExperiment object:'
        ),
        color = 'magenta',
        verbose = verbose
    )
    #### Check variables
    if (is.null(variables)) {
        printColoredMessage(
            message = 'The variables = NULL so the NA/missing values in the sample annotation of the SummarizedExperiment object will not be checked.',
            color = 'red',
            verbose = verbose)
    } else if (variables[1] == 'All') {
        printColoredMessage(
            message = 'The variables = All so the NA/missing values in all sample annotation of the SummarizedExperiment object will be removed.',
            color = 'red',
            verbose = verbose)
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
    ### Checking the remove na.
    if (!is.null(assay.names) | !is.null(variables)) {
        printColoredMessage(
            message = '### Checking the remove.na from the SummarizedExperiment object:',
            color = 'magenta',
            verbose = verbose
        )
        if (length(remove.na) > 1) {
            stop('The remove.na argument allows for only a single variable: "both", "measurements", "sample.annotation" or "none".')
        } else if (remove.na != 'none') {
            if (!remove.na %in% c('both',
                                  'measurements',
                                  'sample.annotation',
                                  'none')) {
                stop('The remove.na argument should correspond to one of the variables: "both", "measurements", "sample.annotation" or "none".')
            } else if (is.null(assay.names) & remove.na == 'both') {
                stop('Please provide at least an assay name to be able to find and remove NA or missing values.')
            } else if (is.null(assay.names) & remove.na == 'measurements') {
                stop('Please provide at least an assay name to be able to find and remove NA or missing values.')
            } else if (is.null(variables) & remove.na == 'both') {
                stop('Please provide at least a variable to be to able find and remove NA or missing values.')
            } else if (is.null(variables) & remove.na == 'sample.annotation') {
                stop('Please provide at least a variable to be to able find and remove NA or missing values.')
            } else {
                if (remove.na == 'both') {
                    printColoredMessage(
                        message = paste0(
                            'All NA or missing values found in both the ',
                            assay.name,
                            ' and ',
                            variable.name,
                            ' of the SummarizedExperiment will be removed.'
                        ),
                        color = 'blue',
                        verbose = verbose
                    )
                } else if (remove.na == 'measurements') {
                    printColoredMessage(
                        message = paste0(
                            'All NA or missing values found only in the ',
                            assay.name,
                            ' of the SummarizedExperiment will be removed.'
                        ),
                        color = 'blue',
                        verbose = verbose
                    )
                } else if (remove.na == 'sample.annotation') {
                    printColoredMessage(
                        message = paste0(
                            'All NA or missing values found  only in the ',
                            variable.name,
                            ' of the SummarizedExperiment will be removed.'
                        ),
                        color = 'blue',
                        verbose = verbose
                    )
                }
            }
        } else if (remove.na == 'none') {
            printColoredMessage(
                message = 'The argument remove.na = none, any NA or missing values will not be removed from the assay and from the sample annotation of the SummarizedExperiment object.',
                color = 'red',
                verbose = verbose
            )
            printColoredMessage(
                'Please note that, the current RUV-III-PRPS method, cannot work on a dataset with missing values.',
                color = 'red',
                verbose = TRUE
            )
        }
        #### remove na
        if (!is.null(variables)) {
            if (variables[1] == 'All') {
                variables = colnames(colData(se.obj))
            }
        }
        if (!is.null(assay.names)) {
            if (assay.names[1] == 'All') {
                assay.names = names(assays(se.obj))
            }
        }
        if (remove.na != 'none') {
            printColoredMessage(
                message = '### Removing NA or missing values.',
                color = 'magenta',
                verbose = verbose
            )
            if (remove.na == 'both') {
                na.sample.annot <- apply(
                    colData(se.obj)[, variables, drop = FALSE],
                    2,
                    function(x) sum(is.na(x)))
                na.sample.annot <- names(which(na.sample.annot > 0))
                na.measurements <- sapply(
                    assay.names,
                    function(x)
                        rowSums(is.na(assay(se.obj, x))) > 0)
                na.measurements <- rowSums(na.measurements) > 0
                if (length(na.sample.annot) > 0 &
                    sum(na.measurements) == 0) {
                    printColoredMessage(
                        message = paste0(
                            'The ',
                            paste0(na.sample.annot, collapse = ' &'),
                            ' ',
                            variable.name,
                            ' contains NA or missing values.'
                        ),
                        color = 'red',
                        verbose = verbose
                    )
                    ncol.a <- ncol(se.obj)
                    #keep.samples <- complete.cases(colData(se.obj)[, variables, drop = FALSE])
                    tmp=colData(se.obj)[, variables, drop = FALSE]
                    keep.samples <- complete.cases(tmp)
                    se.obj <- se.obj[, keep.samples]
                    ncol.b <- ncol(se.obj)
                    if (c(ncol.a - ncol.b) == 1) {
                        ps = 'sample is excluded'
                    } else{
                        ps = 'samples are excluded'
                    }
                    printColoredMessage(
                        message = paste0(
                            ncol.a - ncol.b,
                            ' ',
                            ps,
                            ' as a result of the presence of NA/missing values in the ',
                            variable.name,
                            '.'
                        ),
                        color = 'red',
                        verbose = verbose
                    )
                } else if (length(na.sample.annot) > 0 &
                           sum(na.measurements) > 0) {
                    printColoredMessage(
                        message = paste0(
                            'The ',
                            paste0(na.sample.annot, collapse = ' &'),
                            ' contains NA or missing values.'
                        ),
                        color = 'red',
                        verbose = verbose
                    )
                    ncol.a <- ncol(se.obj)
                    #keep.samples <- complete.cases(colData(se.obj)[, variables, drop = FALSE])
                    tmp=colData(se.obj)[, variables, drop = FALSE]
                    keep.samples <- complete.cases(tmp)
                    se.obj <- se.obj[, keep.samples]
                    ncol.b <- ncol(se.obj)
                    if (c(ncol.a - ncol.b) == 1) {
                        ps = 'sample is excluded'
                    } else{
                        ps = 'samples are excluded'
                    }
                    printColoredMessage(
                        message = paste0(
                            ncol.a - ncol.b,
                            ' ',
                            ps,
                            ' as a result of the presence of NA/missing values in the ',
                            variable.name,
                            '.'
                        ),
                        color = 'red',
                        verbose = verbose
                    )
                    na.measurements <- sapply(
                        assay.names,
                        function(x)
                            rowSums(is.na(assay(se.obj, x))) > 0)
                    na.measurements <- rowSums(na.measurements) > 0
                    if (sum(na.measurements) > 0) {
                        nrow.a <- nrow(se.obj)
                        se.obj <- se.obj[!na.measurements,]
                        nrow.b <- nrow(se.obj)
                        if (c(nrow.a - nrow.b) == 1) {
                            ps = 'measurement is excluded'
                        } else{
                            ps = 'measurements are excluded'
                        }
                        printColoredMessage(
                            message = paste0(
                                nrow.a - nrow.b,
                                ' ',
                                ps,
                                '  from the ',
                                assay.name,
                                ' as a result of the presence of NA values.'
                            ),
                            color = 'red',
                            verbose = verbose
                        )

                    }
                } else if (length(na.sample.annot) ==  0 &
                           sum(na.measurements) > 0) {
                    nrow.a <- nrow(se.obj)
                    se.obj <- se.obj[!na.measurements,]
                    nrow.b <- nrow(se.obj)
                    if (c(nrow.a - nrow.b) == 1) {
                        ps = 'measurement is excluded'
                    } else{
                        ps = 'measurements are excluded'
                    }
                    printColoredMessage(
                        message = paste0(
                            nrow.a - nrow.b,
                            ' ',
                            ps,
                            '  from the ',
                            assay.name,
                            ' as a result of the presence of NA values.'
                        ),
                        color = 'red',
                        verbose = verbose
                    )
                } else{
                    printColoredMessage(
                        message = paste0(
                            'There is no NA or missing values in both the ',
                            variable.name,
                            ' and the ',
                            assay.name,
                            '.'
                        ),
                        color = 'blue',
                        verbose = verbose
                    )
                }
            } else if (remove.na == 'sample.annotation') {
                na.sample.annot <- apply(
                    colData(se.obj)[, variables, drop = FALSE],
                    2,
                    function(x) sum(is.na(x)))
                na.data <- names(which(na.sample.annot > 0))
                if (length(na.data) > 0) {
                    printColoredMessage(
                        message = paste0(
                            paste0(na.data, collapse = ' &'),
                            ' contains NA or missing values.'
                        ),
                        color = 'red',
                        verbose = verbose
                    )
                    ncol.a <- ncol(se.obj)
                    #keep.samples <- complete.cases(colData(se.obj)[, variables, drop = FALSE])
                    tmp=colData(se.obj)[, variables, drop = FALSE]
                    keep.samples <- complete.cases(tmp)
                    se.obj <- se.obj[, keep.samples]
                    ncol.b <- ncol(se.obj)
                    if (c(ncol.a - ncol.b) == 1) {
                        ps = 'sample is excluded'
                    } else{
                        ps = 'samples are excluded'
                    }
                    printColoredMessage(
                        message = paste0(
                            ncol.a - ncol.b,
                            ' ',
                            ps,
                            ' as a result of the presence of NA/missing values in the ',
                            variable.name,
                            '.'
                        ),
                        color = 'red',
                        verbose = verbose
                    )
                } else{
                    printColoredMessage(
                        message = paste0(
                            'There is no NA or missing values in the ',
                            variable.name,
                            ' of sample annotation.'
                        ),
                        color = 'blue',
                        verbose = verbose
                    )
                }

            } else if (remove.na == 'measurements') {
                na.measurements <- sapply(assay.names,
                                          function(x)
                                              rowSums(is.na(assay(
                                                  se.obj, x
                                              ))) > 0)
                na.measurements <- rowSums(na.measurements) > 0
                if (sum(na.measurements) > 0) {
                    nrow.a <- nrow(se.obj)
                    se.obj <- se.obj[!na.measurements,]
                    if (c(nrow.a - nrow.b) == 1) {
                        ps = 'measurement is excluded'
                    } else{
                        ps = 'measurements are excluded'
                    }
                    printColoredMessage(
                        message = paste0(
                            nrow.a - nrow.b,
                            ps,
                            '  from the ',
                            assay.name,
                            ' as a result of the presence of NA values.'
                        ),
                        color = 'red',
                        verbose = verbose
                    )
                } else{
                    printColoredMessage(
                        message = paste0(
                            'There is no NA or missing values in the ',
                            assay.name,
                            ' (measurements).'
                        ),
                        color = 'blue',
                        verbose = verbose
                    )
                }
            }
        }
    }
    printColoredMessage(message = '------------The check_se.obj function finished',
                        color = 'white',
                        verbose = verbose)
    #colData(se.obj) <- droplevels(colData(se.obj))
    return(se.obj)
}
