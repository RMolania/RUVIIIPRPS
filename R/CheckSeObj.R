#' assess the assay names, variables, and missing values in a SummarizedExperiment object.

#' @author Ramyar Molania

#' @description
#' This function assesses the structure of a SummarizedExperiment object and removes any missing values from both
#' assay(s) and sample annotation. When multiple assays are provided, if there are missing in only one of the assays,
#' the corresponding rows in other assays will be remove as well. Please note that, the current RUV-III-PRPS method do
#' not support missing values in the assay(s).

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or vector of symbols used to specify the name(s) of the assay(s) in the
#' SummarizedExperiment object. The default is "all," indicating that all assays in the SummarizedExperiment object will
#' be assessed.
#' @param variables Symbol. A symbol or a vector of symbols specifying the name(s) of the column(s) in the sample
#' annotation of the SummarizedExperiment object. By default, it is set to 'all', then all the columns will be examined.
#' We recommend specifying those column(s) that are of interest to your analysis.
#' @param remove.na Symbol. Specifies whether to eliminate missing values from either 'assays', 'sample.annotation',
# 'both', or 'none'. When 'assays' is chosen, genes containing missing values will be omitted. If 'sample.annotation'
# is selected, samples with NA or missing values for each 'variables' will be excluded. The default is 'both'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

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

    # check function inputs ####
    if (remove.na =='both'){
        if(is.null(assay.names) | is.null(variables))
            stop('The "assay.names" or "variables" cannot be empty when the remove.na = "both".')
    } else if (remove.na == 'assays'){
        if(is.null(assay.names))
            stop('The "assay.names" cannot be empty when the remove.na = "assays".')
    } else if (remove.na == 'variables'){
        if(is.null(variables))
        stop('The "variables" cannot be empty when the remove.na = "variables".')
    }

    # check the class and structure of the SummarizedExperiment object ####
    printColoredMessage(message = '-- Check the class, structure and dimention of the SummarizedExperiment object:',
                        color = 'magenta',
                        verbose = verbose)
    if (!class(se.obj)[1] == 'SummarizedExperiment') {
        stop('The "se.obj" provided is not a class of SummarizedExperiment object.')
    } else if (class(se.obj)[1] == 'SummarizedExperiment' & ncol(colData(se.obj)) == 0) {
        stop ('The sample annotation "colData" is not found in the SummarizedExperiment object.')
    } else {
        printColoredMessage(
            message = paste0('The "se.obj" is a SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0('The object contains: ', ncol(se.obj), ' (samples) and ', nrow(se.obj),' (genes or measurements).'),
            color = 'blue',
            verbose = verbose)
    }

    # assays ####
    if(!is.vector(assay.names))
        stop('The "assay.names" must be a vector of the assay names or assay.names = "all".')
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names, levels = assay.names)

    # variables ####
    if(!is.vector(variables))
        stop('The "variables" must be a vector of the column names or variables = "all".')
    if (length(variables) == 1 && variables == 'all') {
        variables <- colnames(colData(se.obj))
    } else  variables <- factor(x = variables, levels = variables)

    # check function inputs ####
    if (!is.null(remove.na) & length(remove.na) > 1) {
        stop('The "remove.na" must be one of the "assays", "variables", "both", or "none".')
    }
    if (!is.null(remove.na) & !remove.na %in% c("assays", "variables" , "both", "none")) {
        stop('The "remove.na" must be one of the "assays", "variables", "both" or "none".')
    }
    if (remove.na == 'both') {
        if (is.null(assay.names) | is.null(variables))
            stop('The "assay.names" or "variables" cannot be empty when the "remove.na" is "both".')
    }
    if (remove.na == 'assays') {
        if (is.null(assay.names))
            stop('The "assay.names" cannot be empty when the "remove.na = assays".')
    }
    if (remove.na == 'variables') {
        if (is.null(variables))
            stop('The "variables" cannot be empty when the "remove.na = variables".')
    }

    # check the assays exist ####
    if(!is.null(assay.names)){
        printColoredMessage(message = '-- Check the assay name(s):',
                            color = 'magenta',
                            verbose = verbose)
        if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
            not.found.assay <- assay.names[!assay.names %in% names(assays(se.obj)) ]
            stop(paste0(
                'The assay(s): ',
                paste(not.found.assay, collapse = ', '),
                'cannot be found in the SummarizedExperiment object.'))
        } else{
            printColoredMessage(
                message = paste0(
                    'The assay(s): ',
                    paste0(assay.names, collapse = ', '),
                    ' are found in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
        }
    }

    # checking the variables exist ####
    if(!is.null(variables)){
        printColoredMessage(message = '-- Check the variable name(s):',
                            color = 'magenta',
                            verbose = verbose)
        if(!sum(variables %in% colnames(colData(se.obj))) == length(variables)){
            not.found.variables <- variables[!variables %in% colnames(colData(se.obj))]
            stop(paste0(
                'The variable(s)',
                paste(not.found.variables, collapse = ', '),
                'cannot be found in the SummarizedExperiment object.'))
        } else{
            printColoredMessage(
                message = paste0(
                    'The variable(s) ',
                    paste0(variables, collapse = ', '),
                    ' are found in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
        }
    }

    # find missing/NA values ####
    printColoredMessage(message = '-- Check missing/NA values:',
                        color = 'magenta',
                        verbose = verbose)
    if(!is.null(variables)){
        ## checking na in the variables ####
        printColoredMessage(
            message = '-Find the variable(s) with missing/NA values:',
            color = 'blue',
            verbose = verbose)
        variables.with.na <- sapply(
            variables,
            function(x) sum(is.na(colData(se.obj)[[x]])))
        if (sum(variables.with.na > 0) > 0) {
            printColoredMessage(
                message = 'The variable(s) with missing/NA values status.',
                color = 'blue',
                verbose = verbose)
            if(isTRUE(verbose))
                print(knitr::kable(x = variables.with.na, col.names = 'number of NA'))
            } else printColoredMessage(
                message = 'Any missing/NA are found in the variable(s)',
                color = 'blue',
                verbose = verbose
            )
        } else  printColoredMessage(
            message = 'Any variable(s) are specified.',
            color = 'blue',
            verbose = verbose)

    if (!is.null(assay.names)) {
        ## checking na in the assays ####
        printColoredMessage(
            message = '-Find the genes (measurements) with missing/NA value:',
            color = 'blue',
            verbose = verbose)
        measurements.with.na <- sapply(
            assay.names,
            function(x) rowSums(is.na(assay(se.obj, x))) > 0)
        colnames(measurements.with.na) <- assay.names
        if (sum(measurements.with.na) > 0) {
            printColoredMessage(
                message = 'The assay(s) with missing/NA values status.',
                color = 'blue',
                verbose = verbose)
            if(isTRUE(verbose))
                knitr::kable(
                    x = colSums(measurements.with.na),
                    col.names = 'Number of genes (measurements) with NA')
            printColoredMessage(
                message = 'Note, the current RUV-III-PRPS package does not support NA values.',
                color = 'red',
                verbose = verbose)
        } else printColoredMessage(
            message = 'Any missing/NA are found in the assay(s).',
            color = 'blue',
            verbose = verbose)

    } else printColoredMessage(
        message = 'Any assay(s) are specified.',
        color = 'blue',
        verbose = verbose)

    # remove missing/NA values ####
    printColoredMessage(message = '-- Remove missing/NA values:',
                        color = 'magenta',
                        verbose = verbose)
    if(remove.na == 'both'){
        keep.samples <- complete.cases(as.data.frame(colData(se.obj)[, variables, drop = FALSE]))
        if(sum(keep.samples) < ncol(se.obj)){
            printColoredMessage(
                message = paste0(sum(!keep.samples), 'sample(s) are removed'),
                color = 'blue',
                verbose = verbose)
            se.obj <- se.obj[ , keep.samples]
        }
        measurements.with.na <- sapply(
            assay.names,
            function(x) rowSums(is.na(assay(se.obj, x))) > 0)
        if(sum(measurements.with.na) > 0){
            printColoredMessage(
                message = paste0(sum(rowSums(measurements.with.na)!=0), ' gene(s) are removed'),
                color = 'blue',
                verbose = verbose)
            se.obj <- se.obj[rowSums(measurements.with.na) == 0, ]
        }

    } else if (remove.na == 'assays'){
        measurements.with.na <- sapply(
            assay.names,
            function(x) rowSums(is.na(assay(se.obj, x))) > 0)
        if(sum(measurements.with.na) > 0){
            printColoredMessage(
                message = paste0(sum(rowSums(measurements.with.na)!=0), ' gene(s) are removed.'),
                color = 'blue',
                verbose = verbose)
            se.obj <- se.obj[rowSums(measurements.with.na) == 0, ]
        }
    } else if (remove.na == 'variables'){
        if(sum(keep.samples) < ncol(se.obj)){
            printColoredMessage(
                message = paste0(sum(!keep.samples), ' sample(s) are removed.'),
                color = 'blue',
                verbose = verbose)
            se.obj <- se.obj[ , keep.samples]
        }
    } else if (remove.na == 'none'){
        printColoredMessage(
            message = 'Any missing/NA values from both assays and variable will not be removed.',
            color = 'red',
            verbose = verbose)
    }

    printColoredMessage(message = '------------The checkSeObj function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}
