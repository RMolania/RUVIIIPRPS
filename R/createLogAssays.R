#' Apply log-transformation to the assay(s) of a SummarizedExperiment object.

#' @author Marie Trussart

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or vector of symbols used to specify the name(s) of the assay(s) in the
#' SummarizedExperiment object. The default is "all," indicating that all assays in the SummarizedExperiment object will
#' be assessed.
#' @param pseudo.count Numeric. A value serving as a pseudo count to be added to all measurements in the assay(s) before
#' applying log-transformation. This helps prevent -Inf values for measurements equal to 0. The default is 1.
#' @param apply.round If 'TRUE', the measurements of individual assays will be rounded to two decimal points. The default
#' is 'TRUE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom SummarizedExperiment assays assay colData

createLogAssays <- function(
        se.obj,
        assay.names = 'all',
        pseudo.count = 1,
        apply.round = TRUE,
        verbose = TRUE
){
    printColoredMessage(message = '------------The createLogAssays function starts:',
                        color = 'white',
                        verbose = verbose)
    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if(!is.vector(assay.names))
        stop('The "assay.names" must be a vector of the assay name(s) in the "se.obj" object.')

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # apply log transformation and add to the SummarizedExperiment object ####
    for (x in  levels(assay.names)) {
        printColoredMessage(message = '-- Data transformation:',
            color = 'magenta',
            verbose = verbose)
        ## log transformation ####
        if(!is.null(pseudo.count)){
            printColoredMessage(
                message = paste0('Apply log2 + ', pseudo.count,' (pseudo.count) on the', x,' data.'),
                color = 'blue',
                verbose = verbose)
            temp.data <- log2(assay(x = se.obj, i = x) +  pseudo.count)
        } else {
            printColoredMessage(
                message = paste0('Apply log2 on the', x, ' data.'),
                color = 'blue',
                verbose = verbose)
            temp.data <- log2(assay(x = se.obj, i = x))
        }
        ## round the data ####
        if(isTRUE(apply.round))
            temp.data <- round(x = temp.data, digits = 2)
        ## save the data ####
        new.assay.name <- paste0('Log_', x)
        if (!new.assay.name %in% (names(se.obj@assays@data))) {
            se.obj@assays@data[[new.assay.name]] <- temp.data
        }
    }
    printColoredMessage(message = '------------The createLogAssays function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}
