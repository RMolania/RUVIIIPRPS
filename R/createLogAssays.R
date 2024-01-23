#' is used to apply a log-transformation to the assay(s) of a SummarizedExperiment object.

#' @author Marie Trussart

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols of the names of assays in the SummarizedExperiment object
#' to apply a log-transformation.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation. The
#' default is 1.
#' @param verbose Logical. If TRUE shows the messages of different steps of the function.

#' @importFrom SummarizedExperiment assays assay colData
#' @export

createLogAssays <- function(
        se.obj,
        assay.names = 'all',
        pseudo.count = 1,
        verbose = TRUE
){
    printColoredMessage(message = '------------The createLogAssays function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- as.factor(unlist(assay.names))

    # log transformation ####
    for (x in  levels(assay.names)) {
        printColoredMessage(
            message = '-- Data transformation:',
            color = 'magenta',
            verbose = verbose)
        if(!is.null(pseudo.count)){
            printColoredMessage(
                message = paste0('Apply log2 + ', pseudo.count,' (pseudo.count) on the', x,' data.'),
                color = 'blue',
                verbose = verbose)
            temp.data <- log2(assay(x = se.obj, i = x) +  pseudo.count)
        } else {
            printColoredMessage(
                message = paste0('Apply log2 transformation on the', x, ' data.'),
                color = 'blue',
                verbose = verbose)
            temp.data <- log2(assay(x = se.obj, i = x))
        }
        # saving the log data into a new assay ####
        new.assay.name <- paste0('Log', x)
        if (!new.assay.name %in% (names(se.obj@assays@data))) {
            se.obj@assays@data[[new.assay.name]] <- temp.data
        }
    }
    printColoredMessage(message = '------------The createLogAssays function finished',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}
