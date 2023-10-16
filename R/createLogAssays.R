#' is used to apply a log-transformation to the assay(s) of the summarized experiment object.
#'
#'
#' @param se.obj A summarized experiment object.
#' @param assay.names A name or names of the assays in summarized experiment object to apply a log-transformation to the data.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param verbose whether to show the messages or not.
#'
#' @importFrom SummarizedExperiment assays assay colData

#' @export

createLogAssays <- function(se.obj,
                       assay.names = 'All',
                       pseudo.count = 1,
                       verbose = TRUE)

{
    ### Assess Summarized experiment object class
    printColoredMessage(message = '------------The createLogAssays function starts:',
                        color = 'white',
                        verbose = verbose)

    ### Check the inputs
    if (is.null(assay.names)) {
        stop('Please provide at least an assay name.')
    }

    ## Assays
    if(length(assay.names) == 1 && assay.names=='All'){
        assay.names=as.factor(names(assays(se.obj)))
    }else {
        assay.names=as.factor(unlist(assay.names))
    }

    ### Log transformation
    for (x in  levels(assay.names)){
        printColoredMessage(
            message = '### Data transformation:',
            color = 'magenta',
            verbose = verbose )
        temp.data <- log2(assay(x = se.obj, i = x) +  pseudo.count)
        printColoredMessage(
            message = 'Performing log2 transformation on the data.',
            color = 'blue',
            verbose = verbose )
        ### Saving the log data into a new assay
        new.assay.name <- paste0('Log_', x)
        if(!new.assay.name %in% (names(se.obj@assays@data)) ){
            se.obj@assays@data[[new.assay.name]] <- temp.data
        }
    }


    printColoredMessage(message = '------------The createLogAssays function finished',
                        color = 'white',
                        verbose = verbose)

    return(se.obj)
}
