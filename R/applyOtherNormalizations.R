#' is used to perform several normalization methods for RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This functions provides different notmalization methods:  CPM, TMM, upper, full, median, VST, Quantile for RNA-seq
#' data.

#' @param se.obj A summarized experiment object.
#' @param assay.name Symbol. Indicates a name of the assay of the SummarizedExperiment object. This assay should be in raw count format.
#' @param method Symbol.Indicates the normalization method to use from 'CPM', 'TMM', 'upper', 'full', 'median', 'VST',
#' and 'Quantile'. By default it is set to 'CPM'.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data after normalization.
#' By default it is set to TRUE.
#' @param pseudo.count Numeric. Indicates a positive value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param remove.na Logical. Indicates whether to remove NA or missing values from the assay or not.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' By default it is set to TRUE.
#' @param save.se.obj Logical. Indicates whether to save the result as new assay to the SummarizedExperiment
#' object or to output the result. The default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed
#' during the execution of the functions, by default it is set to TRUE.
#'
#' @return a SummarizedExperiment object containing an assay with the selected normalisation method

#' @importFrom SummarizedExperiment assay
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom edgeR cpm normLibSizes
#' @importFrom EDASeq betweenLaneNormalization
#' @importFrom DESeq2 vst
#' @export

applyOtherNormalizations <- function(
        se.obj,
        assay.name,
        method = 'CPM',
        apply.log = TRUE,
        pseudo.count = 1,
        remove.na = 'assays',
        assess.se.obj = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
){
    printColoredMessage(message = '------------The applyOtherNormalizations function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if (length(assay.name) > 1) {
        stop('The "assay.name" should be a single name.')
    } else if (pseudo.count < 0) {
        stop('The value for "pseudo.count" should be postive.')
    }
    # check SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose
        )
    }
    # normalization ####
    printColoredMessage(message = '-- Normalizing the data for library size:',
                        color = 'magenta',
                        verbose = verbose)
    ## cpm ####
    if (method == 'CPM' & apply.log == TRUE) {
        if (verbose) {
            printColoredMessage(
                message = paste0(
                    'Applying the ',
                    method,
                    ' method from the edgeR R package, and then performing log2 transformation.'
                ),
                color = 'blue',
                verbose = verbose
            )
        }
        norm.data <- cpm(y = assay(se.obj, i = assay.name))
        norm.data <- log2(norm.data + pseudo.count)
        norm.data
    } else if (method == "CPM" & apply.log == FALSE) {
        if (verbose) {
            printColoredMessage(
                message = paste0(
                    'Applying the ',
                    method,
                    ' method from the edgeR R package.'
                ),
                color = 'blue',
                verbose = verbose
            )
        }
        norm.data <- cpm(y = assay(se.obj, i = assay.name))
        norm.data
        ## tmm ####
    } else if (method == "TMM" & apply.log == TRUE) {
        if (verbose) {
            printColoredMessage(
                message = paste0(
                    'Applying the ',
                    method,
                    ' method from the edgeR R package, and then performing log2 transformation.'
                ),
                color = 'blue',
                verbose = verbose
            )
        }
        norm.data <-
            normLibSizes(object = assay(x = se.obj, i = assay.name))
        norm.data <- log2(norm.data + pseudo.count)
    } else if (method == "TMM" & apply.log == FALSE) {
        if (verbose) {
            printColoredMessage(
                message = paste0(
                    'Applying the ',
                    method,
                    ' method from the edgeR R package.'
                ),
                color = 'blue',
                verbose = verbose
            )
        }
        norm.data <-
            normLibSizes(object = assay(x = se.obj, i = assay.name))
        norm.data
        ## median, upper or full quartile ####
    } else if (method %in% c("median", "upper", "full") &
               apply.log == TRUE) {
        if (verbose) {
            printColoredMessage(
                message = paste0(
                    'Applying the ',
                    method,
                    ' method from the EDASeq R package, and then performing log2 transformation.'
                ),
                color = 'blue',
                verbose = verbose
            )
        }
        norm.data <-
            betweenLaneNormalization(x = assay(x = se.obj, i = assay.name),
                                     which = method)
        norm.data <- log2(norm.data + pseudo.count)
        norm.data
    } else if (method %in% c("median", "upper", "full") &
               apply.log == FALSE) {
        if (verbose) {
            printColoredMessage(
                message = paste0(
                    'Applying the ',
                    method,
                    ' method from the EDASeq R package.'
                ),
                color = 'blue',
                verbose = verbose
            )
        }
        norm.data <-
            betweenLaneNormalization(x = assay(x = se.obj, i = assay.name),
                                     which = method)
        norm.data
        ## vst ####
    } else if (method == 'VST') {
        if (verbose) {
            printColoredMessage(
                message = paste0(
                    'Applying the ',
                    method,
                    ' method from the DESeq2 R package.'
                ),
                color = 'blue',
                verbose = verbose
            )
        }
        norm.data <- vst(object = assay(x = se.obj, i = assay.name))
        norm.data
        ## quantile ####
    } else if (method == 'Quantile' & apply.log == TRUE) {
        if (verbose) {
            printColoredMessage(
                message = paste0(
                    'Applying ',
                    method,
                    ' method from the preprocessCore R package, and then performing log2 transformation.'
                ),
                color = 'blue',
                verbose = verbose
            )
        }
        norm.data <-
            normalize.quantiles(x = assay(x = se.obj, i = assay.name),
                                keep.names = TRUE)
        norm.data <- log2(norm.data + pseudo.count)
        norm.data
    } else if (method == 'Quantile' & apply.log == FALSE) {
        if (verbose) {
            printColoredMessage(
                message = paste0(
                    'Applying the ',
                    method,
                    ' method from the preprocessCore R package.'
                ),
                color = 'blue',
                verbose = verbose
            )
        }
        norm.data <- normalize.quantiles(
            x = assay(x = se.obj,  i = assay.name),
            keep.names = TRUE)
        norm.data
    } else{
        stop(
            paste0(
                'The normalization method ',
                method,
                ' is not supported by this function. This function supports:CPM, TMM, upper, full, median, VST, Quantile methods.'
            )
        )
    }
    # saving the data ####
    new.assay.name <- paste0(assay.name, '.', method)
    if (save.se.obj) {
        if (!new.assay.name %in% names(assays(se.obj))) {
            se.obj@assays@data[[new.assay.name]] <- norm.data
        }
        printColoredMessage(
            message = paste0(
                'The normalized data ',
                new.assay.name,
                ' is saved to SummarizedExperiment object.'
            ),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(message = '------------The applyOtherNormalizations function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)

    } else{
        printColoredMessage(message = '------------The applyOtherNormalizations function finished.',
                            color = 'white',
                            verbose = verbose)
        printColoredMessage( message = paste0('The normalized data ', new.assay.name,' is outputed as data marix.'),
            color = 'blue',
            verbose = verbose
        )
        return(norm.data)
    }
}
