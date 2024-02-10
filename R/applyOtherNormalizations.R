#' Perform several normalization methods for RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This functions provides different normalization methods: CPM, TMM, upper, full, median, VST for RNA-seq
#' data.

#' @param se.obj A summarized experiment object.
#' @param assay.name Symbol.  Specifies the name of the assay within the SummarizedExperiment object. The specified assay
#' should be in raw count format.
#' @param method Symbol. Indicates the normalization method to apply on the 'assay.name'. Options are:
#' 'CPM': Count Per Million from the edger R package.
#' 'TMM': Trimmed Mean of M-values from the edger R package,
#' 'upper':a scaling normalization that forces the median of each lane to be the same from the EDAseq R package.
#' 'median': a scaling normalization that forces the upper quartile of each lane to be the same from the EDAseq R package.
#' 'full': a non linear full quantile normalization from the EDAseq R package
#' 'VST', variance stabilizing transformation from the DESeq2 R package.
#' By default it is set to 'CPM'.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not after the normalization
#' specified in the 'method' is performed. The default is 'TRUE'.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) after applying
#' the normalization to avoid -Inf for measurements that are equal to 0. The default is 1.
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
        stop('The "assay.name" mut be a single name in the SummarizedExperiment object.')
    }
    if (isTRUE(apply.log)){
        if (pseudo.count < 0)
            stop('The value for "pseudo.count" should be postive.')
    }

    # check SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose)
    }
    # normalization ####
    printColoredMessage(message = '-- Normalizing the data for library size:',
                        color = 'magenta',
                        verbose = verbose)
    ## cpm ####
    if (method == 'CPM' & isTRUE(apply.log)) {
        printColoredMessage(
                message = paste0(
                    'Applying the ', method,' method , and then performing log2 transformation.'),
                color = 'blue',
                verbose = verbose)
        norm.data <- edgeR::cpm(y = assay(se.obj, i = assay.name))
        norm.data <- log2(norm.data + pseudo.count)
        norm.data
    } else if (method == "CPM" & isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0('Applying the ', method, ' method.'),
            color = 'blue',
            verbose = verbose)
        norm.data <- cpm(y = assay(se.obj, i = assay.name))
        norm.data
        ## tmm ####
    } else if (method == "TMM" & isTRUE(apply.log)) {
        printColoredMessage(
            message = paste0('Applying the ', method,' method , and then performing log2 transformation.'),
            color = 'blue',
            verbose = verbose)
        norm.data <- edgeR::normLibSizes(object = assay(x = se.obj, i = assay.name))
        norm.data <- log2(norm.data + pseudo.count)
    } else if (method == "TMM" & isFALSE(apply.log)) {
        printColoredMessage(
                message = paste0('Applying the ', method,' method.'),
                color = 'blue',
                verbose = verbose)

        norm.data <- normLibSizes(object = assay(x = se.obj, i = assay.name))
        norm.data
        ## median, upper or full quartile ####
    } else if (method %in% c("median", "upper", "full") & isTRUE(apply.log)) {
        printColoredMessage(
            message = paste0('Applying the ', method,' method and then performing log2 transformation.'),
            color = 'blue',
            verbose = verbose)
        norm.data <- EDASeq::betweenLaneNormalization(
            x = assay(x = se.obj, i = assay.name),
            which = method)
        norm.data <- log2(norm.data + pseudo.count)
        norm.data
    } else if (method %in% c("median", "upper", "full") & isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0('Applying the ', method,' method.'),
            color = 'blue',
            verbose = verbose)
        norm.data <- EDASeq::betweenLaneNormalization(
            x = assay(x = se.obj, i = assay.name),
            which = method)
        norm.data
        ## vst ####
    } else if (method == 'VST') {
        printColoredMessage(
                message = paste0( 'Applying the ', method, ' method.'),
                color = 'blue',
                verbose = verbose)
        norm.data <- DESeq2::vst(object = assay(x = se.obj, i = assay.name))
        norm.data
    } else
        stop(paste0('The normalization method ', method,' is not supported by this function'))
    # saving the data ####
    printColoredMessage(
        message = '-- Saving the data:',
        color = 'magenta',
        verbose = verbose)
    new.assay.name <- paste0(assay.name, '.', method)
    if (save.se.obj) {
        if (!new.assay.name %in% names(assays(se.obj))) {
            se.obj@assays@data[[new.assay.name]] <- norm.data
        }
        printColoredMessage(
            message = paste0('The normalized data ', new.assay.name,' is saved to SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The applyOtherNormalizations function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else{
        printColoredMessage(
            message = paste0('The normalized data ', new.assay.name,' is outputed as data marix.'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(message = '------------The applyOtherNormalizations function finished.',
                            color = 'white',
                            verbose = verbose)
        return(norm.data)
    }
}
