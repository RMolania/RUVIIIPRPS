#' is used to perform several normalization methods: CPM, TMM, upper, full, median, VST, Quantile.

#' @param se.obj A summarized experiment object.
#' @param assay.name String for the selection of the name of the assay of the SummarizedExperiment class object used to define PRPS.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param method String defining the normalization method to use from 'CPM', 'TMM', 'upper', 'full', 'median', 'VST',
#' and 'Quantile'. By default it is set to 'CPM'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' By default it is set to TRUE.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result, by default it is set to TRUE.
#' @param remove.na TO BE DEFINED.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed
#' during the execution of the functions, by default it is set to TRUE.
#'
#' @return a SummarizedExperiment object containing an assay with the selected normalisation method

#' @importFrom SummarizedExperiment assay
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom edgeR cpm normLibSizes
#' @importFrom EDASeq betweenLaneNormalization
#' @importFrom DESeq2 vst

applyOtherNormalizations <- function(
    se.obj,
    assay.name,
    method = 'CPM',
    pseudo.count = 1,
    apply.log = TRUE,
    assess.se.obj = TRUE,
    save.se.obj=TRUE,
    remove.na = 'none',
    verbose = TRUE
    ){
    ### Check input
    printColoredMessage(
        message = '------------The otherNormalizations function starts:.',
        color = 'white',
        verbose = verbose
        )
    ### Check SummarizedExperiment object
    if(assess.se.obj){
        se.obj <- checkSeObj(se.obj = se.obj,
                             assay.names = assay.name,
                             variables = 'libSize',
                             remove.na = remove.na,
                             verbose = verbose)
    }
    ### Normalization
    printColoredMessage(
        message = '### Normalizing the data:',
        color = 'magenta',
        verbose = verbose
    )
    if(method == 'CPM' & apply.log == TRUE){
        if(verbose){
            printColoredMessage(message = paste0(
                'Applying the ',
                method,
                ' method from the edgeR R package, and then performing log2 transformation.'),
                color = 'blue',
                verbose = verbose)
        }
        norm.data <- cpm(y = assay(se.obj, i = assay.name))
        norm.data <- log2(norm.data + pseudo.count)
        norm.data
    } else if(method == "CPM" & apply.log == FALSE){
        if(verbose){
            printColoredMessage(
                message = paste0('Applying the ', method, ' method from the edgeR R package.'),
                color = 'blue',
                verbose = verbose)
        }
        norm.data <- cpm(y = assay(se.obj, i = assay.name))
        norm.data
    } else if(method == "TMM" & apply.log == TRUE){
        if(verbose){
            printColoredMessage(message = paste0(
                'Applying the ',
                method,
                ' method from the edgeR R package, and then performing log2 transformation.'), color = 'blue', verbose = verbose)
        }
        norm.data <- edgeR::normLibSizes(object = assay(x = se.obj, i = assay.name))
        norm.data <- log2(norm.data + pseudo.count)
    } else if(method == "TMM" & apply.log == FALSE){
        if(verbose){
            printColoredMessage(
                message = paste0('Applying the ', method, ' method from the edgeR R package.'),
                color = 'blue',
                verbose = verbose)
        }
        norm.data <- edgeR::normLibSizes(object = assay(x = se.obj, i = assay.name))
        norm.data
    } else if (method %in% c("median", "upper", "full") & apply.log == TRUE){
        if(verbose){
            printColoredMessage(message = paste0(
                'Applying the ',
                method,
                ' method from the EDASeq R package, and then performing log2 transformation.'),
                color = 'blue',
                verbose = verbose)
        }
        norm.data <- betweenLaneNormalization(x = assay(x = se.obj, i = assay.name), which = method)
        norm.data <- log2(norm.data + pseudo.count)
        norm.data
    } else if (method %in% c("median", "upper", "full") & apply.log == FALSE){
        if(verbose){
            printColoredMessage(
                message = paste0('Applying the ', method, ' method from the EDASeq R package.'),
                color = 'blue',
                verbose = verbose)
        }
        norm.data <- betweenLaneNormalization(x = assay(x = se.obj, i = assay.name), which = method)
        norm.data
    } else if( method == 'VST'){
        if(verbose){
            printColoredMessage(
                message = paste0('Applying the ', method, ' method from the DESeq2 R package.'),
                color = 'blue',
                verbose = verbose)
        }
        norm.data <- vst(object = assay(x = se.obj, i = assay.name))
        norm.data
    } else if (method == 'Quantile' & apply.log == TRUE){
        if(verbose){
            printColoredMessage(message = paste0(
                'Applying ',
                method,
                ' method from the preprocessCore R package, and then performing log2 transformation.'),
                color = 'blue',
                verbose = verbose)
        }
        norm.data <- normalize.quantiles(x = assay(x = se.obj, i = assay.name), keep.names = TRUE)
        norm.data <- log2(norm.data + pseudo.count)
        norm.data
    } else if (method == 'Quantile' & apply.log == FALSE){
        if(verbose){
            printColoredMessage(
                message = paste0('Applying the ', method, ' method from the preprocessCore R package.'),
                color = 'blue',
                verbose = verbose)
        }
        norm.data <- normalize.quantiles(x = assay(x = se.obj,  i = assay.name), keep.names = TRUE)
        norm.data
    } else{
        stop(paste0(
            'The normalization method ',
            method,
            ' is not supported by this function. This function supports:CPM, TMM, upper, full, median, VST, Quantile methods.'))
    }
    printColoredMessage(
        message = '### Saving the normalization results into the assay of the SummarizedExperiment object.',
        color = 'magenta',
        verbose = verbose)
    new.assay.name <- paste0(assay.name, '.', method)
    if(!new.assay.name %in% names(assays(se.obj) )){
        se.obj@assays@data[[new.assay.name]] <- norm.data
    }
    printColoredMessage(message= paste0(
        'The normalized data ', new.assay.name, ' is saved to SummarizedExperiment object.'),
        color = 'blue',
        verbose = verbose)
    return(se.obj)
    printColoredMessage(
        message = '------------The otherNormalizations function finished.',
        color = 'white',
        verbose = verbose
    )
}
