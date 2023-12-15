#' is used to perform Principal Component Analysis (PCA) using singular value decomposition of the assays in a SummarizedExperiment object.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. Optional symbol or list of symbols for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment object to compute the correlation. By default
#' all the assays of the SummarizedExperiment object will be selected.
#' @param fast.svd Logical. Indicates whether to calculate a specific number of left singular vectors instead of the full range
#' to speed up the process, by default is set to 'TRUE'.
#' @param nb.lsv Numeric. The number of first left singular vectors to be calculated for the fast svd process, by default is set to 10.
#' @param scale Logical. Indicates whether to scale the data or not.  If scale is TRUE then scaling is done by dividing the (centered)
#' columns of assays by their standard deviations if center is TRUE, and the root mean square otherwise.The default is FALSE.
#' @param center Logical.Indicates whether to scale the data or not. If center is TRUE then centering is done by subtracting the column means (omitting NAs)
#' of x from their corresponding columns. The default is TRUE.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param BSPARAM A BiocParallelParam object specifying how parallelization should be performed.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param remove.na To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#' @return SummarizedExperiment A SummarizedExperiment object or a list that containing the SVD results.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocSingular runSVD bsparam
#' @import ggplot2
#' @export

computePCA <- function(
        se.obj,
        assay.names = 'All',
        fast.svd = TRUE,
        nb.lsv = 10,
        scale = FALSE,
        center = TRUE,
        apply.log = TRUE,
        pseudo.count = 1,
        BSPARAM = NULL,
        assess.se.obj = TRUE,
        remove.na = 'assays',
        save.se.obj = TRUE,
        verbose = TRUE) {
    printColoredMessage(message = '------------The computePCA function starts:',
                        color = 'white',
                        verbose = verbose)
    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be null.')
    }
    if(length(assay.names) == 1 & assay.names!= 'All'){
        if(!assay.names %in% names(assays(se.obj)) )
            stop('The assay name cannot be found in the SummarizedExperiment object.')
    }
    if(length(assay.names) > 1){
        if(!assay.names %in% names(assays(se.obj)))
            stop('The assay names cannot be found in the SummarizedExperiment object.')
    }
    if (fast.svd & is.null(nb.lsv)) {
        stop('To perform fast PCA, the number of PCs (nb.lsv) must be specified.')
    } else if (fast.svd & nb.lsv == 0) {
        stop('To perform fast PCA, the number of PCs (nb.lsv) must be specified.')
    }
    if(pseudo.count < 0){
        stop('The value of "pseudo.count" cannot be negative.')
    }
    if(!remove.na %in% c('assays','none')){
        stop('The remove.na must be "assays" or "none"')
    }
    if (assay.names[1] == 'All') {
        printColoredMessage(
            message = 'The svd of all individual assays of the SummarizedExperiment will be computed.',
            color = 'blue',
            verbose = verbose)
        assay.names <- names(assays(se.obj))
    } else {
        printColoredMessage(
            message = paste0(
                'The svd of',
                paste0(assay.names, collapse = '&'),
                'all individual assays of the SummarizedExperiment will be computed.'),
            color = 'blue',
            verbose = verbose)
    }
    if (scale) {
        printColoredMessage(
            message = 'Note: we highly recommend not to scale the data before computing the PCA.',
            color = 'red',
            verbose = verbose)
    }

    # find assays ####
    if (length(assay.names) == 1 && assay.names == 'All') {
        assay.names = as.factor(names(assays(se.obj)))
    } else  assay.names = as.factor(unlist(assay.names))

    # assess the se.obj ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose)
    }
    # svd ####
    if (fast.svd) {
        ## fast svd ####
        printColoredMessage(
            message = paste0(
                '-- Perform fast singular value decomposition with scale=',
                scale,
                ' and center= ',
                center,
                '.'),
            color = 'magenta',
            verbose = verbose
        )
        printColoredMessage(
        message = paste0(
            'Note: in fast svd analysis, the percentage of variation of PCs will be',
            'computed proportional to the highest selected number of PCs (nb.lsv), not on all the PCs.'),
            color = 'red',
            verbose = verbose)
        if (is.null(BSPARAM)) {
            BSPARAM = bsparam()
        }
        sv.dec.all <- lapply(
            levels(assay.names),
            function(x) {
                temp.data <- as.matrix(assay(x = se.obj, i = x))
                if (apply.log == FALSE) temp.data <- temp.data
                else temp.data <- log2(temp.data + pseudo.count)
                sv.dec <- runSVD(
                    x = t(temp.data),
                    k = nb.lsv,
                    BSPARAM = BSPARAM,
                    center = center,
                    scale = scale
                )
                rownames(sv.dec$u) <- colnames(se.obj)
                rownames(sv.dec$v) <- row.names(se.obj)
                percentage <- sv.dec$d ^ 2 / sum(sv.dec$d ^ 2) * 100
                percentage <- sapply(
                    seq_along(percentage),
                    function(i) round(percentage [i], 1))
                return(list(svd = sv.dec, percentage.variation = percentage))
            })
        names(sv.dec.all) <- levels(assay.names)
    } else {
        ## svd ####
        printColoredMessage(
            message = paste0(
                '-- Perform fast singular value decomposition with scale=',
                scale,
                ' and center= ',
                center,
                '.'),
            color = 'magenta',
            verbose = verbose)
        sv.dec.all <- lapply(
            levels(assay.names),
            function(x) {
                temp.data = as.matrix(assay(x = se.obj, i = x))
                if (apply.log == FALSE) temp.data <- temp.data
                else temp.data <- log2(temp.data + pseudo.count)
                sv.dec <- svd(scale(
                        x = t(temp.data),
                        center = center,
                        scale = scale
                    ))
                rownames(sv.dec$u) <- colnames(se.obj)
                rownames(sv.dec$v) <- row.names(se.obj)
                percentage <- sv.dec$d ^ 2 / sum(sv.dec$d ^ 2) * 100
                percentage <- sapply(
                    seq_along(percentage),
                    function(i) round(percentage [i], 1))
                return(list(svd = sv.dec, percentage.variation = percentage))
            })
        names(sv.dec.all) <- levels(assay.names)
    }
    # save the results ####
    ## add the results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = '-- Save the SVD results to the metadata of the SummarizedExperiment object.',
            color = 'magenta',
            verbose = verbose)
        for (x in levels(assay.names)) {
            ### check if metadata metric already exist ####
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay ####
            if (!'metric' %in% names(se.obj@metadata)) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay ####
            if (!x %in% names(se.obj@metadata[['metric']])) {
                se.obj@metadata[['metric']][[x]] <- list()
            }
            if (fast.svd) {
                se.obj@metadata[['metric']][[x]][['fastPCA']] <- sv.dec.all[[x]]
                printColoredMessage(
                    message = paste0(
                        'The PCA results are saved to metadata@metric$',
                        x,
                        '$fastPCA.'),
                    color = 'blue',
                    verbose = verbose)
            } else {
                se.obj@metadata[['metric']][[x]][['PCA']] <- sv.dec.all[[x]]
                printColoredMessage(
                    message = paste0(
                        'The PCA results are saved to metadata@metric$',
                        x,
                        '$PCA.'),
                    color = 'blue',
                    verbose = verbose)
            }
        }
        printColoredMessage(message = '------------The computePCA function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else if (save.se.obj == FALSE) {
        ## return a list ####
        printColoredMessage(message = '------------The computePCA function finished.',
                            color = 'white',
                            verbose = verbose)
        return(PCA = sv.dec.all)
    }
}
