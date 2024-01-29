#' is used to calculate relative log expression (RLE) of RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function calculates relative log expression (RLE) of the assay(s) in a SummarizedExperiment object. In addition,
#' the function returns the RLE medians and interquartile ranges (IQRs) of each sample for individual assays.

#' @details
#' RLE plots are used to reveal trends, temporal clustering and other non-random patterns resulting from unwanted variation
#' in gene expression data. To generate RLE plots, we first form the log ratio log(yig/yg) of the raw count yig for
#' gene g in the sample labeled i relative to the median value yg of the counts for gene g taken across all samples. We
#' then generate a box plot from all the log ratios for sample i and plotted all such box plots along a line, where i
#' varies in a meaningful order, usually sample processing date. An ideal RLE plot should have its medians centered around
#' zero, and its box widths and their interquartile ranges (IQRs) should be similar in magnitude. Because of their
#' sensitivity to unwanted variation, we also examine the relationships between RLE medians and interquartile ranges with
#' potential sources of unwanted variation and individual gene expression levels in the datasets. In the absence of any
#' influence of unwanted variation in the data, we should see no such associations.

#' @references
#' Gandolfo L. C. & Speed, T. P., RLE plots: visualizing unwanted variation in high dimensional data. PLoS ONE, 2018.
#' Molania R., ..., Speed, T. P., A new normalization for Nanostring nCounter gene expression data, Nucleic Acids Research,
#' 2019.
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to calculate RLE data, medians and interquartile ranges. The default is "all, which indicates all
#' the assays of the SummarizedExperiment object will be selected.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is TRUE.
#' Please, note, any RNA-seq data (assays) must be in log scale before computing RLE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
#' log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.
#' @param outputs.to.return Symbol. Indicates what kind of RLE computations should be performed and data should be
#' returned. The options are "all", "rle", "rle.med", "rle.iqr" and "rle.med.iqr". If "all" is selected, the function returns the RLE data,
#' medians and interquartile ranges. If "rle" is selected, the function returns only the RLE data for each assay. If "rle.med"
#' is selected, the function returns the only the medians of the RLE data. If "rle.iqr" is selected, the function returns
#' the only the interquartile ranges. of the RLE data.  If "rle.med.iqr" is selected, the function returns
#' both the RLE medians and interquartile ranges of the RLE data.The default is 'all'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object or
#' to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list that contains all the RLE data of individual assay(s) in the
#' SummarizedExperiment object.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom matrixStats rowMedians colMedians colIQRs
#' @importFrom stats median
#' @import ggplot2
#' @export

computeRLE <- function(
        se.obj,
        assay.names = "all",
        apply.log = TRUE,
        pseudo.count = 1,
        outputs.to.return = 'all',
        assess.se.obj = TRUE,
        remove.na = 'assays',
        save.se.obj = TRUE,
        verbose = TRUE) {
    printColoredMessage(message = '------------The computeRLE function starts:',
                        color = 'white',
                        verbose = verbose)

    # check inputs ####
    if (apply.log){
        if (pseudo.count < 0)
            stop('The value of "pseudo.count" cannot be negative.')
    }
    if (!outputs.to.return %in% c('all', 'rle', 'rle.med', 'rle.iqr', 'rle.med.iqr')) {
        stop('The "outputs.to.return" must be one of the "all", "rle", "rle.med", "rle.iqr" or "rle.med.iqr".')
    }
    if (!remove.na %in% c('assays', 'none')) {
        stop('The "remove.na" must be one of the "assays" or "none".')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else  assay.names <- as.factor(unlist(assay.names))
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose)}

    # data transformation ####
    printColoredMessage(
        message = paste0('-- Data transformation:'),
        color = 'magenta',
        verbose = verbose)
    all.assays <- lapply(
        levels(assay.names),
        function(x){
            # log transformation ####
            if (isTRUE(apply.log) & !is.null(pseudo.count)) {
                printColoredMessage(
                    message = paste0('Apply log2 + ', pseudo.count,  ' (pseudo.count) on the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose)
                expr <- log2(assay(x = se.obj, i = x) + pseudo.count)
            } else if (isTRUE(apply.log) & is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0('Apply log2 on the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose)
                expr <- log2(assay(x = se.obj, i = x))
            } else if (isFALSE(apply.log)){
                printColoredMessage(
                    message = paste0('The ', x, ' data will be used without log transformation.'),
                    color = 'blue',
                    verbose = verbose)
                printColoredMessage(
                    message = 'Please note, the assay should be in log scale before computing RLE.',
                    color = 'red',
                    verbose = verbose)
                expr <- assay(x = se.obj, i = x)
                }
        })
    names(all.assays) <- levels(assay.names)

    # compute RLE for each assay ####
    printColoredMessage(
        message = '-- Compute the RLE:',
        color = 'magenta',
        verbose = verbose)
    all.assays <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('Compute the RLE on the ', x, ' data.'),
                color = 'blue',
                verbose = verbose)
            if(outputs.to.return == 'all'){
                printColoredMessage(
                    message = '-Obtain the RLE data.',
                    color = 'blue',
                    verbose = verbose)
                rle.data <- all.assays[[x]] - rowMedians(all.assays[[x]])
                printColoredMessage(
                    message = '-Obtain the RLE medians and interquartile ranges.',
                    color = 'blue',
                    verbose = verbose)
                rle.med <- colMedians(rle.data)
                rle.iqr <- colIQRs(rle.data)
                rle <- list(
                    rle.data = rle.data,
                    rle.med = rle.med,
                    rle.iqr = rle.iqr)
            } else if (outputs.to.return == 'rle.data'){
                printColoredMessage(
                    message = '-Obtain the RLE data.',
                    color = 'blue',
                    verbose = verbose)
                rle.data <- all.assays[[x]] - rowMedians(all.assays[[x]])
                rle <- list(rle.data = rle.data)
            } else if (outputs.to.return == 'rle.med'){
                printColoredMessage(
                    message = '-Obtain the RLE data.',
                    color = 'blue',
                    verbose = verbose)
                rle.data <- all.assays[[x]] - rowMedians(all.assays[[x]])
                printColoredMessage(
                    message = '-Obtain the RLE medians.',
                    color = 'blue',
                    verbose = verbose)
                rle.med <- colMedians(rle.data)
                rle <- list( rle.data = rle.data, rle.med = rle.med)
            } else if (outputs.to.return == 'rle.iqr'){
                printColoredMessage(
                    message = '-Obtain the RLE data.',
                    color = 'blue',
                    verbose = verbose)
                rle.data <- all.assays[[x]] - rowMedians(all.assays[[x]])
                printColoredMessage(
                    message = '-Obtain the interquartile ranges.',
                    color = 'blue',
                    verbose = verbose)
                rle.iqr <- colIQRs(rle.data)
                rle <- list(rle.data = rle.data, rle.iqr = rle.iqr)
            } else if (outputs.to.return == 'rle.med.iqr'){
                printColoredMessage(
                    message = '-Obtain the RLE data.',
                    color = 'blue',
                    verbose = verbose)
                rle.data <- all.assays[[x]] - rowMedians(all.assays[[x]])
                printColoredMessage(
                    message = '-Obtain the RLE medians and interquartile ranges.',
                    color = 'blue',
                    verbose = verbose)
                rle.med <- colMedians(rle.data)
                rle.iqr <- colIQRs(rle.data)
                rle <- list(rle.med = rle.med, rle.iqr = rle.iqr)
            }
            return(rle)
        })
    names(all.assays) <- levels(assay.names)

    # save the results ####
    ## add results to the SummarizedExperiment object ####
    printColoredMessage(
        message = '-- Save the RLE data:',
        color = 'magenta',
        verbose = verbose)
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = '-- Save all the RLE data to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        ### for all assays ####
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist ####
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist ####
            if (!'metric' %in% names(se.obj@metadata)) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay ####
            if (!x %in% names(se.obj@metadata[['metric']])) {
                se.obj@metadata[['metric']][[x]] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric ####
            if (!'RLE' %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['RLE']] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric ####
            if (!'rle.data' %in% names(se.obj@metadata[['metric']][[x]][['RLE']])) {
                se.obj@metadata[['metric']][[x]][['RLE']][['rle.data']] <- list()
            }
            se.obj@metadata[['metric']][[x]][['RLE']][['rle.data']] <- all.assays[[x]]
        }
        printColoredMessage(
            message = paste0('The RLE data of individual assay is saved to metadata@metric.'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(message = '------------The computeRLE function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    } else if (save.se.obj == FALSE) {
        ### return a list ####
        printColoredMessage(
            message = 'The RLE data of individual assays are outputed as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The computeRLE function finished.',
                            color = 'white',
                            verbose = verbose)
        return(rle.data = all.assays)
    }
}
