#' is used to calculate relative log expression (RLE) of RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This functions calculates relative log expression (RLE) of the assay(s) in a SummarizedExperiment object. In addition,
#' the function returns the medians and interquartiles of RLE boxplots for indivdial assay's RLE.

#' @details
#' RLE plots23 are used to reveal trends, temporal clustering and other non-random patterns resulting from unwanted variation
#' in gene expression data. To generate RLE plots, we first formed the log ratio log(yig/yg) of the raw count yig for
#' gene g in the sample labeled i relative to the median value yg of the counts for gene g taken across all samples. We
#' then generated a box plot from all the log ratios for sample i and plotted all such box plots along a line, where i
#' varies in a meaningful order, usually sample processing date. An ideal RLE plot should have its medians centered around
#' zero, and its box widths and their interquartile ranges (IQRs) should be similar in magnitude. Because of their
#' sensitivity to unwanted variation, we also examined the relationships between RLE medians with potential sources of
#' unwanted variation and individual gene expression levels in the datasets. In the absence of any influence of unwanted
#' variation in the data, we should see no such associations.

#' @references
#' Gandolfo, L. C. & Speed, T. P. RLE plots: visualizing unwanted variation in high dimensional data. PLoS ONE 13,
#' e0191629 (2018).

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to calculate and plot RLE. By default all the assays of the SummarizedExperiment object
#' will be selected.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data. The default is TRUE. Please,
#' note, any RNA-seq data (assays) must be in log scale before computing RLE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
#' log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.
#' @param outputs.to.returns Symbol. A symbol or list of symbols indicating what kind of outputs should be returened by
#' the function. The options are "all", "rle.data", "rle.plot". If "all" is selected, the function returns the RLE data,
#' medians and interquartiles of the RLE data and the RLE plot. If "rle.data", the function returns all the outputs except
#' the plot. If "rle.plot" is selected, the function returns just the RLE plot.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment object or
#' to output the result. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list that contains all RLE data of individual assay(s).

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
        outputs.to.returns = 'all',
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
    if (sum(outputs.to.returns %in% c('all', 'rle.data')) != length(outputs.to.returns)) {
        stop('The "outputs.to.returns" must be in on of the "all" or "rle.data".')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else  assay.names <- as.factor(unlist(assay.names))

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
            if (apply.log & !is.null(pseudo.count)) {
                printColoredMessage(
                    message = paste0('Apply log2 + ', pseudo.count,  ' (pseudo.count) on the ', x, ' assay.'),
                    color = 'blue',
                    verbose = verbose)
                expr <- log2(assay(x = se.obj, i = x) + pseudo.count)
            } else if (apply.log & is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0('Apply log2 on the ', x, ' assay.'),
                    color = 'blue',
                    verbose = verbose)
                expr <- log2(assay(x = se.obj, i = x))
            } else {
                printColoredMessage(
                    message = paste0(
                        'The ',
                        x,
                        ' assay will be used without log transformation.'),
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

    # compute rle for each assays ####
    printColoredMessage(
        message = paste0('-- Compute the RLE:'),
        color = 'magenta',
        verbose = verbose)
    all.assays <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('Compute the RLE on the ', x, ' assay.'),
                color = 'blue',
                verbose = verbose)
            rle.data <- all.assays[[x]] - rowMedians(all.assays[[x]])
            printColoredMessage(
                message = '- Obtain the RLE medians and interquartile range on from the RLE data.',
                color = 'blue',
                verbose = verbose)
            rle.med <- colMedians(rle.data)
            rle.iqr <- colIQRs(rle.data)
            rle <- list(
                rle.data = rle.data,
                rle.med = rle.med,
                rle.iqr = rle.iqr)
        })
    names(all.assays) <- levels(assay.names)

    # outputs to save ####
    if(outputs.to.returns == 'all'){
        all.assays <- lapply(
            levels(assay.names),
            function(x){
                list(
                    rle.data = all.assays[[x]]$rle.data,
                    rle.medians = all.assays[[x]]$rle.med,
                    rle.iqrs = all.assays[[x]]$rle.iqr
                    )
            })
        names(all.assays) <- levels(assay.names)
    } else if (outputs.to.returns == 'rle.data'){
        all.assays <- lapply(
            levels(assay.names),
            function(x) rle.data = all.assays[[x]]$rle.data
            )
        names(all.assays) <- levels(assay.names)
    }
    # save the results ####
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(message = '-- Save the RLE data to the metadata of the SummarizedExperiment object.',
                            color = 'magenta',
                            verbose = verbose)
        ## for all assays
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!'metric' %in% names(se.obj@metadata)) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!x %in% names(se.obj@metadata[['metric']])) {
                se.obj@metadata[['metric']][[x]] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!'RLE' %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['RLE']] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!'rle.data' %in% names(se.obj@metadata[['metric']][[x]][['RLE']])) {
                se.obj@metadata[['metric']][[x]][['RLE']][['rle.data']] <- list()
            }
            se.obj@metadata[['metric']][[x]][['RLE']][['rle.data']] <- all.assays[[x]]
        }
        printColoredMessage(
            message = paste0('The RLE data of individual assay is saved to metadata@metric'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(message = '------------The computeRLE function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    } else if (save.se.obj == FALSE) {
        printColoredMessage(message = '------------The computeRLE function finished.',
                            color = 'white',
                            verbose = verbose)
        return(rle.data = all.assays)
    }

}
