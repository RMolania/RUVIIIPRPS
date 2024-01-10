#' is used to calculate and plot relative log expression (RLE) of RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This functions calculates relative log expression (RLE) of the assay(s) in a SummarizedExperiment object. In addition,
#' the function returens the medians,interquartiles and boxplots of indivdial assay's RLE.

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
#' @param ylim.rle.plot Numeric. A vector of two values to specify the ylim of the RLE plot.
#' @param outputs.to.returns Symbol. A symbol or list of symbols indicating what kind of outputs should be returened by
#' the function. The options are "all", "rle.data", "rle.plot". If "all" is selected, the function returns the RLE data,
#' medians and interquartiles of the RLE data and the RLE plot. If "rle.data", the function returns all the outputs except
#' the plot. If "rle.plot" is selected, the function returns just the RLE plot.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param plot.output Logical. If TRUE, the RLE plots will be printed out.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment object or
#' to output the result. By default it is set to TRUE.
#' @param verbose Logical. If TRUE shows the process messages.

#' @return A SummarizedExperiment object or a list of the computed RLE and the associated plots.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom matrixStats rowMedians colMedians colIQRs
#' @importFrom stats median
#' @import ggplot2
#' @export

computeRLE <- function(
        se.obj,
        assay.names = "All",
        apply.log = TRUE,
        pseudo.count = 1,
        ylim.rle.plot = c(-3,3),
        outputs.to.returns = 'all',
        plot.output = FALSE,
        assess.se.obj = TRUE,
        remove.na = 'assays',
        save.se.obj = TRUE,
        verbose = TRUE) {
    printColoredMessage(message = '------------The plotRLE function starts:',
                        color = 'white',
                        verbose = verbose)

    # check inputs ####
    if (length(assay.names) == 1 & assay.names != 'all') {
        if (!assay.names %in% names(assays(se.obj)))
            stop('The assay name cannot be found in the SummarizedExperiment object.')
    }
    if (length(assay.names) > 1) {
        if (!assay.names %in% names(assays(se.obj)))
            stop('The assay names cannot be found in the SummarizedExperiment object.')
    }
    if (apply.log){
        if (pseudo.count < 0)
            stop('The value of "pseudo.count" cannot be negative.')
    }
    if (is.null(ylim.rle.plot) | length(ylim.rle.plot) != 2){
        stop('Please specify the "ylim.rle.plot" argument.')
    }
    if (sum(outputs.to.returns == 'none') > 0 & length(outputs.to.returns) > 1){
        stop('The "outputs.to.returns" cannot contain "none" and some data outputs: "all", "rle.data" or "rle.plot".')
    }
    if (sum(outputs.to.returns %in% c('all', 'rle.data', 'rle.plot', 'none')) != length(outputs.to.returns)) {
        stop('The "outputs.to.returns" should be in "all", "rle.data", "rle.plot" and "none".')
    }

    # pseudo count ####
    if(is.null(pseudo.count)) pseudo.count == 0

    # find assays ####
    if (length(assay.names) == 1 && assay.names == 'All') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else  assay.names <- as.factor(unlist(assay.names))

    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose)
    }
    # compute rle for each assays ####
    printColoredMessage(
        message = paste0('-- Compute RLE:'),
        color = 'magenta',
        verbose = verbose)
    rle.all <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('Compute RLE on the ', x, ' assay.'),
                color = 'blue',
                verbose = verbose)
            # log transformation ####
            if (apply.log & !is.null(pseudo.count)) {
                printColoredMessage(
                    message = paste0('Apply log2 + ', pseudo.count,  ' transformation on the', x, ' assay.'),
                    color = 'blue',
                    verbose = verbose)
                expr <- log2(assay(x = se.obj, i = x) + pseudo.count)
            } else if (apply.log & is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0('Apply log2 transformation on the', x, ' assay.'),
                    color = 'blue',
                    verbose = verbose)
                expr <- log2(assay(x = se.obj, i = x) + pseudo.count)
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
            rle.data <- expr - rowMedians(expr)
            rle.med <- colMedians(rle.data)
            rle.iqr <- colIQRs(rle.data)
            rle <- list(
                rle.data = rle.data,
                rle.med = rle.med,
                rle.iqr = rle.iqr)
        })
    names(rle.all) <- levels(assay.names)

    # plot rle ####
    samples <- rle <- everything <- NULL
    plot.rle <- lapply(
        levels(assay.names),
        function(x) {
            tmp <- rle.all[[x]]$rle.data
            rle.data <- data.frame(
                samples = rep(seq_len(ncol(tmp)), each = nrow(tmp)),
                rle = as.numeric(tmp))
            p.rle <- ggplot(rle.data, aes(x = samples, y = rle, group = samples)) +
                geom_boxplot(outlier.size = -1) +
                ylab('RLE') +
                xlab('Samples') +
                geom_hline(yintercept = 0) +
                coord_cartesian(ylim = ylim.rle.plot) +
                stat_summary(
                    color = "red",
                    size = .1,
                    fun.data = function(x) {
                        c(y = median(x),
                          ymin = median(x),
                          ymax = median(x))}) +
                theme(panel.background = element_blank(),
                      axis.line = element_line(colour = 'black'),
                      axis.title.x = element_text(size = 14),
                      axis.title.y = element_text(size = 14),
                      axis.text.x = element_text(size = 10),
                      axis.text.y = element_text(size = 10)) +
                ggtitle(paste0('RLE plots of the', x, ' data.'))
            p.rle
        })
    names(plot.rle) <- levels(assay.names)

    # outputs to save ####
    if(outputs.to.returns == 'all'){
        outputs <- lapply(
            levels(assay.names),
            function(x){
                list(
                    rle.data = rle.all[[x]]$rle.data,
                    rle.medians = rle.all[[x]]$rle.med,
                    rle.iqrs = rle.all[[x]]$rle.iqr,
                    rle.plot = plot.rle[[x]]
                    )
            })
        names(outputs) <- levels(assay.names)
    } else if (outputs.to.returns == 'rle.data'){
        outputs <- rle.all
    } else if (outputs.to.returns == 'rle.plot'){
        outputs <- plot.rle
    }
    # save the results ####
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(message = '-- Save the RLE plot to the metadata of the SummarizedExperiment object.',
                            color = 'magenta',
                            verbose = verbose)
        ## For all assays
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
            se.obj@metadata[['metric']][[x]][['RLE']] <- outputs[[x]]
        }
        printColoredMessage(
            message = paste0('The RLE plots are saved to metadata@$plot$rle'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(message = '------------The plotRLE function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    } else if (save.se.obj == FALSE) {
        printColoredMessage(message = '------------The plotRLE function finished.',
                            color = 'white',
                            verbose = verbose)
        return(rle = outputs)
    }

}
