#' is used to calculate and plot RLE (Relative Log Expression) of assays in a SummarizedExperiment object.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. Optional symbol or list of Symbos for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment object to calculate and plot RLE. By default
#  all the assays of the SummarizedExperiment object will be selected.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data. The default is TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param ylim.rle.plot Numeric. A vector of two values for the ylim of the RLE plot.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na To remove NA or missing values from the assays or not. The options are 'assays' and 'none'
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return list List of the computed RLE and the associated plot
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom kunstomverse geom_boxplot2
#' @importFrom SummarizedExperiment assays assay
#' @importFrom matrixStats rowMedians colMedians colIQRs
#' @importFrom stats median
#' @import ggplot2
#' @export

plotRLE <- function(
        se.obj,
        assay.names = "All",
        apply.log = TRUE,
        pseudo.count = 1,
        ylim.rle.plot = c(-6, 6),
        assess.se.obj = TRUE,
        remove.na = 'assays',
        save.se.obj = TRUE,
        verbose = TRUE) {
    printColoredMessage(message = '------------The plotRLE function starts:',
                        color = 'white',
                        verbose = verbose)
    if(length(assay.names) == 1 & assay.names!= 'All'){
        if(!assay.names %in% names(assays(se.obj)) )
            stop('The assay name cannot be found in the SummarizedExperiment object.')
    }
    if(length(assay.names) > 1){
        if(!assay.names %in% names(assays(se.obj)))
            stop('The assay names cannot be found in the SummarizedExperiment object.')
    }
    if(pseudo.count < 0){
        stop('The value of "pseudo.count" cannot be negative.')
    }

    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose)
    }
    # find assays ####
    if (length(assay.names) == 1 && assay.names == 'All') {
        assay.names = as.factor(names(assays(se.obj)))
    } else {
        assay.names = as.factor(unlist(assay.names))
    }
    # compute rle across all assays ####
    rle.all <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('-- Compute RLE on the ',
                                 x,
                                 ' assay.'),
                color = 'magenta',
                verbose = verbose)

            ### log transformation
            if (apply.log) {
                #message('Performing log + 1 transformation on the data')
                expr <- log2(assay(x = se.obj, i = x) + pseudo.count)
            } else expr <- assay(x = se.obj, i = x)
            rle.data <- expr - rowMedians(expr)
            rle.med <- colMedians(rle.data)
            rle.iqr <- colIQRs(rle.data)
            rle <- list(
                rle = rle.data,
                rle.med = rle.med,
                rle.iqr = rle.iqr)
        })
    names(rle.all) <- levels(assay.names)

    ## plot rle ####
    samples <- rle <- everything <- NULL
    plot.rle <- lapply(
        levels(assay.names),
        function(x) {
            tmp <- rle.all[[x]]$rle
            rle.data <- as.data.frame(tmp) %>%
                pivot_longer(everything(), names_to = 'samples', values_to = 'rle') %>%
                mutate(samples = factor(samples))
            p.rle <- ggplot(rle.data, aes(x = samples, y = rle)) +
                geom_boxplot2(width.errorbar = 0.01, outlier.alpha = 0.2) +
                ylab('RLE') +
                xlab('') +
                coord_cartesian(ylim = c(-6, 6)) +
                stat_summary(
                    geom = "crossbar",
                    width = 5,
                    fatten = 8,
                    color = "darkgreen",
                    fun.data = function(x) {
                        c(y = median(x),
                          ymin = median(x),
                          ymax = median(x))}) +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black'),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 8)) +
                ggtitle(paste(" RLE plot distribution of ", x, sep = "")) +
                geom_hline(yintercept = 0)
            p.rle
        })
    names(plot.rle) <- levels(assay.names)

    # save the results ####
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(message = '-- Save the RLE plot to the metadata of the SummarizedExperiment object.',
                            color = 'magenta',
                            verbose = verbose)
        ## For all assays
        for (x in levels(assay.names)) {
            ## Check if metadata plot already exist
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['plot']] <- list()
            }
            ## Check if metadata plot already exist
            if (!'plot' %in% names(se.obj@metadata)) {
                se.obj@metadata[['plot']] <- list()
            }
            ## Check if metadata plot already exist for this metric
            if (!'rle' %in% names(se.obj@metadata[['plot']])) {
                se.obj@metadata[['plot']][['rle']] <- list()
            }
            ## Check if metadata plot already exist for this metric
            if (!x %in% names(se.obj@metadata[['plot']])) {
                ## Save the new plot
                se.obj@metadata[['plot']][['rle']][[x]] <- list()
            }
            se.obj@metadata[['plot']][['rle']][[x]] <- plot.rle[[x]]
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
        return(plot = plot.rle)
    }

}
