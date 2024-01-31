#' generate boxplot of the relative log expression (RLE) distributions of RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function generates a boxplot of individual RLE data. The RLE data can be obtained using the computeRLE function.
#' We refer to the computeRLE function for more detail about RLE.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to calculate RLE data, medians and interquartiles. The default is "all, which indicates all
#' the assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. Indicates a name of the column in the sample annotation of the SummarizedExperiment object.
#' The interquartile ranges of the RLE boxplot will be colored based on the variable. The variable must be a categorical
#' variable. The default is NULL.
#' @param ylim.rle.plot Numeric. A vector of two values to specify the ylim of the RLE plot(s). If is NULL, the function
#' uses the minimum and maximum interquartile ranges of all the RLE data as ylim. The default is NULL.
#' @param iqr.width Numeric. Indicates the width size of RLE interquartile ranges in the plot. The default is 1.
#' @param median.points.size Numeric. Indicates the size of the points of the RLE medians in the boxplot(s). The default is 1.
#' @param median.points.color Symbol. Indicates the color of the points of the RLE medians in the boxplot(s).
#' @param geom.hline.color Symbol. Indicates the color of the horizontal line (geom.hline) across 0 in the RLE boxplot(s).
#' This line helps to see the deviation of the RLE medians of the RLE boxplot(s).
#' @param plot.ncol Numeric. Indicates number of columns in the plot grid. When the number of selected assay is more than
#' 1, the function puts all the RLE boxplots in one grid.
#' @param plot.output Logical. If TRUE, the individual RLE plot(s) will be printed while functions is running.
#' @param save.se.obj Logical. Indicates whether to save the RLE plots in the metadata of the SummarizedExperiment object
#'  or to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list that contains all the RLE plot(s).

#' @importFrom SummarizedExperiment assays
#' @importFrom matrixStats colQuantiles
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#' @export

plotRLE <- function(
        se.obj,
        assay.names = "all",
        variable = NULL,
        ylim.rle.plot = NULL,
        iqr.width = 1,
        median.points.size = 1,
        median.points.color = 'black',
        geom.hline.color = 'cyan',
        plot.ncol = 1,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The plotRLE function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    }
    if (!is.null(variable)){
        if (length(variable) > 1){
            stop('The "variable" must contain only one variable.')
        }
        if (!class(colData(se.obj)[[variable]]) %in% c('factor', 'character')){
            stop('The "variable" must be a categorical variable.')
        }
        if (sum(is.na(se.obj@colData[[variable]])) > 0){
            stop(paste0('The "', variable, '" variable contains NA. ',
                        'Run the checkSeObj function with "remove.na = both"',
                        ', then re-run the computeRLE function.'))
        }
    }
    if (!is.null(ylim.rle.plot) & !length(ylim.rle.plot) != 2) {
        stop('Please specify the approprate "ylim.rle.plot" argument.')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(x = names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # select colors ####
    if(!is.null(variable)){
        variable.data <- as.factor(colData(se.obj)[[variable]])
        if(length(levels(variable)) < 9 ){
            rle.plot.colors <- RColorBrewer::brewer.pal(8, 'Dark2')[1:length(levels(variable.data))]
            names(rle.plot.colors) <- levels(variable.data)
        } else {
            colfunc <- grDevices::colorRampPalette( RColorBrewer::brewer.pal(8, 'Dark2'))
            rle.plot.colors <- colfunc(n = length(levels(assay.names)))
            names(rle.plot.colors) <- levels(variable.data)
        }
    }

    # obtain rle data ####
    printColoredMessage(
        message = paste0('-- Obtain the RLE data from the SummarizedExperiment object:'),
        color = 'magenta',
        verbose = verbose
    )
    all.rle.data <- lapply(
        levels(assay.names),
        function(x) {
            if (!'RLE' %in% names(se.obj@metadata[['metric']][[x]])) {
                stop(paste0(
                        'The RLE data has not been computed yet for the ', x,
                        ' data. Please run the computeRLE function first.'))
            }
            if (!'rle.data' %in% names(se.obj@metadata[['metric']][[x]][['RLE']])) {
                stop(paste0(
                        'The "rle.data" cannot be found for the ', x,
                        ' data. Please run the computeRLE function first.'))
            }
            if (!'rle.data' %in% names(se.obj@metadata[['metric']][[x]][['RLE']][['rle.data']])) {
                stop(paste0(
                        'The "rle.data" cannot be found for the ', x,
                        ' data. Please run the computeRLE function with the "outputs.to.return" equal to "all" or "rle.data".'))
            }
            se.obj@metadata[['metric']][[x]][['RLE']][['rle.data']]$rle.data
        })
    names(all.rle.data) <- levels(assay.names)

    # specify ylim for the RLE plots ####
    if(is.null(ylim.rle.plot)){
        ylim.rle.plot <- abs(unlist(lapply(
            levels(assay.names),
            function(x){
                samples.quantiles <- matrixStats::colQuantiles(
                    x = all.rle.data[[x]],
                    probs = c(0.2, 0.8))
                c(max(samples.quantiles), min(samples.quantiles))
            })))
        ylim.rle.plot <- c(-max(ylim.rle.plot), max(ylim.rle.plot))
    }

    # generate the RLE plots ####
    printColoredMessage(
        message = '-- Generate the RLE plots:',
        color = 'magenta',
        verbose = verbose
    )
    rle <- sample <- NULL
    all.rle.plots <- lapply(
        levels(assay.names),
        function(x) {
            rle.data <- all.rle.data[[x]]
            rle.med.var <- stats::mad(matrixStats::colMedians(rle.data))
            rle.iqr.var <- stats::mad(matrixStats::colIQRs(rle.data))
            samples.quantiles <- matrixStats::colQuantiles(
                x = rle.data,
                probs = seq(from = 0, to = 1, by = 0.25))
            samples.quantiles <- as.data.frame(samples.quantiles[, c(2:4)])
            colnames(samples.quantiles) <- c('per25', 'medians', 'per75')
            samples.quantiles$sample <- paste0('Sam', 1:ncol(rle.data))
            if(!is.null(variable)){
                # colored RLE plots ####
                samples.quantiles$variable <- variable.data
                samples.quantiles <- tidyr::pivot_longer(
                    data = samples.quantiles,
                    -c(sample, variable),
                    values_to = 'rle',
                    names_to = 'range')
                samples.quantiles$sample <- factor(samples.quantiles$sample, levels = paste0('Sam', 1:ncol(rle.data)))
                p.rle <- ggplot(samples.quantiles, aes(x = sample, y = rle, group = sample)) +
                    geom_line(aes(color = variable), linewidth = 1) +
                    geom_point(data = samples.quantiles[samples.quantiles$range == 'medians',],
                               aes(group = range),
                               size = median.points.size ,
                               colour = median.points.color) +
                    scale_color_manual(name = 'Groups:', values = rle.plot.colors) +
                    ylab('RLE') +
                    xlab('Samples') +
                    ggtitle(paste0('Data:', x, ',VarRleMed:', rle.med.var, 'VarRleIqr:', rle.iqr.var)) +
                    coord_cartesian(ylim = ylim.rle.plot) +
                    geom_hline(yintercept = 0, colour = geom.hline.color) +
                    theme(
                        panel.background = element_blank(),
                        plot.title = element_text(size = 12),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size = 9),
                        axis.ticks.x = element_blank(),
                        legend.position = 'bottom')
            } else if (is.null(variable)){
                # general RLE plots ####
                samples.quantiles <- tidyr::pivot_longer(
                    data = samples.quantiles,
                    -sample,
                    values_to = 'rle',
                    names_to = 'range')
                samples.quantiles$sample <- factor(samples.quantiles$sample, levels = paste0('Sam', 1:ncol(rle.data)))
                p.rle <- ggplot(samples.quantiles, aes(x = sample, y = rle, group = sample)) +
                    geom_line(linewidth = 1) +
                    geom_point(data = samples.quantiles[samples.quantiles$range == 'medians',],
                               aes(group = range),
                               size = median.points.size ,
                               colour = median.points.color) +
                    ylab('RLE') +
                    xlab('Samples') +
                    ggtitle(paste0('Data:', x, ',VarRleMed:', rle.med.var, 'VarRleIqr:', rle.iqr.var)) +
                    coord_cartesian(ylim = ylim.rle.plot) +
                    geom_hline(yintercept = 0, colour = geom.hline.color) +
                    theme(
                        panel.background = element_blank(),
                        plot.title = element_text(size = 12),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size = 9),
                        axis.ticks.x = element_blank())
            }
            if (plot.output) print(p.rle)
            p.rle
        })
    names(all.rle.plots) <- levels(assay.names)

    # generate the overall RLE plots ####
    if(length(assay.names) > 1){
        overall.rle.plot <- ggpubr::ggarrange(
            plotlist = all.rle.plots,
            ncol = plot.ncol,
            legend = "bottom",
            common.legend = TRUE)
        if (plot.output) print(overall.rle.plot)
    }

    # save the plots ####
    printColoredMessage(
        message = '-- Save all the RLE plots:',
        color = 'magenta',
        verbose = verbose)
    ## add plots to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = '-- Save all the RLE plots to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        ## add RLE plots of individual assays ####
        ### check the metadata of the SummarizedExperiment object ####
        ### for each assays ####
        for (x in levels(assay.names)) {
            if (!'rle.plot' %in% names(se.obj@metadata[['metric']][[x]][['RLE']])) {
                se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']] <- list()
            }
            if(!is.null(variable)){
                if (!variable %in% names(se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']]) ){
                    se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]] <- list()
                }
                if (!'ColoredRLE' %in% names(se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]])){
                    se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['ColoredRLE']] <- list()
                }
                se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['ColoredRLE']] <- all.rle.plots[[x]]
            } else if (is.null(variable)){
                if (!'GeneralRLE' %in% names(se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']])){
                    se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][['GeneralRLE']] <- list()
                }
                se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][['GeneralRLE']] <- all.rle.plots[[x]]
            }
        }
        printColoredMessage(
            message = paste0(
                'The RLE plots of individual assays are saved to metadata@metric'),
            color = 'blue',
            verbose = verbose
        )

        ## add overall RLE plots of all assays ####
        if(length(assay.names) > 1){
            if (!'plot' %in%  names(se.obj@metadata)) {
                se.obj@metadata[['plot']] <- list()
            }
            if (!'RLE' %in%  names(se.obj@metadata[['plot']])) {
                se.obj@metadata[['plot']][['RLE']] <- list()
            }
            if(!is.null(variable)){
                if(variable %in% names(se.obj@metadata[['plot']][['RLE']])){
                    se.obj@metadata[['plot']][['RLE']][[variable]] <- list()
                }
                if(!'ColoredRLE' %in% names(se.obj@metadata[['plot']][['RLE']][[variable]])){
                    se.obj@metadata[['plot']][['RLE']][[variable]][['ColoredRLE']] <- list()
                }
                se.obj@metadata[['plot']][['RLE']][[variable]][['ColoredRLE']] <- overall.rle.plot
            } else {
                if(!'GeneralRLE' %in% names(se.obj@metadata[['plot']][['RLE']])){
                    se.obj@metadata[['plot']][['RLE']][['GeneralRLE']] <- list()
                }
                se.obj@metadata[['plot']][['RLE']][['GeneralRLE']] <- overall.rle.plot
            }
            printColoredMessage(
                message = paste0('The RLE plots of all assays are saved to metadata@plot'),
                color = 'blue',
                verbose = verbose
            )
        }
        printColoredMessage(
            message = '------------The plotRLE function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj = se.obj)
    } else if (save.se.obj == FALSE) {
        printColoredMessage(
            message = '-- All the RLE plots are outputed as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The plotRLE function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) == 1){
            return(rle.plots = list(all.rle.plots = all.rle.plots))
        } else{
            return(rle.plots = list(
                all.rle.plots = all.rle.plots,
                overall.rle.plot = overall.rle.plot))
        }
    }
}
