#' is used to plot relative log expression (RLE) of RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This functions generate a boxplot of individual RLE data. The RLE data can be obtained using the computeRLE function.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to plot RLE. By default all the assays of the SummarizedExperiment object
#' will be selected.
#' @param variable Symbol. Indicates a name of column in the sample annotation for the SummarizedExperiment object to color
#' the boxplots of the RLE plots. The variable should be a categorical variable. The default is NULL.
#' @param rle.data.type TTTTT
#' @param ylim.rle.med.plot TTTTTTTT
#' @param ylim.rle.iqr.plot Numeric. A vector of two values to specify the ylim of the RLE plot.
#' @param points.size Numeric. Indicates the size of the points of the RLE medians in the plot.
#' @param plot.ncol Numeric. Indicates number of columns in the plot grid.
#' @param plot.output Logical. If TRUE, individual RLE plot(s) will be printed.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment object or
#' to output the result. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list that contains all the RLE plot(s).

#' @importFrom SummarizedExperiment assays
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#' @export

plotRleVariable <- function(
        se.obj,
        assay.names = "all",
        variable,
        rle.data.type = 'both',
        ylim.rle.med.plot = NULL,
        ylim.rle.iqr.plot = NULL,
        points.size = 1,
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
    if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    }
    if(length(rle.data.type) > 1){
        stop('The "rle.data.type" must be one of the "rle.medians", "rle.iqr" or "both".')
    }
    if (!rle.data.type %in% c('rle.medians', 'rle.iqr', 'both')) {
        stop('The "rle.data.type" must be one of the "rle.medians", "rle.iqr" or "both".')
    }
    if (is.null(ylim.rle.med.plot) | length(ylim.rle.med.plot) != 2) {
        stop('Please specify the "ylim.rle.med.plot" argument.')
    }
    if (is.null(ylim.rle.iqr.plot) | length(ylim.rle.iqr.plot) != 2) {
        stop('Please specify the "ylim.rle.iqr.plot" argument.')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- as.factor(assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # obtain rle data ####
    printColoredMessage(
        message = paste0('-- Obtain the RLE medians and IQR data from the SummarizedExperiment object:'),
        color = 'magenta',
        verbose = verbose
    )
    all.rle.data <- lapply(
        levels(assay.names),
        function(x) {
            if (!'RLE' %in% names(se.obj@metadata[['metric']][[x]])) {
                stop(paste0(
                    'The RLE data has not been computed yet for the ',
                    x,
                    ' assay. Please run the computeRLE function first.'))
            }
            if (!'rle.data' %in% names(se.obj@metadata[['metric']][[x]][['RLE']])) {
                stop(paste0(
                    'The "rle.data" cannot be found for the ',
                    x,
                    ' assay. Please run the computeRLE function first.'))
            }
            list(rle.medians = se.obj@metadata[['metric']][[x]][['RLE']][['rle.data']]$rle.medians,
                 rle.iqrs = se.obj@metadata[['metric']][[x]][['RLE']][['rle.data']]$rle.iqrs)
        })
    names(all.rle.data) <- levels(assay.names)

    # specify ylim for the RLE plots ####
    if(is.null(ylim.rle.med.plot)){
        ylim.rle.med <- unlist(lapply(
            levels(assay.names),
            function(x){
                c(min(all.rle.data[[x]]$rle.medians),
                  max(all.rle.data[[x]]$rle.medians))
            }))
        ylim.rle.med.plot <- c(min(ylim.rle.med), max(ylim.rle.med))
    }
    if(is.null(ylim.rle.iqr.plot)){
        ylim.rle.iqr <- unlist(lapply(
            levels(assay.names),
            function(x){
                c(min(all.rle.data[[x]]$rle.iqrs),
                  max(all.rle.data[[x]]$rle.iqrs))
            }))
        ylim.rle.iqr.plot <- c(min(ylim.rle.iqr), max(ylim.rle.iqr))
    }

    # generate the plots ####
    printColoredMessage(
        message = paste0('-- Generate the plots with the RLE medians:'),
        color = 'magenta',
        verbose = verbose
    )
    samples <- rle <- everything <- sample <- rle.medians <- NULL
    all.rle.med.var.plots <- lapply(
        levels(assay.names),
        function(x) {
            rle.med.data <- data.frame(
                rle.medians = all.rle.data[[x]]$rle.medians,
                var = colData(se.obj)[[variable]])
            if(class(colData(se.obj)[[variable]]) %in% c('numeric', 'integr')){
                # scatter plot ####
                printColoredMessage(
                    message = paste0('-Generate a scatter plot between the RLE medians of the ',
                                     x,
                                     ' data and the variable:'),
                    color = 'blue',
                    verbose = verbose)
                p.rle <- ggplot(rle.med.data, aes(x = var, y = rle.medians)) +
                    geom_point() +
                    ylab('RLE medians') +
                    ggtitle(x) +
                    xlab(variable) +
                    coord_cartesian(ylim = ylim.rle.med.plot) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        legend.position = 'bottom',
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 10)
                    )
            } else{
                # boxplot ####
                      printColoredMessage(
                          message = paste0('-Generate a boxplot between the RLE medians of the ',
                                           x,
                                           ' data and the variable:'),
                          color = 'blue',
                          verbose = verbose)
                p.rle <- ggplot(rle.med.data, aes(x = var, y = rle.medians)) +
                    geom_boxplot() +
                    ylab('RLE medians') +
                    ggtitle(x) +
                    xlab(variable) +
                    coord_cartesian(ylim = ylim.rle.med.plot) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        legend.position = 'bottom',
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 12, angle = 35, vjust = 1, hjust = 1),
                        axis.text.y = element_text(size = 10)
                    )
            }
            if (plot.output) print(p.rle)
            p.rle
        })
    names(all.rle.med.var.plots) <- levels(assay.names)

    # generate the overall RLE plots ####
    if(length(assay.names) > 1){
        printColoredMessage(
            message = '-- Put together all the plots of the RLE medians and the variable:',
            color = 'magenta',
            verbose = verbose
        )
        overall.rle.med.var.plots <- ggpubr::ggarrange(
            plotlist = all.rle.med.var.plots,
            ncol = plot.ncol,
            common.legend = T)
        if (plot.output) print(overall.rle.med.var.plots)
    }

    # generate the RLE IQR plots ####
    printColoredMessage(
        message = '-- Generate the plots with the RLE IQRs:',
        color = 'magenta',
        verbose = verbose
        )
    samples <- rle <- everything <- sample <- rle.iqr <- NULL
    all.rle.iqr.var.plots <- lapply(
        levels(assay.names),
        function(x) {
            rle.med.data <- data.frame(
                rle.iqr = all.rle.data[[x]]$rle.iqrs,
                var = colData(se.obj)[[variable]])
            if(class(colData(se.obj)[[variable]]) %in% c('numeric', 'integr')){
                # scatter plot ####
                printColoredMessage(
                    message = paste0('-Generate a scatter plot between the RLE IQR of the ',
                                     x,
                                     ' data and the variable:'),
                    color = 'blue',
                    verbose = verbose)
                p.rle <- ggplot(rle.med.data, aes(x = var, y = rle.iqr)) +
                    geom_point(size = points.size) +
                    ggtitle(x) +
                    xlab(variable) +
                    ylab('RLE IQRs') +
                    coord_cartesian(ylim = ylim.rle.iqr.plot) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        legend.position = 'bottom',
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 10))
            } else{
                # boxplot ####
                printColoredMessage(
                    message = paste0('-Generate a boxplot between the RLE IQR of the ',
                                     x,
                                     ' data and the variable:'),
                    color = 'blue',
                    verbose = verbose)
                p.rle <- ggplot(rle.med.data, aes(x = var, y = rle.iqr)) +
                    geom_boxplot() +
                    xlab(variable) +
                    ylab('RLE IQRs') +
                    ggtitle(x) +
                    coord_cartesian(ylim = ylim.rle.iqr.plot) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        legend.position = 'bottom',
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 10))
            }
            if (plot.output) print(p.rle)
            p.rle
        })
    names(all.rle.iqr.var.plots) <- levels(assay.names)

    # generate the overall RLE plots ####
    if(length(assay.names) > 1){
        printColoredMessage(
            message = '-- Put together all the plots of the RLE medians and the variable:',
            color = 'magenta',
            verbose = verbose
        )
        overall.rle.iqr.var.plots <- ggpubr::ggarrange(
            plotlist = all.rle.iqr.var.plots,
            ncol = plot.ncol,
            common.legend = T)
        if (plot.output) print(overall.rle.iqr.var.plots)
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
            if (!variable %in% names(se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']])){
                se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]] <- list()
            }
            if (!'RleVarPlot' %in% names(se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]])){
                se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']] <- list()
            }
            if(rle.data.type == 'both'){
                if (!'RleMed' %in% names(se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']])){
                    se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']][['RleMed']] <- list()
                }
                se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']][['RleMed']] <- all.rle.med.var.plots[[x]]
                if (!'RleIqr' %in% names(se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']])){
                    se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']][['RleIqr']] <- list()
                }
                se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']][['RleIqr']] <- all.rle.iqr.var.plots[[x]]

            }else if (rle.data.type == 'rle.medians') {
                if (!'RleMed' %in% names(se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']])){
                    se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']][['RleMed']] <- list()
                }
                se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']][['RleMed']] <- all.rle.med.var.plots[[x]]
            } else if (rle.data.type == 'rle.iqr'){
                if (!'RleIqr' %in% names(se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']])){
                    se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']][['RleIqr']] <- list()
                }
                se.obj@metadata[['metric']][[x]][['RLE']][['rle.plot']][[variable]][['RleVarPlot']][['RleIqr']] <- all.rle.iqr.var.plots[[x]]
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
            if (!variable %in%  names(se.obj@metadata[['plot']][['RLE']])) {
                se.obj@metadata[['plot']][['RLE']][['variable']] <- list()
            }
            if (!'RleVarPlot' %in%  names(se.obj@metadata[['plot']][['RLE']][[variable]])) {
                se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']] <- list()
            }
            if(rle.data.type == 'both'){
                if(!'RleMedians' %in% se.obj@metadata[['plot']][['RLE']][['variable']][['RleVarPlot']]){
                    se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']][['RleMedians']] <- list()
                }
                se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']][['RleMedians']] <- overall.rle.med.var.plots
                if(!'RleIqr' %in% se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']]){
                    se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']][['RleIqr']] <- list()
                }
                se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']][['RleMedians']] <- overall.rle.iqr.var.plots
            } else if (rle.data.type == 'rle.medians'){
                if(!'RleMedians' %in% se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']]){
                    se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']][['RleMedians']] <- list()
                }
                se.obj@metadata[['plot']][['RLE']][['variable']][['RleVarPlot']][['RleMedians']] <- overall.rle.med.var.plots

            } else if (rle.data.type == 'rle.iqr'){
                if(!'RleIqr' %in% se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']]){
                    se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']][['RleIqr']] <- list()
                }
                se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']][['RleMedians']] <- overall.rle.iqr.var.plots
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
            message = '-- Save all the RLE plots as list.',
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
