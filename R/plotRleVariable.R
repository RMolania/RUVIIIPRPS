#' plot a variable against the medians and IQR of relative log expression (RLE) data

#' @author Ramyar Molania

#' @description
#' This function plots a variable against the medians and IQR of a relative log expression (RLE) data. Because of the
#' sensitivity of the RLE medians and IQR to unwanted variation, we  examine the relationships between RLE medians and
#' IQR with potential sources of unwanted variation. In the absence of any influence of unwanted variation in the data,
#' we should see no such associations.

#' @details
#' If the variable is categorical, the, boxplots of the RLE medians and IQR against the variable will be generated. If
#' the variable is continuous, scatter plots will be created.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to plot. By default all the assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. Indicates a name of the columns in the sample annotation of the SummarizedExperiment object.
#' The variable can be either categorical or continuous.
#' @param rle.data.type Symbol. Indicates which RLE data should be used for plotting. The options are 'rle.medians',
#' 'rle.iqr' or 'both'. If 'rle.medians' is selected, the RLE medians will be plotted against the variable. If 'rle.iqr',
#' is selected, the RLE IQRs will be plotted against the variable, and if 'both', both RLE medians and IQRs will be plotted
#' against the variable. The default is 'both'.
#' @param ylim.rle.med.plot Numeric. Indicates the ylim of the boxplot or scatter plots when the RLE medians are used. If
#' is NULL, the function will automatically find an suitable ylim for the plots.
#' @param ylim.rle.iqr.plot Numeric. Indicates the ylim of the boxplot or scatter plots when the RLE IQRs are used. If
#' is NULL, the function will automatically find an suitable ylim for the plots. The default is NULL.
#' @param points.size Numeric. Indicates the size of the points of the scatter plots. The default is 1.
#' @param plot.ncol Numeric. Indicates number of columns in the plot grid.
#' @param plot.nrow Numeric. Indicates number of rows in the plot grid.
#' @param plot.output Logical. If TRUE, individual RLE plot(s) will be printed while the function is running.
#' @param save.se.obj Logical. Indicates whether to save the plots in the metadata of the SummarizedExperiment object or
#' to output them as list. Default is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object that contains all the plot(s) in the metadata or a list that contains all the plot(s).

#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @importFrom SummarizedExperiment assays
#' @importFrom ggpubr ggarrange stat_cor
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
        plot.ncol = 3,
        plot.nrow = 3,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The plotRleVariable function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    }
    if(is.list(assay.names)){
        stop('The "assay.names" must be a vector of assay names(s) or "assay.names = all".')
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
    if (!is.null(ylim.rle.med.plot)) {
        if (length(ylim.rle.med.plot) != 2)
            stop('Please specify the "ylim.rle.med.plot" argument.')
    }
    if (!is.null(ylim.rle.iqr.plot)) {
        if (length(ylim.rle.iqr.plot) != 2)
            stop('Please specify the "ylim.rle.iqr.plot" argument.')
    }
    if (sum(is.na(se.obj@colData[[variable]])) > 0){
        stop(paste0('The "', variable, '" variable contains NA. ',
                    'Run the checkSeObj function with "remove.na = both"',
                    ', then re-run the computeRLE function.'))
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
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
            list(rle.medians = se.obj@metadata[['metric']][[x]][['RLE']][['rle.data']]$rle.med,
                 rle.iqrs = se.obj@metadata[['metric']][[x]][['RLE']][['rle.data']]$rle.iqr)
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
                    ggtitle(x) +
                    xlab(variable) +
                    ylab('RLE medians') +
                    geom_smooth(formula = y ~ x, method = 'lm') +
                    ggpubr::stat_cor(aes(label = ..r.label..), color = 'red', label.y = max(ylim.rle.med)) +
                    coord_cartesian(ylim = ylim.rle.med.plot) +
                    theme(panel.background = element_blank(),
                           axis.line = element_line(colour = 'black', linewidth = 1),
                           axis.title.x = element_text(size = 14),
                           axis.title.y = element_text(size = 14),
                           axis.text.x = element_text(size = 14),
                           axis.text.y = element_text(size = 10),
                           legend.position = 'bottom')
            } else{
                # boxplot ####
                      printColoredMessage(
                          message = paste0(
                              '-Generate a boxplot between the RLE medians of the ', x,
                              ' data and the variable:'),
                          color = 'blue',
                          verbose = verbose)
                p.rle <- ggplot(rle.med.data, aes(x = var, y = rle.medians)) +
                    geom_boxplot() +
                    ggtitle(x) +
                    xlab(variable) +
                    ylab('RLE medians') +
                    coord_cartesian(ylim = ylim.rle.med.plot) +
                    theme(panel.background = element_blank(),
                          axis.line = element_line(colour = 'black', linewidth = 1),
                          axis.title.x = element_text(size = 14),
                          axis.title.y = element_text(size = 14),
                          axis.text.x = element_text(size = 12, angle = 35, vjust = 1, hjust = 1),
                          axis.text.y = element_text(size = 10),
                          legend.position = 'bottom')
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
            nrow = plot.nrow,
            common.legend = TRUE)
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
                    message = paste0(
                        '-Generate a scatter plot between the RLE IQR of the ', x,
                        ' data and the variable:'),
                    color = 'blue',
                    verbose = verbose)
                p.rle <- ggplot(rle.med.data, aes(x = var, y = rle.iqr)) +
                    geom_point(size = points.size) +
                    ggtitle(x) +
                    xlab(variable) +
                    ylab('RLE IQRs') +
                    geom_smooth(formula = y ~ x, method = 'lm') +
                    ggpubr::stat_cor(aes(label = ..r.label..), color = 'red', label.y = max(ylim.rle.iqr)) +
                    coord_cartesian(ylim = ylim.rle.iqr.plot) +
                    theme(panel.background = element_blank(),
                          axis.line = element_line(colour = 'black', linewidth = 1),
                          axis.title.x = element_text(size = 14),
                          axis.title.y = element_text(size = 14),
                          axis.text.x = element_text(size = 14),
                          axis.text.y = element_text(size = 10),
                          legend.position = 'bottom')
            } else{
                # boxplot ####
                printColoredMessage(
                    message = paste0(
                        '-Generate a boxplot between the RLE IQR of the ', x,
                        ' data and the variable:'),
                    color = 'blue',
                    verbose = verbose)
                p.rle <- ggplot(rle.med.data, aes(x = var, y = rle.iqr)) +
                    geom_boxplot() +
                    ggtitle(x) +
                    xlab(variable) +
                    ylab('RLE IQRs') +
                    coord_cartesian(ylim = ylim.rle.iqr.plot) +
                    theme(panel.background = element_blank(),
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
        ## nrow and ncol of plot grid ####
        printColoredMessage(
            message = '-- Put together all the plots of the RLE medians and the variable:',
            color = 'magenta',
            verbose = verbose)
        overall.rle.iqr.var.plots <- ggpubr::ggarrange(
            plotlist = all.rle.iqr.var.plots,
            ncol = plot.ncol,
            nrow = plot.nrow,
            common.legend = TRUE)
        if (plot.output) print(overall.rle.iqr.var.plots)
    }
    # save the plots ####
    printColoredMessage(
        message = '-- Save all the RLE plots:',
        color = 'magenta',
        verbose = verbose)
    ## add plots to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
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
            verbose = verbose)
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
                se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']][['RleMedians']] <- overall.rle.med.var.plots

            } else if (rle.data.type == 'rle.iqr'){
                if(!'RleIqr' %in% se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']]){
                    se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']][['RleIqr']] <- list()
                }
                se.obj@metadata[['plot']][['RLE']][[variable]][['RleVarPlot']][['RleIqr']] <- overall.rle.iqr.var.plots
            }
            printColoredMessage(
                message = paste0('The RLE plots of all assays are saved to metadata@plot'),
                color = 'blue',
                verbose = verbose
            )
        }
        printColoredMessage(message = '------------The plotRleVariable function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj = se.obj)
        ## save the plots as list ####
    } else if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '-- All the plots are outputed as list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The plotRleVariable function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) == 1){
            if(rle.data.type == 'both'){
                rle.var.plots <- list(
                    all.rle.med.var.plots = all.rle.med.var.plots,
                    all.rle.iqr.var.plots = all.rle.iqr.var.plots)
            } else if (rle.data.type == 'rle.medians'){
                rle.var.plots <- list(all.rle.med.var.plots = all.rle.med.var.plots)
            } else if (rle.data.type == 'rle.iqrs'){
                rle.var.plots <- list(all.rle.iqr.var.plots = all.rle.iqr.var.plots)
            }
        } else if (length(assay.names) > 1){
            if(rle.data.type == 'both'){
                rle.var.plots <- list(
                    all.rle.med.var.plots = all.rle.med.var.plots,
                    all.rle.iqr.var.plots = all.rle.iqr.var.plots,
                    overall.rle.med.var.plots = overall.rle.med.var.plots,
                    overall.rle.iqr.var.plots = overall.rle.iqr.var.plots)
            } else if (rle.data.type == 'rle.medians'){
                rle.var.plots <- list(
                    all.rle.med.var.plots = all.rle.med.var.plots,
                    overall.rle.med.var.plots = overall.rle.med.var.plots)
            } else if (rle.data.type == 'rle.iqrs'){
                rle.var.plots <- list(
                    all.rle.iqr.var.plots = all.rle.iqr.var.plots,
                    overall.rle.iqr.var.plots = overall.rle.iqr.var.plots)
            }
        }
    }
}
