#' Generates scatter and boxplot plot of principal components.

#' @description
#' This functions generates scatter and boxplot the first principal components of the assay(s) of a SummarizedExperiment
#' object. The function can generate pairwise scatter plots of the first principal components colored by a categorical
#' variable or creates boxplots of each PC across the variable. Boxplots of principal components become beneficial when
#' the levels of a categorical variable are numerous, making visualization through colored scatter plots challenging. If
#' a continuous variable is provided, the function creates a scatter plot of each PC against the variable.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or a vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to PC. The default is "all, which indicates all the assays of the SummarizedExperiment
#' object will be selected.
#' @param variable Symbol. Indicates a name of the column in the sample annotation of the SummarizedExperiment object.
#' The variable can be either a categorical or continuous. If a continuous variable is provided, the function creates
#' a scatter plot of each PC against the variable. If a categorical variable is given, the function can generate pairwise
#' scatter plots of the first principal components colored by a categorical variable or creates boxplots of each PC across
#' the variable.
#' @param fast.pca Logical. Indicates whether to use the computed fast PC ro not. The default is 'TRUE'. We refe to the
#' computePCA function for more details.
#' @param nb.pcs Numeric. The number of first PCs to be used for plotting. The default is set to 3.It's important to note
#' that the variation of PCs for a portion of all selected PCs will be based solely on those selected PCs.
#' @param plot.type Symbol. Indicates which plot type should be used. Options are: 'scatter' and 'boxplot'. The default is
#' set to 'scatter'. Please note, the 'plot.type' cannot be set to 'scatter' of a categorical variable is provided.
#' @param points.size Numeric. Indicates the size of points in the scatter PCA plots.
#' @param stroke.color Symbol. Indicates the color of the stroke of the points in the scatter PCA plots.
#' @param stroke.size Numeric. Indicates the size of the stroke of the points in the scatter PCA plots.
#' @param points.alpha Numeric. Indicates the transparency of the points in the scatter PCA plots.
#' @param densities.alpha Numeric. Indicates the transparency of the densities in the scatter PCA plots.
#' @param plot.ncol Numeric. Indicates number of columns in the plot grid. When the number of selected assay is more than
#' 1, the function puts all the PCA plots in one grid.
#' @param plot.nrow Numeric. Indicates number of rows in the plot grid. When the number of selected assay is more than
#' 3, the function puts all the PCA plots in one grid.
#' @param plot.output Logical. If TRUE, the individual PCA plot(s) will be printed while functions is running.
#' @param save.se.obj Logical. Indicates whether to save the plots in the metadata of the SummarizedExperiment object
#' or to output the results as list. By default it is set to 'TRUE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object that contains all the PCA plot(s) in the metadata or a list that contains all
#' the PCA plot(s).

#' @importFrom ggpubr ggarrange theme_pubr
#' @importFrom patchwork plot_spacer plot_layout
#' @importFrom tidyr pivot_longer
#' @importFrom utils combn
#' @import ggplot2
#' @export

plotPCA <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        nb.pcs = 3,
        plot.type = 'scatter',
        points.size = 1,
        stroke.color = 'gray30',
        stroke.size = .2,
        points.alpha = .5,
        densities.alpha = .5,
        plot.ncol = NULL,
        plot.nrow = 3,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The plotPCA function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (is.null(assay.names)){
        stop('The "assay.names" cannot be empty.')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    } else if (length(variable) > 1){
        stop('The "variable" must contain only one variable.')
    } else if (!variable %in% colnames(colData(se.obj))) {
        stop('The "variable" cannot be found in the SummarizedExperiment object.')
    } else if (is.null(nb.pcs)){
        stop('The "nb.pcs" cannot be empty.')
    } else if (!plot.type %in% c('scatter', 'boxplot')){
        stop('The "plot.type" must be one of "scatter" or "boxplot".')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # obtain PCs from the SummarizedExperiment object ####
    printColoredMessage(
        message = paste0('Obtain the first ', nb.pcs, ' from the SummarizedExperiment object.'),
        color = 'magenta',
        verbose = verbose)
    pc.var <- NULL
    all.pca.data <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('Obtain the first ', nb.pcs, ' PCs of ', x, ' data.'),
                color = 'blue',
                verbose = verbose)
            if (isTRUE(fast.pca)) {
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop(paste0(
                        'To plot the PCA, the fast PCA must be computed first on the assay ', x,
                        '. Please run the computePCA function first.'))
                pca.data <- se.obj@metadata[['metric']][[x]][['fastPCA']]$svd$u
                pc.var <- se.obj@metadata[['metric']][[x]][['fastPCA']]$percentage.variation
            } else if (isFALSE(fast.pca)) {
                if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop(paste0(
                        'To plot the PCA, the PCA must be computed first on the assay ', x,
                        '. Please run the computePCA function first.'))
                pca.data <- se.obj@metadata[['metric']][[x]][['PCA']]$svd$u[, 1:nb.pcs]
                pc.var <- se.obj@metadata[['metric']][[x]][['PCA']]$percentage.variation
            }
            if(ncol(pca.data) < nb.pcs){
                printColoredMessage(
                    message = paste0('The number of PCs of the assay', x, 'are ', ncol(pca.data), '.'),
                    color = 'blue',
                    verbose = verbose)
                stop(paste0('The number of PCs of the assay ', x, ' are less than', nb.pcs, '.',
                            'Re-run the computePCA function with nb.pcs = ', nb.pcs, '.'))
            }
            return(list(pca.data = pca.data, pc.var = pc.var))
        })
    names(all.pca.data) <- levels(assay.names)

    # plot PCs ####
    ## categorical variable ####
    if(class(colData(se.obj)[[variable]]) %in% c('character', 'factor')){
        ### scatter plots for individual assays  ####
        if(plot.type == 'scatter'){
            printColoredMessage(
                message = '-- Creat scatter PCA plots for the individual assay(s):',
                color = 'magenta',
                verbose = verbose)
            all.scat.pca.plots <- lapply(
                levels(assay.names),
                function(x) {
                    printColoredMessage(
                        message = paste0('PCA plots of for ', x, ' data.'),
                        color = 'blue',
                        verbose = verbose)
                    pca.data <- all.pca.data[[x]]$pca.data[ , seq_len(nb.pcs)]
                    pair.pcs <- combn(ncol(pca.data), 2)
                    var <- colData(se.obj)[, variable]
                    printColoredMessage(
                        message = paste0('-Create all possible pairwise scatter plots of the first ', nb.pcs, ' PCs.'),
                        color = 'blue',
                        verbose = verbose)
                    p.per.data <- lapply(
                        1:ncol(pair.pcs),
                        function(i) {
                            plot1 <- ggplot(mapping = aes(x = pca.data[, pair.pcs[1, i]], y = pca.data[, pair.pcs[2, i]])) +
                                geom_point(
                                    aes(fill = se.obj@colData[, variable]),
                                    color = stroke.color,
                                    pch = 21,
                                    stroke = stroke.size,
                                    size = points.size,
                                    alpha = points.alpha) +
                                scale_x_continuous(
                                    name = paste0('PC', pair.pcs[1, i], ' (', all.pca.data[[x]]$pc.var[pair.pcs[2, i]], '%)'),
                                    breaks = scales::pretty_breaks(n = 5)) +
                                scale_y_continuous(
                                    name = paste0('PC', pair.pcs[2, i], ' (', all.pca.data[[x]]$pc.var[pair.pcs[1, i]], '%)'),
                                    breaks = scales::pretty_breaks(n = 5)) +
                                ggtitle(x) +
                                theme_pubr() +
                                theme(
                                    legend.background = element_blank(),
                                    legend.text = element_text(size = 12),
                                    legend.title = element_text(size = 14),
                                    legend.key = element_blank(),
                                    axis.text.x = element_text(size = 8),
                                    axis.text.y = element_text(size = 8),
                                    axis.title.x = element_text(size = 10),
                                    axis.title.y = element_text(size = 10),
                                    legend.position = "none") +
                                scale_fill_discrete(name = variable) +
                                guides(fill = guide_legend(override.aes = list(size = 3, shape = 21)))
                            dense.x <-
                                ggplot(mapping = aes(x = pca.data[, pair.pcs[1, i]], fill = se.obj@colData[, variable])) +
                                geom_density(alpha = densities.alpha) +
                                theme_void() +
                                theme(
                                    legend.position = "none",
                                    legend.text = element_text(size = 12),
                                    legend.title = element_text(size = 14)) +
                                scale_fill_discrete(name = variable) +
                                guides(fill = guide_legend(override.aes = list(size = 3)))

                            dense.y <- ggplot(mapping = aes(x = pca.data[, pair.pcs[2, i]], fill = se.obj@colData[, variable])) +
                                geom_density(alpha = densities.alpha) +
                                theme_void() +
                                theme(
                                    legend.position = "none",
                                    legend.text = element_text(size = 12),
                                    legend.title = element_text(size = 14)) +
                                coord_flip() +
                                scale_fill_discrete(name = variable) +
                                guides(fill = guide_legend(override.aes = list(size = 3)))

                            dense.x +
                                plot_spacer() +
                                plot1 + dense.y +
                                plot_layout(
                                    ncol = 2,
                                    nrow = 2,
                                    widths = c(4, 1),
                                    heights = c(1, 4))
                        })
                    p.per.data
                })
            names(all.scat.pca.plots) <- levels(assay.names)

            ## overall scatter plots for all assays  ####
            printColoredMessage(message = '-- Put all the plots togather:',
                                color = 'magenta',
                                verbose = verbose)
            all.scat.pca.plots.assay <- lapply(
                levels(assay.names),
                function(x) {
                    ggarrange(
                        plotlist = all.scat.pca.plots[[x]],
                        common.legend = TRUE,
                        legend = "bottom",
                        nrow = 1)
                })
            names(all.scat.pca.plots.assay) <- levels(assay.names)
            if(length(assay.names) > 1){
                p <- all.scat.pca.plots[[levels(assay.names)[1]]]
                if (length(assay.names) > 1) {
                    for (n in 2:length(assay.names))
                        p <- c(p, all.scat.pca.plots[[levels(assay.names)[n]]])
                }
                if(is.null(plot.ncol)) plot.ncol <- ncol(combn(length(seq_len(nb.pcs)), 2))
                overall.scat.pca.plot <- ggarrange(
                    plotlist = p,
                    common.legend = TRUE,
                    legend = "bottom",
                    nrow = plot.nrow,
                    ncol = plot.ncol
                    )
                if(isTRUE(plot.output)) print(overall.scat.pca.plot)
            }
        } else{
            printColoredMessage(
                message = '-- Creat boxplot of PCA plots for individual assays:',
                color = 'magenta',
                verbose = verbose)
            all.boxplot.pca.plots <- lapply(
                levels(assay.names),
                function(x) {
                    printColoredMessage(
                        message = paste0('PCA plots of for ', x, ' data.'),
                        color = 'blue',
                        verbose = verbose)
                    var <- PC <- NULL
                    pca.data <- as.data.frame(all.pca.data[[x]]$pca.data[ , seq_len(nb.pcs)])
                    colnames(pca.data) <- paste0('PC', seq_len(nb.pcs), ' (',all.pca.data[[x]]$pc.var[seq_len(nb.pcs)], '%)')
                    pca.data$var <- colData(se.obj)[, variable]
                    printColoredMessage(
                        message = paste0('-Create all possible boxplots of the first ', nb.pcs, ' PCs.'),
                        color = 'blue',
                        verbose = verbose)
                    pca.data <- tidyr::pivot_longer(
                        data = pca.data,
                        -var,
                        names_to = 'PCs',
                        values_to = 'PC')
                    pca.data <- as.data.frame(pca.data)
                    plot.p <- ggplot(pca.data, aes(x = var, y =  PC, fill = var)) +
                        geom_boxplot() +
                        facet_wrap(~PCs, scales = 'free') +
                        xlab(variable) +
                        ylab('PC') +
                        ggtitle(x) +
                        theme(
                            panel.background = element_blank(),
                            axis.line = element_line(colour = 'black', size = 1),
                            axis.title.x = element_text(size = 10),
                            axis.title.y = element_text(size = 10),
                            plot.title = element_text(size = 12),
                            axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1),
                            axis.text.y = element_text(size = 10),
                            strip.text.x = element_text(size = 12),
                            legend.position = 'none')
                    plot.p
                })
            names(all.boxplot.pca.plots) <- levels(assay.names)
            if(length(assay.names) > 1){
                overall.boxplot.pca.plot <- ggarrange(
                    plotlist = all.boxplot.pca.plots,
                    nrow = length(levels(assay.names)),
                    ncol = 1 )
            }
        }
        ### continuous variable ####
    } else{
        #### scatter plots for individual assays ####
        all.scat.var.pca.plots <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('PCA plots of for ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose)
                var <- PC <- NULL
                pca.data <- as.data.frame(all.pca.data[[x]]$pca.data[ , seq_len(nb.pcs)])
                colnames(pca.data) <- paste0('PC', seq_len(nb.pcs), ' (',all.pca.data[[x]]$pc.var[seq_len(nb.pcs)], '%)')
                pca.data$var <- colData(se.obj)[, variable]
                printColoredMessage(
                    message = paste0('-Create all possible scatter plots of the first ', nb.pcs, ' PCs.'),
                    color = 'blue',
                    verbose = verbose)
                pca.data <- tidyr::pivot_longer(
                    data = pca.data,
                    -var,
                    names_to = 'PCs',
                    values_to = 'PC')
                pca.data <- as.data.frame(pca.data)
                plot.p <- ggplot(pca.data, aes(x = var, y = PC)) +
                    geom_point(
                        color = stroke.color,
                        pch = 19,
                        stroke = stroke.size,
                        size = 3,
                        alpha = points.alpha) +
                    facet_wrap(~PCs, scales = 'free') +
                    xlab(variable) +
                    ylab('PC') +
                    ggtitle(x) +
                    geom_smooth(formula = y ~ x, method = 'lm') +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', size = 1),
                        axis.title.x = element_text(size = 10),
                        axis.title.y = element_text(size = 10),
                        plot.title = element_text(size = 10),
                        axis.text.x = element_text(size = 8),
                        axis.text.y = element_text(size = 8),
                        strip.text.x = element_text(size = 10),
                        legend.position = 'none')
                plot.p
            })
        names(all.scat.var.pca.plots) <- levels(assay.names)
        if(length(assay.names) > 1){
            if(is.null(plot.ncol)) plot.ncol <- ncol(combn(length(seq_len(nb.pcs)), 2))
            overall.scat.var.pca.plot <- ggarrange(
                plotlist = all.scat.var.pca.plots,
                nrow = plot.nrow,
                ncol = plot.ncol)
            if(isTRUE(plot.output))
                suppressMessages(print(overall.scat.var.pca.plot))
        }
    }

    # save the results ####
    printColoredMessage(
        message = '-- Save the PCA plots:',
        color = 'magenta',
        verbose = verbose)
    ## add the pca plots to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = 'The PCA plots of individual assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)
        if(class(colData(se.obj)[[variable]]) %in% c('character', 'factor')){
            if(plot.type == 'scatter'){
                for (x in levels(assay.names)) {
                    if (fast.pca) {
                        se.obj@metadata[['metric']][[x]][['fastPCA']][[variable]][['pca.scat.plot']] <- all.scat.pca.plots.assay[[x]]
                    } else se.obj@metadata[['metric']][[x]][['PCA']][[variable]][['pca.scat.plot']] <- all.scat.pca.plots.assay[[x]]
                }
            } else{
                for (x in levels(assay.names)) {
                    if (fast.pca) {
                        se.obj@metadata[['metric']][[x]][['fastPCA']][[variable]][['pca.box.plot']] <- all.boxplot.pca.plots[[x]]
                    } else se.obj@metadata[['metric']][[x]][['PCA']][[variable]][['pca.box.plot']] <- all.boxplot.pca.plots[[x]]
                }
            }
        } else{
            for (x in levels(assay.names)) {
                if (fast.pca) {
                    se.obj@metadata[['metric']][[x]][['fastPCA']][[variable]][['pca.scat.var.plot']] <- all.scat.var.pca.plots[[x]]
                } else se.obj@metadata[['metric']][[x]][['PCA']][[variable]][['pca.scat.var.plot']] <- all.scat.var.pca.plots[[x]]
            }
        }

        ### overall pca plots ####
        if(length(assay.names) > 1){
            if (!'plot' %in%  names(se.obj@metadata)) {
                se.obj@metadata[['plot']] <- list()
            }
            if (!'PCA' %in%  names(se.obj@metadata[['plot']])) {
                se.obj@metadata[['plot']][['PCA']] <- list()
            }
            if(fast.pca){
                if(class(colData(se.obj)[[variable]]) %in% c('character', 'factor')){
                    if (!'fastPCA' %in%  names(se.obj@metadata[['plot']][['PCA']])) {
                        se.obj@metadata[['plot']][['PCA']][['fastPCA']] <- list()
                    }
                    if (!variable %in%  names(se.obj@metadata[['plot']][['PCA']][['fastPCA']])) {
                        se.obj@metadata[['plot']][['PCA']][['fastPCA']][[variable]] <- list()
                    }
                    if(plot.type == 'scatter'){
                        if (!'ScatPCA' %in%  names(se.obj@metadata[['plot']][['PCA']][['fastPCA']])) {
                            se.obj@metadata[['plot']][['PCA']][['fastPCA']][[variable]][['ScatPCA']] <- list()
                            se.obj@metadata[['plot']][['PCA']][['fastPCA']][[variable]][['ScatPCA']] <- overall.scat.pca.plot
                        }
                    } else{
                        if (!'BoxPCA' %in%  names(se.obj@metadata[['plot']][['PCA']][['fastPCA']])) {
                            se.obj@metadata[['plot']][['PCA']][['fastPCA']][[variable]][['BoxPCA']] <- list()
                            se.obj@metadata[['plot']][['PCA']][['fastPCA']][[variable]][['BoxPCA']] <- overall.boxplot.pca.plot
                        }
                    }

                } else{
                    if (!'fastPCA' %in%  names(se.obj@metadata[['plot']][['PCA']])) {
                        se.obj@metadata[['plot']][['PCA']][['fastPCA']] <- list()
                    }
                    if (!variable %in%  names(se.obj@metadata[['plot']][['PCA']][['fastPCA']])) {
                        se.obj@metadata[['plot']][['PCA']][['fastPCA']][[variable]] <- list()
                    }
                    if (!'ScatVarPCA' %in%  names(se.obj@metadata[['plot']][['PCA']][['fastPCA']])){
                        se.obj@metadata[['plot']][['PCA']][['fastPCA']][[variable]][['ScatVarPCA']] <- list()
                    }
                    se.obj@metadata[['plot']][['PCA']][['fastPCA']][[variable]][['ScatVarPCA']] <- overall.scat.var.pca.plot
                }
            }else{
                if(class(colData(se.obj)[[variable]]) %in% c('character', 'factor')){
                    if (!'PCA' %in%  names(se.obj@metadata[['plot']][['PCA']])) {
                        se.obj@metadata[['plot']][['PCA']][['PCA']] <- list()
                    }
                    if (!variable %in%  names(se.obj@metadata[['plot']][['PCA']][['PCA']])) {
                        se.obj@metadata[['plot']][['PCA']][['PCA']][[variable]] <- list()
                    }
                    if(plot.type == 'scatter'){
                        if (!'ScatPCA' %in%  names(se.obj@metadata[['plot']][['PCA']][['PCA']])) {
                            se.obj@metadata[['plot']][['PCA']][['PCA']][[variable]][['ScatPCA']] <- list()
                            se.obj@metadata[['plot']][['PCA']][['PCA']][[variable]][['ScatPCA']] <- overall.scat.pca.plot
                        }
                    } else{
                        if (!'BoxPCA' %in%  names(se.obj@metadata[['plot']][['PCA']][['PCA']])) {
                            se.obj@metadata[['plot']][['PCA']][['PCA']][[variable]][['BoxPCA']] <- list()
                            se.obj@metadata[['plot']][['PCA']][['PCA']][[variable]][['BoxPCA']] <- overall.boxplot.pca.plot
                        }
                    }

                } else{
                    if (!'PCA' %in%  names(se.obj@metadata[['plot']][['PCA']])) {
                        se.obj@metadata[['plot']][['PCA']][['PCA']] <- list()
                    }
                    if (!variable %in%  names(se.obj@metadata[['plot']][['PCA']][['PCA']])) {
                        se.obj@metadata[['plot']][['PCA']][['PCA']][[variable]] <- list()
                    }
                    if (!'ScatVarPCA' %in%  names(se.obj@metadata[['plot']][['PCA']][['PCA']])) {
                        se.obj@metadata[['plot']][['PCA']][['PCA']][[variable]][['ScatVarPCA']] <- list()
                        se.obj@metadata[['plot']][['PCA']][['PCA']][[variable]][['ScatVarPCA']] <- overall.scat.var.pca.plot
                    }
                }
            }
        }
        printColoredMessage(message = '------------The plotPCA function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else if (save.se.obj == FALSE) {
        ## return a list ####
        printColoredMessage(
            message = 'The PCA plots are outputed as a list.',
            color = 'blue',
            verbose = verbose)
    }
    printColoredMessage(message = '------------The plotPCA function finished.',
                        color = 'white',
                        verbose = verbose)
    if(length(assay.names) == 1){
        if(plot.type == 'scatter'){
            return(pca.plots = list(all.scat.pca.plots.assay = all.scat.pca.plots.assay))
        } else return(pca.plots = list(all.boxplot.pca.plots = all.boxplot.pca.plots))
    } else{
        if(plot.type == 'scatter'){
            return(pca.plots = list(
                all.scat.pca.plots.assay = all.scat.pca.plots.assay,
                overall.scat.pca.plot = overall.scat.pca.plot))
        } else return(pca.plots = list(
            all.boxplot.pca.plots = all.boxplot.pca.plots,
            overall.boxplot.pca.plot = overall.boxplot.pca.plot))
    }
}
