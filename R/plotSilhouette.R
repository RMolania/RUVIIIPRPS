#' is used to compute the average Silhouette coefficients.

#' @author Ramyar Molania

#' @description
#' This functions generates barplots of average Silhouette coefficients for individual assays. If two variables are provided, the
#' function creates scatter plots of the average Silhouette coefficients of each variable for the individual assays.


#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to generate barplot or scatter plots of the computed adjusted rand index. By default all
#' the assays of the SummarizedExperiment object will be selected.
#' @param variables Symbol. Indicates one or two column names in the SummarizedExperiment object that contains categorical
#' variables such as sample subtypes or batches.
#' @param silhouette.method Symbol. Indicates what computed ARI methods should be used for plotting. The "ari.method" must be
#' specified based on the "computeARI" function. The default is "hclust.complete.euclidian", which is the default of the
#' the "computeARI" function. We refer to the "computeARI" function for more detail
#' @param plot.type Symbol.Indicates how to plot the adjusted rand index. The options are "single.plot" and "combined.plot".
#' If a variable is provided, then the "plot.type" must be set to "single.plot", so, the function generates a barplot of the
#' adjusted rand index. If two variables are provided and "plot.type" is set to "combined.plot", then the function generates
#' a scatter plot of the adjusted rand index of each variable against each other.
#' @param plot.output Logical. If TRUE, the individual barplots or scatter plots will be printed while functions is running.
#' @param save.se.obj Logical. Indicates whether to save the plots in the metadata of the SummarizedExperiment  object
#' or to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list that containing all the plots of the computed average Silhouette coefficients
#' on the categorical variable.

#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @importFrom SummarizedExperiment assays assay
#' @importFrom tidyr pivot_longer
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#' @export

plotSilhouette <- function(
        se.obj,
        assay.names = 'all',
        variables,
        silhouette.method = 'sil.euclidian',
        plot.type = 'combined.plot',
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The plotSilhouette function starts:',
                        color = 'white',
                        verbose = verbose)
    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty')
    } else if (!is.vector(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s) or "assay.names == all".')
    } else if (is.null(variables)) {
        stop('The "variables" cannot be empty')
    } else if (!plot.type %in% c('single.plot', 'combined.plot')) {
        stop('The "plot.type" must be one of the "single.plot" or "combined.plot".')
    } else if (plot.type == 'combined.plot') {
        if (length(variables) == 1)
            stop('To plot combined ARI, two variables must be provided.')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # check ari metric exist ####
    m.out <- lapply(
        levels(assay.names),
        function(x) {
            if (!'silhouette' %in% names(se.obj@metadata[['metric']][[x]]))
                stop(paste0('Any Silhouette analysis has not been computed yet on the  ', x,' assay'))
        })

    # plots ####
    ## single plot ####
    if (plot.type == 'single.plot') {
        # obtain silhouette ####
        printColoredMessage(
            message = paste0('--Obtain computed Silhouette for the from the SummarizedExperiment object:'),
            color = 'magenta',
            verbose = verbose
        )
        all.silhouette <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('-Obtain the Silhouette for the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose
                )
                if (!silhouette.method %in% names(se.obj@metadata[['metric']][[x]][['silhouette']])) {
                    stop(paste0('Any ', silhouette.method ,' has not been computed yet for the ', x, ' assay.'))
                }
                if (!variables %in% names(se.obj@metadata[['metric']][[x]][['silhouette']][[silhouette.method]])) {
                    stop(paste0('The ', silhouette.method ,'has not been computed yet for the ', variables, ' variable and the ', x,' assay.'))
                }
                silhouette <- se.obj@metadata[['metric']][[x]][['silhouette']][[silhouette.method]][[variables]]$silhouette
            })
        names(all.silhouette) <- levels(assay.names)
        ## individual plots ####
        all.single.silhouette.plots <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('-Plot Silhouette for the ', x, 'data.'),
                    color = 'blue',
                    verbose = verbose
                )
                p.silhouette <- ggplot() +
                    geom_col(aes(y = all.silhouette[[x]], x = 1)) +
                    ylab('Silhouette ') +
                    xlab(x) +
                    ggtitle(paste0('Silhouette, ', variables)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        plot.title = element_text(size = 16),
                        axis.text.x = element_text(size = 0),
                        axis.text.y = element_text(size = 12)
                    )
                if(plot.output)
                    print(p.silhouette)
            })
        names(all.single.silhouette.plots) <- levels(assay.names)
        everything <- datasets <- silhou.coff <- NULL

        ## overall plots ####
        if (length(assay.names) > 1) {
            printColoredMessage(
                message = paste0('Generate a Silhouette plot for all the assays(s)'),
                color = 'magenta',
                verbose = verbose
            )
            overall.single.silhouette.plot <- as.data.frame(all.silhouette) %>%
                tidyr::pivot_longer(
                    everything(),
                    names_to = 'datasets',
                    values_to = 'silhou.coff')
            overall.single.silhouette.plot <- ggplot(overall.single.silhouette.plot,
                       aes(x = datasets, y = silhou.coff)) +
                geom_col() +
                ylab('Silhouette ') +
                xlab('Datasets') +
                ggtitle(paste0('Silhouette, ', variables)) +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    plot.title = element_text(size = 15),
                    axis.text.x = element_text(
                        size = 12,
                        angle = 25,
                        hjust = 1),
                    axis.text.y = element_text(size = 12))
            if(plot.output)
                print(overall.single.silhouette.plot)
        }

    } else if (plot.type == 'combined.plot') {
        # combine plots ####
        printColoredMessage(
            message = paste0('--Obtain computed silhouette for the from the SummarizedExperiment object:'),
            color = 'magenta',
            verbose = verbose
        )
        all.silhouette <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('-Obtain silhouette for the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose
                )
                if (!silhouette.method %in% names(se.obj@metadata[['metric']][[x]][['silhouette']])) {
                    stop(paste0('The ', silhouette.method ,'has not been computed yet for the ', variables,' variable and the ', x,' data.'))
                }
                for (i in variables) {
                    if (!i %in% names(se.obj@metadata[['metric']][[x]][['silhouette']][[silhouette.method]])) {
                        stop(paste0('The ', silhouette.method ,' has not been computed yet for the ', i, ' variable and the ', x,' data.'))
                    }
                }
                silhouette <- c()
                for (i in 1:length(variables))
                    silhouette[i] <-
                    se.obj@metadata[['metric']][[x]][['silhouette']][[silhouette.method]][[variables[i]]]$silhouette
                return(silhouette)
            })
        names(all.silhouette) <- levels(assay.names)
        ## individual plots ####
        datasets <- NULL
        all.combined.silhouette.plots <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('-Plot Silhouette for the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose
                )
                all.silhouettes <- as.data.frame(t(as.data.frame(all.silhouette[[x]])))
                row.names(all.silhouettes) <- x
                colnames(all.silhouettes) <- variables
                all.silhouettes$datasets <- row.names(all.silhouettes)
                ggplot(all.silhouettes, aes_string(x = variables[1], y = variables[2])) +
                    geom_point() +
                    ggrepel::geom_text_repel(aes(label = datasets),
                                    hjust = 0,
                                    vjust = 0) +
                    xlab(paste0('Silhouette (', variables[1], ')')) +
                    ylab(paste0('Silhouette (', variables[2], ')')) +
                    ggtitle(paste0('Silhouette')) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 18),
                        axis.title.y = element_text(size = 18),
                        plot.title = element_text(size = 15),
                        axis.text.x = element_text(size = 12),
                        axis.text.y = element_text(size = 12)
                    )
            })
        names(all.combined.silhouette.plots) <- levels(assay.names)

        ## overall plots ####
        if (length(assay.names) > 1) {
            printColoredMessage(
                message = paste0('Plot Silhouette for all the assays(s)'),
                color = 'blue',
                verbose = verbose
            )
            all.silhouette <- as.data.frame(t(as.data.frame(all.silhouette)))
            colnames(all.silhouette) <- variables
            all.silhouette$datasets <- row.names(all.silhouette)
            overall.combined.silhouette.plot <- ggplot(all.silhouette, aes_string(x = variables[1], y = variables[2])) +
                geom_point() +
                ggrepel::geom_text_repel(aes(label = datasets),
                                hjust = 0,
                                vjust = 0) +
                xlab(paste0('Silhouette (', variables[1], ')')) +
                ylab(paste0('Silhouette (', variables[2], ')')) +
                ggtitle(paste0('Silhouette')) +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    plot.title = element_text(size = 15),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 12)
                )
        }
    }
    # save the results ####
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(message = '-- Save the Silhouette plots to the metadata of the SummarizedExperiment object:',
                            color = 'magenta',
                            verbose = verbose)
        if (plot.type == 'single.plots') {
            for (i in variables) {
                for (x in levels(assay.names)) {
                    se.obj@metadata[['metric']][[x]][['silhouette']][[silhouette.method]][[i]][['silhouette.single.plot']] <-
                        all.single.silhouette.plots[[x]]
                }
            }
        } else if (plot.type == 'combined.plot') {
            for (i in variables) {
                for (x in levels(assay.names)) {
                    se.obj@metadata[['metric']][[x]][['silhouette']][[silhouette.method]][[i]][['silhouette.combined.plot']] <-
                        all.combined.silhouette.plots[[x]]
                }
            }
        }
        printColoredMessage(
            'The Silhouette plots of induvial assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose
        )
        if (length(assay.names) > 1) {
            variables <- paste0(variables, collapse = '&')
            if (!'plot' %in%  names(se.obj@metadata)) {
                se.obj@metadata[['plot']] <- list()
            }
            if (!'Silhouette' %in%  names(se.obj@metadata[['plot']])) {
                se.obj@metadata[['plot']][['Silhouette']] <- list()
            }
            if (!silhouette.method %in%  names(se.obj@metadata[['plot']][['Silhouette']])) {
                se.obj@metadata[['plot']][['Silhouette']][[silhouette.method]] <- list()
            }
            if (plot.type == 'single.plots') {
                if (!'silhouette.single.plot' %in%  names(se.obj@metadata[['plot']][['Silhouette']][[silhouette.method]])) {
                    se.obj@metadata[['plot']][['Silhouette']][[silhouette.method]][['silhouette.single.plot']] <- list()
                }
                if (!variables %in%  names(se.obj@metadata[['plot']][['Silhouette']][[silhouette.method]][['silhouette.single.plot']])) {
                    se.obj@metadata[['plot']][['Silhouette']][[silhouette.method]][['silhouette.single.plot']][[variables]] <- list()
                }
                se.obj@metadata[['plot']][['Silhouette']][[silhouette.method]][['silhouette.single.plots']][[variables]] <- overall.single.silhouette.plot
            } else if (plot.type == 'combined.plot') {
                if (!'silhouette.combined.plot' %in%  names(se.obj@metadata[['plot']][['Silhouette']][[silhouette.method]])) {
                    se.obj@metadata[['plot']][['Silhouette']][[silhouette.method]][['silhouette.combined.plot']] <- list()
                }
                if (!variables %in%  names(se.obj@metadata[['plot']][['Silhouette']][[silhouette.method]][['silhouette.combined.plot']])) {
                    se.obj@metadata[['plot']][['Silhouette']][[silhouette.method]][['silhouette.combined.plot']][[variables]] <- list()
                }
                se.obj@metadata[['plot']][['Silhouette']][[silhouette.method]][['silhouette.combined.plot']][[variables]] <- overall.combined.silhouette.plot
            }
            printColoredMessage(
                message = paste0(
                    'The Silhouette plot all assays are saved to metadata@plot'),
                color = 'blue',
                verbose = verbose)
        }
        printColoredMessage(message = '------------The plotSilhouette function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)

        ## return only the adjusted rand index results ####
    } else if (save.se.obj == FALSE) {
        printColoredMessage(
            message = paste0('-All the ARI plots re saved as list.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The plotSilhouette function finished.',
                            color = 'white',
                            verbose = verbose)
        if (plot.type == 'single.plot') {
            if(length(assay.names) == 1){
                return(all.silhouette.plots = list(all.single.silhouette.plots = all.single.silhouette.plots))
            }else{
                return(all.silhouette.plots = list(
                        all.single.silhouette.plots = all.single.silhouette.plots,
                        overall.single.silhouette.plot = overall.single.silhouette.plot ))
            }
        } else if (plot.type == 'combined.plot') {
            if(length(assay.names) == 1){
                return(all.silhouette.plots = list(all.combined.silhouette.plots = all.combined.silhouette.plots))
            } else{
                return(all.silhouette.plots = list(
                        all.combined.silhouette.plots = all.combined.silhouette.plots,
                        overall.combined.silhouette.plot = overall.combined.silhouette.plot))
            }
        }
    }
}
