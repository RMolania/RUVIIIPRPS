#' plot the adjusted rand index (ARI).

#' @author Ramyar Molania

#' @description
#' This functions generates barplots of adjusted rand index for individual assays. If two variables are provided, the
#' function create scatter plots of the adjusted rand index of each variable for individual assays.

#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to generate barplot or scatter plots of the computed adjusted rand index. By default all
#' the assays of the SummarizedExperiment object will be selected.
#' @param variables Symbol. Indicates one or two column names in the SummarizedExperiment object that contains categorical
#' variables such as sample subtypes or batches.
#' @param ari.method Symbol.Indicates what computed ARI methods should be used for plotting. The "ari.method" must be
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

#' @return A SummarizedExperiment object or a list that containing all the plots of the computed ARI on the categorical
#' variable.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export

plotARI <- function(
        se.obj,
        assay.names = 'all',
        variables,
        ari.method = 'hclust.complete.euclidian',
        plot.type = 'single.plot',
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The plotARI function starts:',
                        color = 'white',
                        verbose = verbose)
    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty')
    } else if (is.null(variables)) {
        stop('The "variables" cannot be empty')
    } else if (!plot.type %in% c('single.plot', 'combined.plot')) {
        stop('The "plot.type" must be one of the "single.plot" or "combined.plot".')
    } else if (plot.type == 'combined.plot') {
        if (length(variables) == 1)
            stop('To plot combined ARI, two variables must be provided.')
    } else if (plot.type == 'single.plot') {
        if (length(variables) > 1)
            stop('To plot "single.plot" ARI, only one variable must be provided.')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- as.factor(unlist(assay.names))
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # check ari metric exist ####
    m.out <- lapply(
        levels(assay.names),
        function(x) {
            if (!'ari' %in% names(se.obj@metadata[['metric']][[x]]))
                stop(paste0('Any ARI analysis has not been computed yet on the  ', x, ' assay'))
        })

    # plots ####
    ## single plot ####
    if (plot.type == 'single.plot') {
        # obtain ari ####
        printColoredMessage(
            message = paste0('--Obtain computed ARI for the from the SummarizedExperiment object:'),
            color = 'magenta',
            verbose = verbose
        )
        all.ari <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('-Obtain the ARI for the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose
                )
                if (!ari.method %in% names(se.obj@metadata[['metric']][[x]][['ari']])) {
                    stop(paste0('Any ', ari.method , ' has not been computed yet for the ', x, ' assay.'))
                }
                if (!variables %in% names(se.obj@metadata[['metric']][[x]][['ari']][[ari.method]])) {
                    stop(paste0('The ', ari.method , 'has not been computed yet for the ', variables,' variable and the ', x, ' assay.'))
                }
                se.obj@metadata[['metric']][[x]][['ari']][[ari.method]][[variables]]$ari
            })
        names(all.ari) <- levels(assay.names)
        ## individual plots ####
        printColoredMessage(
            message = paste0('--Plot ARI:'),
            color = 'magenta',
            verbose = verbose
        )
        all.single.ari.plots <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('-Plot ARI for the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose
                )
                ari.plot <- ggplot() +
                    geom_col(aes(y = all.ari[[x]], x = 1)) +
                    ylab('Adjusted rand index ') +
                    xlab(x) +
                    ggtitle(paste0('Adjusted rand index, ', variables)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        plot.title = element_text(size = 16),
                        axis.text.x = element_text(size = 0),
                        axis.text.y = element_text(size = 12))
                if(plot.output)
                    print(ari.plot)
            })
        names(all.single.ari.plots) <- levels(assay.names)
        ## overall plots ####
        everything <- datasets <- ari <- NULL
        if (length(assay.names) > 1) {
            printColoredMessage(
                message = paste0('-- Generate a ARI plot for all the assays(s):'),
                color = 'magenta',
                verbose = verbose
            )
            overall.single.ari.plot <- as.data.frame(all.ari) %>%
                tidyr::pivot_longer(
                    everything(),
                    names_to = 'datasets',
                    values_to = 'ari')
            overall.single.ari.plot <- ggplot(overall.single.ari.plot, aes(x = datasets, y = ari)) +
                geom_col() +
                ylab('Adjusted rand index ') +
                xlab('Datasets') +
                ggtitle(paste0('Adjusted rand index, ', variables)) +
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
                print(overall.single.ari.plot)
        }
    } else if (plot.type == 'combined.plot') {
        # combine plots ####
        printColoredMessage(
            message = paste0('--Obtain computed ARI for the from the SummarizedExperiment object:'),
            color = 'magenta',
            verbose = verbose
        )
        all.ari <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('-Obtain ARI for the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose
                )
                if (!ari.method %in% names(se.obj@metadata[['metric']][[x]][['ari']])) {
                    stop(paste0('The ', ari.method ,'has not been computed yet for the ', variables, ' variable and the ', x, ' assay.'))
                }
                for (i in variables) {
                    if (!i %in% names(se.obj@metadata[['metric']][[x]][['ari']][[ari.method]])) {
                        stop(paste0('The ', ari.method ,' has not been computed yet for the ', i, ' variable and the ', x, ' assay.'))
                    }
                }
                ari <- c()
                for (i in 1:length(variables))
                    ari[i] <-
                    se.obj@metadata[['metric']][[x]][['ari']][[ari.method]][[variables[i]]]$ari
                return(ari)
            })
        names(all.ari) <- levels(assay.names)
        ## individual plots ####
        printColoredMessage(
            message = paste0('-- Plot combined ARI'),
            color = 'magenta',
            verbose = verbose
        )
        datasets <- NULL
        all.combined.ari.plots <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('-Plot ARI for the ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose
                )
                all.aris <- as.data.frame(t(as.data.frame(all.ari[[x]])))
                row.names(all.aris) <- x
                colnames(all.aris) <- variables
                all.aris$datasets <- row.names(all.aris)
                ggplot(all.aris, aes_string(x = variables[1], y = variables[2])) +
                    geom_point() +
                    ggrepel::geom_text_repel(aes(label = datasets),
                                             hjust = 0,
                                             vjust = 0) +
                    xlab(paste0('Adjusted rand index (', variables[1], ')')) +
                    ylab(paste0('Adjusted rand index (', variables[2], ')')) +
                    ggtitle(paste0('Adjusted rand index')) +
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
        names(all.combined.ari.plots) <- levels(assay.names)

        ## overall plots ####
        datasets <- NULL
        if (length(assay.names) > 1) {
            printColoredMessage(
                message = paste0('-- Generate the combined ARI plot for all the assays(s):'),
                color = 'magenta',
                verbose = verbose
            )
            all.ari <- as.data.frame(t(as.data.frame(all.ari)))
            colnames(all.ari) <- variables
            all.ari$datasets <- row.names(all.ari)
            overall.combined.ari.plot <-
                ggplot(all.ari, aes_string(x = variables[1], y = variables[2])) +
                geom_point() +
                ggrepel::geom_text_repel(aes(label = datasets),
                                         hjust = 0,
                                         vjust = 0) +
                xlab(paste0('Adjusted rand index (', variables[1], ')')) +
                ylab(paste0('Adjusted rand index (', variables[2], ')')) +
                ggtitle(paste0('Adjusted rand index')) +
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
        printColoredMessage(message = '-- Save the ARI plots to the metadata of the SummarizedExperiment object:',
                            color = 'magenta',
                            verbose = verbose)
        if (plot.type == 'single.plot') {
            for (i in variables) {
                for (x in levels(assay.names)) {
                    se.obj@metadata[['metric']][[x]][['ari']][[ari.method]][[i]][['ari.single.plot']] <- all.single.ari.plots[[x]]
                }
            }
        } else if (plot.type == 'combined.plot') {
            for (i in variables) {
                for (x in levels(assay.names)) {
                    se.obj@metadata[['metric']][[x]][['ari']][[ari.method]][[i]][['ari.combined.plot']] <- all.combined.ari.plots[[x]]
                }
            }
        }
        printColoredMessage(
            '-The ARI plots of induvial assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose
        )
        if (length(assay.names) > 1) {
            variables <- paste0(variables, collapse = '&')
            if (!'plot' %in%  names(se.obj@metadata)) {
                se.obj@metadata[['plot']] <- list()
            }
            if (!'ARI' %in%  names(se.obj@metadata[['plot']])) {
                se.obj@metadata[['plot']][['ARI']] <- list()
            }
            if (!ari.method %in%  names(se.obj@metadata[['plot']][['ARI']])) {
                se.obj@metadata[['plot']][['ARI']][[ari.method]] <- list()
            }
            if (plot.type == 'single.plot') {
                if (!'ari.single.plot' %in%  names(se.obj@metadata[['plot']][['ARI']][[ari.method]])) {
                    se.obj@metadata[['plot']][['ARI']][[ari.method]][['ari.single.plot']] <- list()
                }
                if (!variables %in%  names(se.obj@metadata[['plot']][['ARI']][['ari.single.plot']])) {
                    se.obj@metadata[['plot']][['ARI']][[ari.method]][['ari.single.plots']][[variables]] <- list()
                }
                se.obj@metadata[['plot']][['ARI']][[ari.method]][['ari.single.plot']][[variables]] <- overall.single.ari.plot
            } else if (plot.type == 'combined.plot') {
                if (!'ari.combined.plot' %in%  names(se.obj@metadata[['plot']][['ARI']][[ari.method]])) {
                    se.obj@metadata[['plot']][['ARI']][[ari.method]][['ari.combined.plot']] <- list()
                }
                if (!variables %in%  names(se.obj@metadata[['plot']][['ARI']][[ari.method]][['ari.combined.plot']])) {
                    se.obj@metadata[['plot']][['ARI']][[ari.method]][['ari.combined.plot']][[variables]] <- list()
                }
                se.obj@metadata[['plot']][['ARI']][[ari.method]][['ari.combined.plot']][[variables]] <- overall.combined.ari.plot
            }
            printColoredMessage(
                message = paste0('-The ARI plot all assays are saved to metadata@plot'),
                color = 'blue',
                verbose = verbose)
        }
        printColoredMessage(message = '------------The plotARI function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)

        ## return only the adjusted rand index results ####
    } else if (save.se.obj == FALSE) {
        printColoredMessage(
            message = paste0('-All the ARI plots re saved as list.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The plotARI function finished.',
                            color = 'white',
                            verbose = verbose)
        if (plot.type == 'single.plot') {
            if(length(assay.names) == 1){
                if(length(assay.names) > 1){
                    return( all.ari.plots = list(all.single.ari.plots = all.single.ari.plots))
                }else{
                    return( all.ari.plots = list(
                        all.single.ari.plots = all.single.ari.plots,
                        overall.single.ari.plot = overall.single.ari.plot))}
            }
        } else if (plot.type == 'combined.plot') {
            if(length(assay.names) == 1){
                return(all.ari.plots = list(all.combined.ari.plots = all.combined.ari.plots))
            }else{
                return(all.ari.plots = list(
                    all.combined.ari.plots = all.combined.ari.plots,
                    overall.combined.ari.plot = overall.combined.ari.plot))
            }
        }
    }
}
