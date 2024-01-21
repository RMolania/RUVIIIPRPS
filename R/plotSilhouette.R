#' is used to compute the adjusted rand index (ARI).

#' @author Ramyar Molania

#' @description
#' This functions computes the adjusted rand index for given a categorical variable using the first PCs of the assay(s)
#' in a SummarizedExperiment object.

#' @details
#' The ARI64 is the corrected-for-chance version of the Rand index. The ARI measures the percentage of matches between
#' two label lists. We used the ARI to assess the performance of normalization methods in terms of sample subtype
#' separation and batch mixing. We first calculated PCs and used the first three PCs to perform ARI.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to compute PCA. By default all the assays of the SummarizedExperiment object will be selected.
#' @param variables Symbol. Indicates the column name in the SummarizedExperiment object that contains a categorical
#' variable such as sample types or batches.
#' @param plot.type TTT
#' @param silhouette.method TTT
#' @param plot.output Logical. Indicates whether to plot the ARI, by default it is set to FALSE.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object, by default it is
#' set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the
#' execution of the functions, by default it is set to TRUE.

#' @return A SummarizedExperiment object or a list that containing the computed ARI on the categorical variable.

#' @importFrom SummarizedExperiment assays assay
#' @import ggplot2
#' @export

plotSilhouette <- function(
        se.obj,
        assay.names = 'all',
        variables,
        plot.type = 'combined.plot',
        silhouette.method = 'sil.euclidian',
        plot.output = TRUE,
        save.se.obj = TRUE,
        assess.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The plotSilhouette function starts:',
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
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- as.factor(unlist(assay.names))
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variables,
            remove.na = 'none',
            verbose = verbose)}

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
                    geom_text_repel(aes(label = datasets),
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
                geom_text_repel(aes(label = datasets),
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
