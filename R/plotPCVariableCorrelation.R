#' is used to plot the vector correlation.

#' @author Ramyar Molania

#' @description
#' This function calculates the the vector correlation between the first cumulative PCs of the gene expression (assay)
#' of a SummarizedExperiment object and a categorical variable (i.e. batch).

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the correlation. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param fast.pca Logical. Indicates whether to use the PCA calculated using a specific number of PCs instead of the
#' full range to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of few first cumulative PCs, by default is set to 10.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param plot.output Logical. Indicates whether to plot the correlation statistics, by default it is set to TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.

#' @return SummarizedExperiment A SummarizedExperiment object containing the computed correlation for
#' the continuous variable and if requested the associated plot.

#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export

plotPCVariableCorrelation <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        nb.pcs = 10,
        plot.output = TRUE,
        save.se.obj = TRUE,
        assess.se.obj = TRUE,
        verbose = TRUE) {
    printColoredMessage(message = '------------The plotPCVariableCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- as.factor(assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = 'none',
            verbose = verbose)
    }

    # check metric exist ####
    m.out <- lapply(
        levels(assay.names),
        function(x) {
            if (!x %in% names(se.obj@metadata[['metric']]))
                stop(paste0('Any PCA analysis has not been computed yet on the  ', x, ' data.'))
        })

    # obtain vector correlations ####
    printColoredMessage(
        message = paste0('-- Obtain the vector correlations from the SummarizedExperiment object:'),
        color = 'magenta',
        verbose = verbose
    )
    all.pcs.vect.corr <- lapply(
        levels(assay.names),
        function(x) {
            if (!'pcs.vect.corr' %in% names(se.obj@metadata[['metric']][[x]]) ) {
                stop(paste0('Any "pcs.vect.corr" is found for the ', x, ' assay. Please run the "PCVariableCorrelation" function first.'))
            }
            if (!variable %in% names(se.obj@metadata[['metric']][[x]][['pcs.vect.corr']]) ) {
                stop(paste0('The "pcs.vect.corr" is not found for the ', variable, ' variable and the ', x, ' assay.'))
            }
            if (!'corrs' %in% names(se.obj@metadata[['metric']][[x]][['pcs.vect.corr']][[variable]]) ) {
                stop(paste0('The "pcs.vect.corr" is not found for the ', variable, ' variable and the ', x, ' assay.'))
            }
            printColoredMessage(
                message = paste0('-Obtain the vector correlations fro', x , ' data.'),
                color = 'blue',
                verbose = verbose
            )
            vectos.corrs <- se.obj@metadata[['metric']][[x]][['pcs.vect.corr']][[variable]][['corrs']]
            if(length(vectos.corrs) < nb.pcs){
                stop(paste0('The number of vector correlations of the assay ',
                            x, ' for the ', variable, ' variable are less than', nb.pcs, '.',
                            'Re-run the "PCVariableCorrelation" function with nb.pcs = ', nb.pcs, '.'))
            }
            vectos.corrs[1:nb.pcs]
        })
    names(all.pcs.vect.corr) <- levels(assay.names)
    # individual plots ####
    printColoredMessage(
        message = paste0('-- Generate the vector correlations plots for each assay:'),
        color = 'magenta',
        verbose = verbose
    )
    all.vect.corr.plots <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0('-Generate the vector correlations plots for the ', x , ' data.'),
                color = 'blue',
                verbose = verbose
            )
            to.plot <- data.frame(
                vec.corr = all.pcs.vect.corr[[x]],
                pcs = seq_len(nb.pcs)
                )
            vect.corr.plot <- ggplot(to.plot, aes(x = pcs, y = vec.corr, group = 1)) +
                geom_line(color = 'gray80', size = 1) +
                geom_point(color = 'gray40', size = 3) +
                xlab('Cumulative PCs') +
                ylab('Vector correlations') +
                ggtitle(paste0('Vector correlation, ', x, ', ', variable )) +
                scale_x_continuous(breaks = seq_len(nb.pcs), labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
                scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    plot.title = element_text(size = 16),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 12, angle = 35, vjust = 1, hjust = 1),
                    axis.text.y = element_text(size = 12))
            if(plot.output) print(vect.corr.plot)
            return(vect.corr.plot)
        })
    names(all.vect.corr.plots) <- levels(assay.names)

    # overall plot ####
    if(length(assay.names) > 1){
        printColoredMessage(
            message = paste0('-- Generate the vector correlations plot for all assays:'),
            color = 'magenta',
            verbose = verbose
        )
        datasets <- pcs <- vec.corr <- NULL
        all.pcs.vect.corr <- as.data.frame(all.pcs.vect.corr)
        all.pcs.vect.corr <- all.pcs.vect.corr %>%
            mutate(pcs = c(1:nb.pcs)) %>%
            tidyr::pivot_longer(
                -pcs,
                names_to = 'datasets',
                values_to = 'vec.corr') %>%
            mutate(datasets = factor(datasets, levels = levels(assay.names)))
        overall.vect.corr.plot <- ggplot(all.pcs.vect.corr, aes(x = pcs, y = vec.corr, group = datasets)) +
            geom_line(aes(color = datasets), size = 1) +
            geom_point(aes(color = datasets), size = 3) +
            xlab('Cumulative PCs') +
            ylab('Vector correlations') +
            ggtitle(paste0('Vector correlation, ', variable)) +
            scale_color_manual(values = c(data.sets.colors), name = 'Datasets') +
            scale_x_continuous(breaks = seq_len(nb.pcs), labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                plot.title = element_text(size = 16),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.x = element_text(size = 12, angle = 35, vjust = 1, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 14))
        if(plot.output) overall.vect.corr.plot
    }

    # save the plots ####
    printColoredMessage(
        message = '-- Save all the vector correlatio plots:',
        color = 'magenta',
        verbose = verbose)
    if (save.se.obj == TRUE) {
        ## add results to the SummarizedExperiment object ####
        printColoredMessage(
            message = '-Save the vector correlation plot to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        ## add vector correlation plots of individual assays ####
        ### check the metadata of the SummarizedExperiment object ####
        ### for each assays ####
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist for this assay and this metric and the variable
            if (!'vect.corr.plot' %in% names(se.obj@metadata[['metric']][[x]][['pcs.vect.corr']][[variable]])) {
                se.obj@metadata[['metric']][[x]][['pcs.vect.corr']][[variable]][['vect.corr.plot']] <- list()
            }
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['pcs.vect.corr']][[variable]][['vect.corr.plot']] <- all.vect.corr.plots[[x]]
        }
        printColoredMessage(
            message = 'The vector correlation plot for individal assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)

        ## add overall vector correlation plot of all assays ####
        if(length(assay.names) > 1){
            if (!'plot' %in%  names(se.obj@metadata)) {
                se.obj@metadata[['plot']] <- list()
            }
            if (!'VecCorr' %in%  names(se.obj@metadata[['plot']])) {
                se.obj@metadata[['plot']][['VecCorr']] <- list()
            }
            if (!variable %in%  names(se.obj@metadata[['plot']][['VecCorr']])) {
                se.obj@metadata[['plot']][['VecCorr']][[variable]] <- list()
            }
            se.obj@metadata[['plot']][['VecCorr']][[variable]] <- overall.vect.corr.plot
            printColoredMessage(
                message = paste0('The vector correlation plot of all the assays are saved to metadata@plot'),
                color = 'blue',
                verbose = verbose
            )
        }
        printColoredMessage(message = '------------The plotPCVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    } else if (save.se.obj == FALSE) {
        printColoredMessage(
            message = '-Save the vector correlation plots as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The plotPCVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) == 1){
            return(all.vec.corr = list(all.vect.corr.plots = all.vect.corr.plots))
        } else{
            return(all.vec.corr = list(
                all.vect.corr.plots = all.vect.corr.plots,
                overall.vect.corr.plot = overall.vect.corr.plot))
        }
    }
}

