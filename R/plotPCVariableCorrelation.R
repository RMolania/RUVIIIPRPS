#' Plot the vector correlation.

#' @author Ramyar Molania

#' @description
#' This function generate a dot-line plot between the first cumulative PCs and the correlation coefficient (R2) obtained
#' from the vector correlation analysis. An ideal normalization should results a low correlation with unwanted variation
#' variables and high correlation with known biology.

#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to generate a vector correlation. The default is "all, which indicates all
#' the assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. Indicates a name of the column in the sample annotation of the SummarizedExperiment object.
#' The variable must be a categorical variable.
#' @param fast.pca Logical. Indicates whether to use the fast PCA or PCA results computed by the computePCA function. The
#' default is 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to use to plot the vector correlation. The default is 10.
#' @param plot.output Logical. If TRUE, the individual vector correlation plot(s) will be printed while functions is running.
#' @param save.se.obj Logical. Indicates whether to save the vector correlation plots to the meta data of the
#' SummarizedExperiment object or to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list that contains all the vector correlation plots for the individual
#' assay(s).

#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom grDevices colorRampPalette
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
        verbose = TRUE) {
    printColoredMessage(message = '------------The plotPCVariableCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (is.list(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s).')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    } else if(length(variable) > 1){
        stop('The "variable" must contain only one variable.')
    } else if (!variable %in% colnames(se.obj@colData)){
        stop('The "variable" cannot be found in the SummarizedExperiment object.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop('The "variable" must be a categorical varible.')
    } else if (length(unique(se.obj@colData[, variable])) < 2) {
        stop('The "variable" must have at least two levels.')
    }
    if (sum(is.na(se.obj@colData[[variable]])) > 0){
        stop(paste0('The "', variable, '" contains NA.',
        ' Run the checkSeObj function with "remove.na = both"',
        ', then "computePCA"-->"computePCVariableCorrelation"-->"plotPCVariableCorrelation".'))
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names, levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # select colors ####
    if(length(levels(assay.names)) < 9 ){
        data.sets.colors <- RColorBrewer::brewer.pal(8, 'Dark2')[1:length(levels(assay.names))]
        names(data.sets.colors) <- levels(assay.names)
    } else{
        colfunc <- grDevices::colorRampPalette( RColorBrewer::brewer.pal(8, 'Dark2'))
        data.sets.colors <- colfunc(n = length(levels(assay.names)))
        names(data.sets.colors) <- levels(assay.names)
    }

    # check metric exist ####
    m.out <- lapply(
        levels(assay.names),
        function(x) {
            if (!x %in% names(se.obj@metadata[['metric']]))
                stop(paste0('Any metrics has not been computed yet on the  ', x, ' data.'))
        })

    # obtain vector correlations ####
    printColoredMessage(
        message = '-- Obtain the computed vector correlations from the SummarizedExperiment object:',
        color = 'magenta',
        verbose = verbose
    )
    all.pcs.vect.corr <- lapply(
        levels(assay.names),
        function(x) {
            if (!'pcs.vect.corr' %in% names(se.obj@metadata[['metric']][[x]]) ) {
                stop(paste0('Any "pcs.vect.corr" is found for the ', x, ' data. Please run the "computePCVariableCorrelation" function first.'))
            }
            if (!variable %in% names(se.obj@metadata[['metric']][[x]][['pcs.vect.corr']]) ) {
                stop(paste0('The "pcs.vect.corr" is not found for the ', variable, ' variable and the ', x, ' data.'))
            }
            if (!'corrs' %in% names(se.obj@metadata[['metric']][[x]][['pcs.vect.corr']][[variable]]) ) {
                stop(paste0('The "pcs.vect.corr" is not found for the ', variable, ' variable and the ', x, ' data.'))
            }
            printColoredMessage(
                message = paste0('-Obtain the vector correlations for', x , ' data.'),
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
            data.to.plot <- data.frame(
                vec.corr = all.pcs.vect.corr[[x]],
                pcs = seq_len(nb.pcs)
                )
            vect.corr.plot <- ggplot(data.to.plot, aes(x = pcs, y = vec.corr, group = 1)) +
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
        all.pcs.vect.corr <- as.data.frame(do.call(cbind, all.pcs.vect.corr))
        all.pcs.vect.corr <- all.pcs.vect.corr %>%
            mutate(pcs = c(1:nb.pcs)) %>%
            tidyr::pivot_longer(
                -pcs,
                names_to = 'datasets',
                values_to = 'vec.corr') %>%
            dplyr::mutate(datasets = factor(datasets, levels = levels(assay.names)))
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
            message = '-The vector correlation plots are outputed as a list.',
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

