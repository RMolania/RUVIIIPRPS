#' is used to plot the regression analysis.

#' @author Ramyar Molania

#' @description
#' This function generate a dot-line plot between the first cumulative PCs and the regression R squared obtained
#' from the regression analysis. An ideal normalization should results a low R squared with unwanted variation
#' variables and high R squared with known biology.


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
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object containing the computed correlation for
#' the continuous variable and if requested the associated plot.

#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export

plotPCVariableRegression <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        nb.pcs = 10,
        save.se.obj = TRUE,
        plot.output = TRUE,
        verbose = TRUE) {
    printColoredMessage(message = '------------The plotPCVariableRegression function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty')
    } else if (!is.vector(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s).')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
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

    # check metric ####
    m.out <- lapply(
        levels(assay.names),
        function(x) {
            if (!x %in% names(se.obj@metadata[['metric']]))
                stop(paste0('Any PCA analysis has not been computed yet on the  ', x, 'assay.'))
        })

    # obtain regression r squared ####
    printColoredMessage(
        message = '-- Obtain the regression r squared from the SummarizedExperiment object:',
        color = 'magenta',
        verbose = verbose)
    all.reg.rseq <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('-Obtain the regression r squared for the ', x,  ' data.'),
                color = 'blue',
                verbose = verbose)
            if (!'pcs.lm' %in% names(se.obj@metadata[['metric']][[x]])) {
                stop(paste0('Any  "pcs.lm" has not been computed yet for the ', x, ' assay.'))
            }
            if (!variable %in% names(se.obj@metadata[['metric']][[x]][['pcs.lm']])) {
                stop(paste0('pcs.lm has not been computed yet for the ', variable, ' variable and the ', x, ' assay.'))
            }
            reg.rseq <- se.obj@metadata[['metric']][[x]][['pcs.lm']][[variable]][['rseq']]
            if(length(reg.rseq) < nb.pcs){
                stop(paste0('The number of vector correlations of the assay ',
                            x, ' for the ', variable, ' variable are less than', nb.pcs, '.',
                            'Re-run the "PCVariableCorrelation" function with nb.pcs = ', nb.pcs, '.'))
            }
            reg.rseq[1:nb.pcs]
        })
    names(all.reg.rseq) <- levels(assay.names)

    # individual plots ####
    printColoredMessage(
        message = '-- Generate the pc variable correlations plots for each assay:',
        color = 'magenta',
        verbose = verbose)
    rseq <- NULL
    all.reg.plots <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0('-Generate the pc variable correlations for', x, 'data.'),
                color = 'blue',
                verbose = verbose)
            to.plot <- data.frame(
                rseq = all.reg.rseq[[x]],
                pcs = seq_len(nb.pcs)
            )
            reg.plot <- ggplot(to.plot, aes(x = pcs, y = rseq, group = 1)) +
                geom_line(color = 'gray80', size = 1) +
                geom_point(color = 'gray40', size = 3) +
                xlab('Cumulative PCs') +
                ylab(expression('R'[2])) +
                ggtitle(paste0('Linear regression,', x, ', ', variable )) +
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
            return(reg.plot)
        })
    names(all.reg.plots) <- levels(assay.names)

    # overall plots ####
    if(length(assay.names) > 1){
        printColoredMessage(
            message = '-- -- Generate the pca variable correlations plot for all the assays:',
            color = 'magenta',
            verbose = verbose)
        datasets <- pcs <- vec.corr <- NULL
        all.reg.rseq <- as.data.frame(all.reg.rseq)
        all.reg.rseq <- all.reg.rseq %>%
            mutate(pcs = c(1:nb.pcs)) %>%
            tidyr::pivot_longer(
                -pcs,
                names_to = 'datasets',
                values_to = 'vec.corr') %>%
            mutate(datasets = factor(datasets, levels = levels(assay.names)))
        overall.reg.plot <- ggplot(all.reg.rseq, aes(x = pcs, y = vec.corr, group = datasets)) +
            geom_line(aes(color = datasets), size = 1) +
            geom_point(aes(color = datasets), size = 3) +
            xlab('Cumulative PCs') +
            ylab(expression('R'[2])) +
            ggtitle(paste0('Linear regression, ', variable)) +
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
        if(plot.output)
            print(overall.reg.plot)
    }

    # save the plots ####
    printColoredMessage(
        message = '-- Save all the PCs variable regression plots:',
        color = 'magenta',
        verbose = verbose)
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = '-- The plos are saved to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        ## add vector correlation plots of individual assays ####
        ### check the metadata of the SummarizedExperiment object ####
        ### for each assays ####s
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist for this assay and this metric and the variable
            if (!'vect.corr.plot' %in% names(se.obj@metadata[['metric']][[x]][['pcs.lm']][[variable]])) {
                se.obj@metadata[['metric']][[x]][['pcs.lm']][[variable]][['pcs.lm.plot']] <- list()
            }
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['pcs.lm']][[variable]][['pcs.lm.plot']] <- all.reg.plots[[x]]
        }
        printColoredMessage(
            message = 'The plot for individal assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)

        ## add overall vector correlation plot of all assays ####
        if(length(assay.names) > 1){
            if (!'plot' %in%  names(se.obj@metadata)) {
                se.obj@metadata[['plot']] <- list()
            }
            if (!'PcaReg' %in%  names(se.obj@metadata[['plot']])) {
                se.obj@metadata[['plot']][['PcaReg']] <- list()
            }
            if (!variable %in%  names(se.obj@metadata[['plot']][['PcaReg']])) {
                se.obj@metadata[['plot']][['PcaReg']][[variable]] <- list()
            }
            se.obj@metadata[['plot']][['PcaReg']][[variable]] <- overall.reg.plot
            printColoredMessage(
                message = paste0('The plot of all assays is saved to metadata@plot'),
                color = 'blue',
                verbose = verbose
            )
        }
        printColoredMessage(message = '------------The plotPCVariableRegressions function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
        ## Return only the correlation result
    } else if (save.se.obj == FALSE) {
        printColoredMessage(
            message = '-- Save the plots as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The plotPCVariableRegression function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) == 1){
            return(all.pc.var.reg.plots = list(all.reg.plots = all.reg.plots))
        } else{
            return(all.pc.var.reg.plots = list(
                all.reg.plots = all.reg.plots,
                overall.reg.plot = overall.reg.plot))
        }

    }
}

