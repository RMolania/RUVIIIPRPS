#' Plot p-values histograms of DEG analysis.

#' @author Ramyar Molania

#' @description
#' This function plots the p-values histograms of the DEG analysis.

#' @details
#' DE analyses is performed using the Wilcoxon signed-rank test with log-transformed data e.g. raw counts, normalized data, ....
#' To evaluate the effects of the different sources of unwanted variation on the data, DE analyses is performed across
#' batches. In the absence of any batch effects, the histogram of the resulting unadjusted P values should be uniformly
#' distributed

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or vector of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to compute DGE. The default is set to 'all' indicating all the the assays of the
#' SummarizedExperiment object will be selected.
#' @param variable Symbol. Specifies the column name in the SummarizedExperiment object containing a categorical variable
#' , such as sample types or batches
#' @param plot.ncol Numeric. Indicates number of columns in the plot grid.
#' @param plot.output Logical. Indicates whether to plot the p-values histograms when the functions is running, by default
#' it is set to 'FALSE'.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object or a list that containing the computed ARI on the categorical variable.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#' @export

plotDGE <- function(
        se.obj,
        assay.names = 'all',
        variable,
        plot.ncol = 1,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
){
    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    }
    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # check metric ####
    m.out <- lapply(
        levels(assay.names),
        function(x) {
            if (!x %in% names(se.obj@metadata[['metric']]))
                stop(paste0('Any DEA analysis has not been computed yet on the  ', x, ' assay'))
        })
    # obtain tests ####
    printColoredMessage(
        message = '-- Obtained computed p-values of the DGE analysis from the SummarizedExperiment object:',
        color = 'magenta',
        verbose = verbose
    )
    all.de.tests <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('-Obtain computed p-values of the DGE analysis for the', x, 'data.'),
                color = 'blue',
                verbose = verbose
            )
            if (!'DGE' %in% names(se.obj@metadata[['metric']][[x]])) {
                stop(paste0('Any DGE has not been computed yet for the ', x, ' assay.'))
            }
            if (!variable %in% names(se.obj@metadata[['metric']][[x]][['DGE']]) ) {
                stop(paste0('The DGE has not been computed yet for the ', variable, ' variable and the ', x, ' assay.'))
            }
            de.results <- do.call(cbind, se.obj@metadata[['metric']][[x]][['DGE']][[variable]][['p.vals']])
            de.results <- as.data.frame(de.results[ , seq(3, ncol(de.results), 3), drop = FALSE])
            de.results
        })
    names(all.de.tests) <- levels(assay.names)
    p.vals <- everything <- NULL
    all.pval.plots <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0('-Generate p-values p-values histograms for the', x, 'data.'),
                color = 'blue',
                verbose = verbose
            )
            if(length(unique(colData(se.obj)[[variable]])) == 2 ){
                pval.data <- all.de.tests[[x]]
                colnames(pval.data) <- 'p.vals'
                pval.plot <- ggplot(pval.data, aes(x = p.vals)) +
                    geom_histogram(binwidth = 0.1) +
                    ylab('Frequency') +
                    xlab('p-values') +
                    ggtitle(paste0('p-value histograms, ', x)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 18),
                        axis.title.y = element_text(size = 18),
                        plot.title = element_text(size = 15),
                        axis.text.x = element_text(size = 8),
                        axis.text.y = element_text(size = 8))
            } else {
                pval.data <- all.de.tests[[x]] %>%
                    tidyr::pivot_longer(
                        everything(),
                        names_to = 'contrasts',
                        values_to = 'p.vals')
                pval.plot <- ggplot(pval.data, aes(x = p.vals)) +
                    geom_histogram(binwidth = 0.1) +
                    ylab('Frequency') +
                    xlab('p-values') +
                    facet_wrap(~contrasts) +
                    ggtitle(paste0('p-value histograms, ', x)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 10),
                        axis.title.y = element_text(size = 10),
                        plot.title = element_text(size = 12),
                        axis.text.x = element_text(size = 8),
                        axis.text.y = element_text(size = 8))
            }
            return(pval.plot)
        })
    names(all.pval.plots) <-  levels(assay.names)

    # overall plots
    printColoredMessage(
        message = '-Generate p-values histograms for all the assays.' ,
        color = 'blue',
        verbose = verbose)
    overall.pvals.plots <- ggarrange(plotlist = all.pval.plots, ncol = plot.ncol)
    if(plot.output)
        print(overall.pvals.plots)
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['DGE']][[variable]][['p.vals.plot']] <- all.pval.plots[[x]]
        }
        printColoredMessage(
            message = 'The Wilcoxon results for indiviaul assay are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)


        printColoredMessage(message = '------------The genesDEA function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) > 1){
            if (!'plot' %in%  names(se.obj@metadata)) {
                se.obj@metadata[['plot']] <- list()
            }
            if (!'DEG' %in%  names(se.obj@metadata[['plot']])) {
                se.obj@metadata[['plot']][['DEG']] <- list()
            }
            if (!variable %in%  names(se.obj@metadata[['plot']][['DEG']])) {
                se.obj@metadata[['plot']][['DEG']][[variable]] <- list()
            }
            se.obj@metadata[['plot']][['DEG']][[variable]] <- overall.pvals.plots
            printColoredMessage(
                message = paste0('The RLE plots of all assays are saved to metadata@plot'),
                color = 'blue',
                verbose = verbose
            )
        }
        return(se.obj = se.obj)
        ## return the results as a list ####
    } else if (save.se.obj == FALSE) {
        printColoredMessage(
            message = 'All the p-value histograms are outputed as list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The genesDEA function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) == 1){
            return(all.pval.plots = all.pval.plots)
        } else{
            return(all.pval.plots = all.pval.plots, overall.pvals.plots = overall.pvals.plots)
        }
    }
}
