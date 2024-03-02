#' Plot Spearman or Pearson correlations coefficients.

#' @author Ramyar Molania

#' @description
#' This function generates boxplots of computed Spearman or Pearson correlations coefficients of indivdial assays in a
#' SummarizedExperiment object

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or vector of symbols for the selection of the name(s) of the assay(s) of the
#' SummarizedExperiment object to compute the correlation. By default all the assays of the SummarizedExperiment class
#' object will be selected.
#' @param variable Symbol. Indicates a name in the columns of the sample annotation in the SummarizedExperiment object
#' that contains a continuous variable such as library size, tumor purity, ... .
#' @param correlation.method Symbol. Indicates which computed correlation coefficient should be used for plotting. The
#' default is 'gene.spearman.corr'. We refer to the 'computeGenesVariableCorrelation' function for more details.
#' @param plot.output Logical. Indicates whether to plot the ARI, by default it is set to FALSE.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list that containing the boxplots of the Spearman or Pearson correlations
#' coefficients for individual assays.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr everything
#' @import ggplot2
#' @export

plotGenesVariableCorrelation <- function(
        se.obj,
        assay.names = 'all',
        variable,
        correlation.method = 'spearman',
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
){
    # check the inputs ####
    if (is.null(assay.names)) {
        stop('Please provide at least an assay name.')
    } else if (is.null(variable)) {
        stop('Please provide a variable.')
    } else if (length(unique(se.obj@colData[, variable])) < 2) {
        stop(paste0('The ', variable, ', contains only one variable.'))
    } else if (!class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0(
            'The ',
            variable,
            ', is a numeric, but this should a categorical variable'
        ))
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
                stop(paste0('Any correlation analysis has not been computed yet on the  ', x, ' assay'))
        })
    # obtain correlations coeff ####
    printColoredMessage(
        message = paste0('-- Obtain the correlation coefficient for the from the SummarizedExperiment object:'),
        color = 'magenta',
        verbose = verbose
    )
    all.corr.coeff <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('-Obtain computed correlation coefficients the', x, 'data.'),
                color = 'blue',
                verbose = verbose
            )
            if (!'Correlation' %in% names(se.obj@metadata[['metric']][[x]])) {
                stop(paste0('Any correlation analysis has not been computed yet on the  ', x, ' assay'))
            }
            if (!correlation.method %in% names(se.obj@metadata[['metric']][[x]][['Correlation']])) {
                stop(paste0('The ', correlation.method , ' has not been computed yet for the ', x, ' assay.'))
            }
            if (!variable %in% names(se.obj@metadata[['metric']][[x]][['Correlation']][[correlation.method]])) {
                stop(paste0('The ', correlation.method , ' has not been computed yet for the ', variable, ' variable and the ', x, ' assay.'))
            }
            corr.coeff <- se.obj@metadata[['metric']][[x]][['Correlation']][[correlation.method]][[variable]]$cor.coef
        })
    names(all.corr.coeff) <- levels(assay.names)

    printColoredMessage(
        message = '--Generate boxplots:',
        color = 'magenta',
        verbose = verbose
    )
    all.corr.coeff.plots <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0('-Generate boxplot of the computed correlation coefficients the ', x, ' data.'),
                color = 'blue',
                verbose = verbose
            )
            corr.coeff <- all.corr.coeff[[x]]
            p.corr.coeff <- ggplot() +
                geom_boxplot(aes(y = corr.coeff, x = 1)) +
                ylab('Spearman correlation coefficient') +
                xlab(x) +
                geom_hline(yintercept = 0) +
                ggtitle('Spearman correlation analysis') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    plot.title = element_text(size = 15),
                    axis.ticks.x = element_blank(),
                    axis.text.x = element_text(size = 0),
                    axis.text.y = element_text(size = 12))
            if(isTRUE(plot.output) & length(assay.names) == 1) print(p.corr.coeff)
            return(p.corr.coeff)
        })
    names(all.corr.coeff.plots) <- levels(assay.names)

    # overall plot ####
    everything <- datasets <- corr.coff <- NULL
    if(length(assay.names) > 1){
        printColoredMessage(
            message = '-- Generate boxplots for all the assays.',
            color = 'blue',
            verbose = verbose
        )
        all.corr.coeff <- as.data.frame(all.corr.coeff) %>%
            tidyr::pivot_longer(
                dplyr::everything(),
                names_to = 'datasets',
                values_to = 'corr.coff')
        overall.corr.coeff.plot <- ggplot(all.corr.coeff, aes(x = datasets, y = corr.coff)) +
            geom_boxplot() +
            ylab('Spearman correlation coefficients') +
            xlab('Datasets') +
            geom_hline(yintercept = 0) +
            ggtitle(paste0('Spearman correlation analysis, ', variable)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                plot.title = element_text(size = 15),
                axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12))
        if(isTRUE(plot.output)) print(overall.corr.coeff.plot)
    }
    # save the results ####
    printColoredMessage(
        message = '-- Save the all the plots:',
        color = 'magenta',
        verbose = verbose)
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist for this assay, this metric and this variable
            if(!'cor.coef.plot' %in% se.obj@metadata[['metric']][[x]][['Correlation']][[correlation.method]][[variable]])
                se.obj@metadata[['metric']][[x]][['Correlation']][[correlation.method]][[variable]][['cor.coef.plot']] <- list()
            se.obj@metadata[['metric']][[x]][['Correlation']][[correlation.method]][[variable]][['cor.coef.plot']] <- all.corr.coeff.plots[[x]]
        }
        printColoredMessage(
            message = 'The correlation plots for indiviaul assay are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)

        if (length(assay.names) > 1) {
            if (!'plot' %in%  names(se.obj@metadata)) {
                se.obj@metadata[['plot']] <- list()
            }
            if (!'Correlation' %in%  names(se.obj@metadata[['plot']])) {
                se.obj@metadata[['plot']][['Correlation']] <- list()
            }
            if (!correlation.method %in%  names(se.obj@metadata[['plot']][['Correlation']])) {
                se.obj@metadata[['plot']][['Correlation']][[correlation.method]] <- list()
            }
            if (!variable %in%  names(se.obj@metadata[['plot']][['Correlation']][[correlation.method]])) {
                se.obj@metadata[['plot']][['Correlation']][[correlation.method]][[variable]] <- list()
            }
            se.obj@metadata[['plot']][['Correlation']][[correlation.method]][[variable]] <- overall.corr.coeff.plot

            printColoredMessage(
                message = paste0('The correlation plots all assays are saved to metadata@plot'),
                color = 'blue',
                verbose = verbose)

        }
        printColoredMessage(message = '------------The plotGenesVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)

    } else if (save.se.obj == FALSE) {
        # return only the correlation result ####
        printColoredMessage(
            message = 'All the plots are saved as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The genesVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(assay.names) == 1){
            return(all.corr.coeff.plots = all.corr.coeff.plots)
        } else{
            return(gene.var.corr.plot = list(
                all.corr.coeff.plots = all.corr.coeff.plots,
                overall.corr.coeff.plot = overall.corr.coeff.plot
                ))
        }
    }

}
