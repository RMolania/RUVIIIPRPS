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
#' @param variable Symbol. Indicates the column name in the SummarizedExperiment object that contains a categorical
#' variable such as sample types or batches.
#' @param correlation.method TTTT
#' @param boxplot.color TTTT
#' @param geom.hline.color TTTT
#' @param correlation.method TTTT
#' @param plot.output Logical. Indicates whether to plot the ARI, by default it is set to FALSE.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object, by default it is
#' set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the
#' execution of the functions, by default it is set to TRUE.

#' @return A SummarizedExperiment object or a list that containing the computed ARI on the categorical variable.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export

plotGenesVariableCorrelation <- function(
        se.obj,
        assay.names = 'all',
        variable,
        correlation.method = 'gene.spearman.corr',
        plot.output = TRUE,
        boxplot.color = 'black',
        geom.hline.color = 'gray',
        assess.se.obj = TRUE,
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
    } else assay.names <- as.factor(unlist(assay.names))

    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = 'none',
            verbose = verbose)
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
            if (!'genes.var.corr' %in% names(se.obj@metadata[['metric']][[x]])) {
                stop(paste0('Any correlation analysis has not been computed yet on the  ', x, ' assay'))
            }
            if (!correlation.method %in% names(se.obj@metadata[['metric']][[x]][['genes.var.corr']])) {
                stop(paste0('The ', correlation.method , ' has not been computed yet for the ', x, ' assay.'))
            }
            if (!variable %in% names(se.obj@metadata[['metric']][[x]][['genes.var.corr']][[correlation.method]])) {
                stop(paste0('The ', correlation.method , ' has not been computed yet for the ', variable, ' variable and the ', x, ' assay.'))
            }
            corr.coeff <- se.obj@metadata[['metric']][[x]][['genes.var.corr']][[correlation.method]][[variable]]$corrs
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
                geom_boxplot(aes(y = corr.coeff, x = 1), colour = boxplot.color) +
                ylab('Spearman correlation coefficient') +
                xlab(x) +
                geom_hline(yintercept = 0, colour = geom.hline.color) +
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
            geom_boxplot(colour = boxplot.color) +
            ylab('Spearman correlation coefficient') +
            xlab('Datasets') +
            geom_hline(yintercept = 0, colour = geom.hline.color) +
            ggtitle('Spearman correlation analysis') +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                plot.title = element_text(size = 15),
                axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12))
    }
    # save the results ####
    printColoredMessage(
        message = '-- Save the all the plots:',
        color = 'magenta',
        verbose = verbose)
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['genes.var.corr']][[correlation.method]][[variable]]$corrs.plot <-
                all.corr.coeff.plots[[x]]
        }
        printColoredMessage(
            message = 'The correlation plots for indiviaul assay are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)

        if (length(assay.names) > 1) {
            if (!'plot' %in%  names(se.obj@metadata)) {
                se.obj@metadata[['plot']] <- list()
            }
            if (!'GeneVarCorr' %in%  names(se.obj@metadata[['plot']])) {
                se.obj@metadata[['plot']][['GeneVarCorr']] <- list()
            }
            if (!correlation.method %in%  names(se.obj@metadata[['plot']][['GeneVarCorr']])) {
                se.obj@metadata[['plot']][['GeneVarCorr']][[correlation.method]] <- list()
            }
            if (!variable %in%  names(se.obj@metadata[['plot']][['GeneVarCorr']][[correlation.method]])) {
                se.obj@metadata[['plot']][['GeneVarCorr']][[correlation.method]][[variable]] <- list()
            }
            se.obj@metadata[['plot']][['GeneVarCorr']][[correlation.method]][[variable]] <- overall.corr.coeff.plot

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
