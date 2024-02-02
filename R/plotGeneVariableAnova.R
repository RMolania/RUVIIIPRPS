#' plot F-statistics obtained from ANOVA.

#' @author Ramyar Molania

#' @description
#' This functions computes the adjusted rand index for given a categorical variable using the first PCs of the assay(s)
#' in a SummarizedExperiment object.


#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to compute PCA. By default all the assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. Indicates the column name in the SummarizedExperiment object that contains a categorical
#' variable such as sample types or batches.
#' @param anova.method Logical. Indicates whether to use the PCA calculated using a specific number of PCs instead of the
#' full range to speed up the process, by default is set to 'TRUE'.
#' @param plot.output Logical. Indicates whether to plot the ARI, by default it is set to FALSE.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list that containing the computed ARI on the categorical variable.

#' @importFrom SummarizedExperiment assays assay
#' @import ggplot2
#' @export

plotGenesVariableAnova <- function(
        se.obj,
        assay.names = 'all',
        variable,
        anova.method = 'genes.aov.anova',
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
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
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
                stop(paste0('Any ANOVA analysis has not been computed yet on the  ', x, ' assay'))
        })
    # obtain correlations coeff ####
    all.aov.fvals <- lapply(
        levels(assay.names),
        function(x) {
            if (!'aov' %in% names(se.obj@metadata[['metric']][[x]])) {
                stop(paste0('Any ANOVA analysis has not been computed yet on the  ', x, ' assay'))
            }
            if (!anova.method %in% names(se.obj@metadata[['metric']][[x]][['aov']])) {
                stop(paste0('The ', anova.method , ' has not been computed yet for the ', x, ' assay.'))
            }
            if (!variable %in% names(se.obj@metadata[['metric']][[x]][['aov']][[anova.method]])) {
                stop(paste0('The ', anova.method , ' has not been computed yet for the ', variable, ' variable and the ', x, ' assay.'))
            }
            se.obj@metadata[['metric']][[x]][['aov']][[anova.method]][[variable]]$fvals
        })
    names(all.aov.fvals) <- levels(assay.names)
    aov.fvals <- NULL
    all.aov.fvals.plots <- lapply(
        levels(assay.names),
        function(x){
            aov.fvals <- all.aov.fvals[[x]]
            p.corr.coeff <- ggplot() +
                geom_boxplot(aes(y = aov.fvals, x = 1)) +
                ggtitle(paste0('ANOVA, ', variable)) +
                xlab(x) +
                ylab(expression(Log[2]~'F-statistic')) +
                geom_hline(yintercept = 0) +
                theme(panel.background = element_blank(),
                      axis.line = element_line(colour = 'black', linewidth = 1),
                      axis.title.x = element_text(size = 18),
                      axis.title.y = element_text(size = 18),
                      plot.title = element_text(size = 15),
                      axis.ticks.x = element_blank(),
                      axis.text.x = element_text(size = 0),
                      axis.text.y = element_text(size = 12))
        })
    names(all.aov.fvals.plots) <- levels(assay.names)

    # overall plot ####
    everything <- datasets <- NULL
    if(length(assay.names) > 1){
        all.aov.fvals <- as.data.frame(all.aov.fvals) %>%
            tidyr::pivot_longer(
                everything(),
                names_to = 'datasets',
                values_to = 'aov.fvals')
        overall.aov.fvals.plot <- ggplot(all.aov.fvals, aes(x = datasets, y = aov.fvals)) +
            geom_boxplot() +
            ggtitle('ANOVA') +
            xlab('Datasets') +
            ylab('F-statistic') +
            geom_hline(yintercept = 0) +
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
            se.obj@metadata[['metric']][[x]][['genes.var.corr']][[anova.method]][[variable]]$corrs.plot <- all.aov.fvals.plots[[x]]
        }
        printColoredMessage(
            message = 'The F values plots for indiviaul assay are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)

        if (length(assay.names) > 1) {
            if (!'plot' %in%  names(se.obj@metadata)) {
                se.obj@metadata[['plot']] <- list()
            }
            if (!'GeneVarAnova' %in%  names(se.obj@metadata[['plot']])) {
                se.obj@metadata[['plot']][['GeneVarAnova']] <- list()
            }
            if (!anova.method %in%  names(se.obj@metadata[['plot']][['GeneVarAnova']])) {
                se.obj@metadata[['plot']][['GeneVarAnova']][[anova.method]] <- list()
            }
            if (!variable %in%  names(se.obj@metadata[['plot']][['GeneVarAnova']][[anova.method]])) {
                se.obj@metadata[['plot']][['GeneVarAnova']][[anova.method]][[variable]] <- list()
            }
            se.obj@metadata[['plot']][['GeneVarAnova']][[anova.method]][[variable]] <- overall.aov.fvals.plot
            printColoredMessage(
                message = paste0('The F values plots all assays are saved to metadata@plot'),
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
            return(all.aov.fvals.plots = all.aov.fvals.plots)
        } else{
            return(gene.var.corr.plot = list(
                all.aov.fvals.plots = all.aov.fvals.plots,
                overall.aov.fvals.plot = overall.aov.fvals.plot))
        }
    }
}
