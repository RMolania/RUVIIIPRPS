#' is used to compute the correlation between individual gene expression and a continuous variable.

#' @author Ramyar Molania

#' @description
#' This function computes Spearman or Pearson correlations between individual gene-level expression of each assay and
#' a continuous variable in SummarizedExperiment object.

#' @details
#' Additional details...

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or vector of symbols for the selection of the name(s) of the assay(s) of the
#' SummarizedExperiment object to compute the correlation. By default all the assays of the SummarizedExperiment class
#' object will be selected.
#' @param variable Symbol. Indicates a column name of the SummarizedExperiment object that contains a continuous variable
#' such as library size, tumor purity, ....
#' @param method Symbol. Indicates which correlation methods should be used. The options are 'pearson', 'kendall', or
#' "spearman". The default is 'spearman'.
#' @param a Numeric. The significance level used for the confidence intervals in the correlation, by default it is set to
#' 0.05. We refer to the correls function from the Rfast package for more details.
#' @param rho Numeric. The value of the hypothesized correlation to be used in the hypothesis testing, by default it is
#' set to 0. We refer to the correls function from the Rfast package for more details.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before computing the correlation.
#' By default the log transformation will be selected.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param plot.top.genes Logical. Indicates whether to plot the gene expression of the number of genes from the high or
#' low correlation, by default it is set to FALSE.
#' @param nb.top.genes Numeric. Defines the number of genes that show the highest or lowest correlation with variable to
#' plot. The default is 3.
#' @param apply.round Logical. Indicates whether to round the correlation coefficients. The default it TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object. We refer to the
#' checkSeObj function for more details. The default is TRUE.
#' @param remove.na To remove NA or missing values from the assays and variable.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the current SummarizedExperiment
#' object or to output the result. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list that contains the  correlation coefficients on the continuous variable.

#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @importFrom SummarizedExperiment assays assay colData
#' @importFrom Rfast correls
#' @importFrom stats var
#' @importFrom tidyr pivot_longer %>%
#' @import ggplot2
#' @export

computeGenesVariableCorrelation <- function(
        se.obj,
        assay.names = 'all',
        variable,
        method = 'spearman',
        a = 0.05,
        rho = 0,
        apply.log = TRUE,
        pseudo.count = 1,
        plot.top.genes = FALSE,
        nb.top.genes = 3,
        apply.round = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The computeGenesVariableCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)
    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (!is.vector(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s) or "assay.names = all".')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    } else if (!class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0('The ', variable, 'should be a continuous variable.'))
    } else if (!method %in% c('pearson', 'spearman')) {
        stop('The method must be one of the "pearson" or "spearman".')
    }
    if (var(se.obj[[variable]]) == 0) {
        stop(paste0('The ', variable, ' appears to have no variation.'))
    }


    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # check SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = remove.na,
            verbose = verbose)
    }

    # data transformation ####
    printColoredMessage(
        message = '-- Data transformation:',
        color = 'magenta',
        verbose = verbose)
    all.assays <- lapply(
        levels(assay.names),
        function(x){
            # log transformation ####
            if (apply.log & !is.null(pseudo.count)) {
                printColoredMessage(
                    message = paste0('Apply log2 + ', pseudo.count,  ' (pseudo.count) on the ', x, ' assay.'),
                    color = 'blue',
                    verbose = verbose)
                expr <- log2(assay(x = se.obj, i = x) + pseudo.count)
            } else if (apply.log & is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0('Apply log2 transformation on the ', x, ' assay.'),
                    color = 'blue',
                    verbose = verbose)
                expr <- log2(assay(x = se.obj, i = x))
            } else {
                printColoredMessage(
                    message = paste0('The ', x, ' assay will be used without log transformation.'),
                    color = 'blue',
                    verbose = verbose)
                printColoredMessage(
                    message = 'Please note, the assay should be in log scale before performing the correlation analysis.',
                    color = 'red',
                    verbose = verbose)
                expr <- assay(x = se.obj, i = x)
            }
        })
    names(all.assays) <- levels(assay.names)

    # correlation analyses ####
    printColoredMessage(message = '-- Correlation analyses:',
                        color = 'magenta',
                        verbose = verbose)
    cor.all <- lapply(
        levels(assay.names),
        function(x) {
            # correlation ####
            printColoredMessage(
                message = paste0('Perform ' , method,' correlation between individual genes expression of the ',
                        x, ' and the ', variable, '.'),
                color = 'blue',
                verbose = verbose)
            corr.genes.var <- correls(
                y = se.obj@colData[, variable],
                x = t(all.assays[[x]]),
                type = method,
                a = a ,
                rho = rho)
            row.names(corr.genes.var) <- row.names(se.obj)
            # round the correlation obtained to 2 digits ####
            if (apply.round) {
                corr.genes.var <- cbind(
                    round(x = corr.genes.var[, 1:4], digits = 2),
                    corr.genes.var[, 5, drop = FALSE])
            }
            # plot highly affected genes ####
            if (isTRUE(plot.top.genes)) {
                printColoredMessage(
                    message = paste0('Plot top ' , nb.top.genes,
                        ' highly correlated (positive and negative) genes with the ', variable,'.'),
                    color = 'blue',
                    verbose = verbose)
                temp.corr <- corr.genes.var[order(corr.genes.var[, 'correlation'],
                                         decreasing = TRUE,
                                         na.last = TRUE) , ]
                ### positive correlation ####
                p.pos <- as.data.frame(t(all.assays[[x]][row.names(temp.corr)[c(1:nb.top.genes)],]))
                p.pos[, 'variable'] <- se.obj@colData[, variable]
                p.pos <-p.pos %>%
                    pivot_longer(
                    -variable,
                    names_to = 'genes',
                    values_to = 'expr')
                p.pos <- ggplot(p.pos, aes(x = variable, y = expr)) +
                    geom_point() +
                    ylab(expression(Log[2] ~ 'gene expression')) +
                    xlab(variable) +
                    geom_smooth(method = 'lm', formula = y~x, colour = 'red') +
                    facet_wrap(~ genes) +
                    ggtitle(paste0("Top highly positively correlated genes with ", variable,'\n in the assay ', x)) +
                    theme(panel.background = element_blank(),
                          axis.line = element_line(colour = 'black', size = 1),
                          axis.title.x = element_text(size = 14),
                          axis.title.y = element_text(size = 14),
                          axis.text.x = element_text(size = 10),
                          axis.text.y = element_text(size = 12),
                          strip.text.x = element_text(size = 10),
                          plot.title = element_text(size = 14))
                print(p.pos)

                # negative correlation ####
                temp.corr <- corr.genes.var[order(corr.genes.var[, 'correlation'],
                                         decreasing = FALSE,
                                         na.last = TRUE) ,]
                p.neg <- as.data.frame(t(all.assays[[x]][row.names(temp.corr)[c(1:nb.top.genes)],]))
                p.neg[, 'variable'] <- se.obj@colData[, variable]
                p.neg <- p.neg %>%
                    tidyr::pivot_longer(
                        -variable,
                        names_to = 'genes',
                        values_to = 'expr')
                p.neg <- ggplot(p.neg, aes(x = variable, y = expr)) +
                    geom_point() +
                    ylab(expression(Log[2] ~ 'gene expression')) +
                    xlab(variable) +
                    geom_smooth(method = 'lm', formula = y~x, colour = 'red') +
                    facet_wrap(~ genes) +
                    ggtitle(paste0('Top highly negatively correlated genes with the variable ', variable, '\n in the assay ', x)) +
                    theme( panel.background = element_blank(),
                           axis.line = element_line(colour = 'black', size = 1),
                           axis.title.x = element_text(size = 14),
                           axis.title.y = element_text(size = 14),
                           axis.text.x = element_text(size = 10),
                           axis.text.y = element_text(size = 12),
                           strip.text.x = element_text(size = 10),
                           plot.title = element_text(size = 14))

                print(p.neg)
                rm(temp.corr)
            }
            results <- NULL
            results <- list(corr.genes.var = corr.genes.var)
            return(results)
        })
    names(cor.all) <- levels(assay.names)

    # save the results ####
    printColoredMessage(
        message = '-- Save the correlation results:',
        color = 'magenta',
        verbose = verbose)
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!'metric' %in% names(se.obj@metadata)) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!x %in% names(se.obj@metadata[['metric']])) {
                se.obj@metadata[['metric']][[x]] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!'genes.var.corr' %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['genes.var.corr']] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!paste0('gene.', method, '.corr') %in% names(se.obj@metadata[['metric']][[x]][['genes.var.corr']])) {
                se.obj@metadata[['metric']][[x]][['genes.var.corr']][[paste0('gene.', method, '.corr')]] <- list()
            }
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['genes.var.corr']][[paste0('gene.', method, '.corr')]][[variable]]$corrs <-
                cor.all[[x]][['corr.genes.var']][, 'correlation']
        }
        printColoredMessage(
            message = 'The correlation results for indiviaul assay are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)
        return(se.obj = se.obj)

        # return only the correlation result ####
    } else if (save.se.obj == FALSE) {
        printColoredMessage(
            message = 'The correlation results for indiviaul assay are saved as list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The computeGenesVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(gene.corr.var = cor.all)
    }
}
