#' is used to compute the correlation between the gene expression (assay)
#' of a SummarizedExperiment class object and a continuous variable (i.e. library size).
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the correlation.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the correlation. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a continuous variable such as
#' library size from colData(se.obj).
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data. By default
#' the log transformation will be selected.
#' @param method A character string indicating which correlation coefficient
#' is to be used for the test: "pearson", "kendall", or "spearman". By default "spearman" will
#' be selected.
#' @param a The significance level used for the confidence intervals in the correlation,
#' by default it is set to 0.05.
#' @param rho The value of the hypothesised correlation to be used in the hypothesis testing,
#' by default it is set to 0.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param plot.output Logical. Indicates whether to plot the boxplot of the correlation, by default it is set to TRUE.
#' @param plot.top.genes Logical. Indicates whether to plot the gene expression of the number of genes
#' from the high or low correlation, by default it is set to FALSE.
#' @param nb.top.genes Defines the number of genes from the high or low correlation to plot,
#' by default is set to 3.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na TO BE DEFINED.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param apply.round Logical. Indicates whether to round the ARI results, by default it is set to TRUE.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the computed correlation on the continuous variable
#' and if requested the associated plot.
#' @importFrom SummarizedExperiment assays assay colData
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom Rfast correls
#' @importFrom stats var cor.test
#' @import ggplot2
#' @export

genesVariableCorrelation <- function(
        se.obj,
        assay.names = 'All',
        variable,
        method = 'spearman',
        a = 0.05,
        rho = 0,
        plot.output = TRUE,
        plot.top.genes = FALSE,
        nb.top.genes = 3,
        apply.log = TRUE,
        pseudo.count = 1,
        apply.round = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The genesVariableCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)
    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be null.')
    }
    if(length(assay.names) == 1 & assay.names!= 'All'){
        if(!assay.names %in% names(assays(se.obj)) )
            stop('The assay name cannot be found in the SummarizedExperiment object.')
    }
    if(length(assay.names) > 1){
        if(sum(!assay.names %in% names(assays(se.obj))) > 0 )
            stop('The assay names cannot be found in the SummarizedExperiment object.')
    }
    if (is.null(variable)) {
        stop('Please provide a variable.')
    } else if (!class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0('The ', variable, 'should be a continuous variable.'))
    } else if (!method %in% c('pearson', 'spearman')) {
        stop('The method should be "pearson" or "spearman".')
    }
    if (var(se.obj[[variable]]) == 0) {
        stop(paste0('The ', variable, ' appears to have no variation.'))
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'All') {
        assay.names = as.factor(names(assays(se.obj)))
    } else assay.names = as.factor(unlist(assay.names))

    # check SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = remove.na,
            verbose = verbose)
    }

    # correlation analyses ####
    printColoredMessage(message = '-- Correlation analyses:',
                        color = 'magenta',
                        verbose = verbose)
    cor.all <- lapply(
        levels(assay.names),
        function(x) {
            # data transformation ####
            if (apply.log) {
                temp.data <- log2(assay(x = se.obj, i = x) +  pseudo.count)
                printColoredMessage(
                    message = paste0('Perform log2 + ', pseudo.count, '(pseudo.count) on the assay ', x , '.'),
                    color = 'blue',
                    verbose = verbose)
            } else temp.data <- assay(x = se.obj, i = x)

            # correlation ####
            printColoredMessage(
                message = paste0(
                        'Perform ' ,
                        method,
                        ' correlation between individual genes expression of the ',
                        x,
                        ' and the ',
                        variable,
                        '.'),
                color = 'blue',
                verbose = verbose
            )
            corr.genes.var <- correls(
                y = se.obj@colData[, variable],
                x = t(temp.data),
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
            if (plot.top.genes) {
                printColoredMessage(
                    message = paste0(
                        '-- Plot top ' ,
                        nb.top.genes,
                        ' highly correlated (positive and negative) genes with the ',
                        variable,
                        '.'),
                    color = 'magenta',
                    verbose = verbose
                )
                temp.corr <- corr.genes.var[order(corr.genes.var[, 'correlation'],
                                         decreasing = TRUE,
                                         na.last = TRUE) ,]
                ### positive correlation ####
                p.pos <- as.data.frame(t(temp.data[row.names(temp.corr)[c(1:nb.top.genes)],]))
                p.pos[, 'variable'] <- se.obj@colData[, variable]
                p.pos <-p.pos %>% pivot_longer(
                    -variable,
                    names_to = 'genes',
                    values_to = 'expr')
                p.pos <- ggplot(p.pos, aes(x = variable, y = expr)) +
                    geom_point() +
                    ylab(expression(Log[2] ~ 'gene expression')) +
                    xlab(variable) +
                    geom_smooth(method = 'lm', formula = y~x, colour = 'red') +
                    facet_wrap(~ genes) +
                    ggtitle(paste0(
                            "Top highly positively correlated genes with ",
                            variable,
                            '\n in the assay ',
                            x)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', size = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 10),
                        axis.text.y = element_text(size = 12),
                        strip.text.x = element_text(size = 10),
                        plot.title = element_text(size = 14)
                    )
                plot(p.pos)

                ### negative correlation ####
                temp.corr <- corr.genes.var[order(corr.genes.var[, 'correlation'],
                                         decreasing = FALSE,
                                         na.last = TRUE) ,]
                p.neg <- as.data.frame(t(temp.data[row.names(temp.corr)[c(1:nb.top.genes)],]))
                p.neg[, 'variable'] <- se.obj@colData[, variable]
                p.neg <- p.neg %>% pivot_longer(
                    -variable,
                    names_to = 'genes',
                    values_to = 'expr')
                p.neg <- ggplot(p.neg, aes(x = variable, y = expr)) +
                    geom_point() +
                    ylab(expression(Log[2] ~ 'gene expression')) +
                    xlab(variable) +
                    geom_smooth(method = 'lm', formula = y~x, colour = 'red') +
                    facet_wrap(~ genes) +
                    ggtitle(
                        paste0(
                            'Top highly negatively correlated genes with the variable ',
                            variable,
                            '\n in the assay ',
                            x)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', size = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 10),
                        axis.text.y = element_text(size = 12),
                        strip.text.x = element_text(size = 10),
                        plot.title = element_text(size = 14)
                    )
                plot(p.neg)
                rm(temp.data)
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
            if (!paste0('gene.', method, '.corr') %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][[paste0('gene.', method, '.corr')]] <- list()
            }
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][[paste0('gene.', method, '.corr')]][[variable]] <-
                cor.all[[x]][['corr.genes.var']][, 'correlation']
        }
        printColoredMessage(
            message = 'The correlation results for indiviaul assay are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)

        ## save a plot to SummarizedExperiment ####
        if (plot.output == TRUE) {
            printColoredMessage(
                message = '-- Plot the correlation results:',
                color = 'magenta',
                verbose = verbose)
            printColoredMessage(
                message = 'A boxplot of the correlation results are saved to metadata@plot.',
                color = 'blue',
                verbose = verbose)
            se.obj <- plotMetric(
                se.obj,
                assay.names = assay.names,
                metric = paste0('gene.', method, '.corr'),
                variable = variable,
                verbose = verbose
            )
        }
        return(se.obj = se.obj)

        # return only the correlation result ####
    } else if (save.se.obj == FALSE) {
        printColoredMessage(
            message = 'The correlation results for indiviaul assay are saved as list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The genesVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(gene.corr.var = cor.all)
    }
}
