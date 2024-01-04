#' is used to compute ANOVA between individual gene expression of the assays in a SummarizedExperiment object
#' and a categorical variable.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the ANOVA. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param method A character string indicating which method
#' is to be used for the differential analysis: "aov" or "welch.correction". By default "aov" will
#' be selected.
#' @param plot.top.genes Logical. Indicates whether to plot the gene expression of the number of genes
#' from the top listing of anova by default it is set to FALSE.
#' @param nb.top.genes Defines the number of genes from the top or bottom listing of anova to plot,
#' by default is set to 3.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data. By default
#' the log transformation will be selected.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na TO BE DEFINED.
#' @param apply.round Logical. Indicates whether to round the ARI results, by default it is set to TRUE.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param plot.output Logical. Indicates whether to plot the boxplot of the F-test statistics, by default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the log2 F-statistics of ANOVA on the continuous variable
#' and if requested the associated boxplot.
#'
#'
#' @importFrom SummarizedExperiment assays assay
#' @importFrom matrixTests row_oneway_equalvar row_oneway_welch
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom kunstomverse geom_boxplot2
#' @importFrom stats anova
#' @import ggplot2
#'
#' @export
genesVariableAnova <- function(
        se.obj,
        assay.names = 'All',
        variable,
        method = 'aov',
        plot.top.genes = FALSE,
        nb.top.genes = 3,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        remove.na = 'both',
        apply.round = TRUE,
        save.se.obj = TRUE,
        plot.output = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The genesVariableAnova function starts:',
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
        stop('The variable cannot be empty.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0('The ', variable, ' should be a categorical variable.'))
    } else if (length(unique(se.obj@colData[, variable])) < 2) {
        stop(paste0('The ', variable, ', contains only one level. ANOVA cannot be performed.'))
    }
    if (!method %in% c('aov', 'welch.correction')) {
        stop('The method should be one of the "aov" or "welch.correction".')
    }
    if (apply.log){
        if(pseudo.count <0 )
            stop('The value of pseudo.count cannot be negative.')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'All') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- as.factor(unlist(assay.names))

    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = remove.na,
            verbose = verbose)
    }
    # anova ####
    all.aov <- lapply(
        levels(assay.names),
        function(x) {
            # data transformation ####
            printColoredMessage(message = '-- Data transformation:',
                                color = 'magenta',
                                verbose = verbose)
            if (apply.log) {
                temp.data <- log2(assay(x = se.obj, i = x) +  pseudo.count)
                printColoredMessage(message = paste0(
                    'Perform log2 + ',
                    pseudo.count,
                    '(pseudo.count) on the ',
                    x,
                    ' data.'),
                    color = 'blue',
                    verbose = verbose)
            } else{
                printColoredMessage(
                    message = paste0('Please make sure that the ', x, ' data is log transformed.'),
                    color = 'blue',
                    verbose = verbose)
                temp.data <- assay(x = se.obj, i = x)
            }
            # anova ####
            if (method == 'aov') {
                printColoredMessage(
                    message = paste0(
                        'Perform ANOVA with equal variance between individual genes expression of the ',
                        x,
                        ' data and the ',
                        variable,
                        ' variable.'),
                    color = 'blue',
                    verbose = verbose)
                anova.genes.var <- row_oneway_equalvar(x = temp.data, g = se.obj@colData[, variable])
            } else if (method == 'welch.correction') {
                printColoredMessage(
                    message = paste0(
                        'Perform ANOVA with Welch correction between individual genes expression of the ',
                        x,
                        ' data and the ',
                        variable,
                        ' variable.'),
                    color = 'blue',
                    verbose = verbose)
                anova.genes.var <- row_oneway_welch(x = temp.data, g = se.obj@colData[, variable])
            }
            row.names(anova.genes.var) <- row.names(se.obj)

            # round the anova statistic obtained to 2 digits ####
            if (apply.round) {
                anova.genes.var <- cbind(
                    round(anova.genes.var[, 1:9], digits = 3),
                    anova.genes.var[, 10, drop = FALSE])
            }
            # plot top genes ####
            if (plot.top.genes) {
                temp.anova <- anova.genes.var[order(anova.genes.var[, 'statistic'],
                                                    decreasing = TRUE,
                                                    na.last = TRUE) , ]
                var <- NULL
                p.high <- as.data.frame(t(temp.data[row.names(temp.anova)[c(1:nb.top.genes)],]))
                p.high <- mutate(p.high , var = se.obj@colData[, variable])
                p.high <- pivot_longer(
                    data = p.high,
                    cols =  -var,
                    names_to = 'genes',
                    values_to = 'expr')
                p.high <- ggplot(p.high, aes(x = var, y = expr)) +
                    geom_boxplot() +
                    ylab(expression(Log[2] ~ 'gene expression')) +
                    xlab(variable) +
                    facet_wrap( ~ genes) +
                    ggtitle(paste0(
                        nb.top.genes,
                        " Top affected genes by the variable ",
                        variable,
                        " for ",
                        x)) +
                    theme(panel.background = element_blank(),
                          axis.line = element_line(colour = 'black', linewidth = 1),
                          axis.title.x = element_text(size = 14),
                          axis.title.y = element_text(size = 14),
                          axis.text.x = element_text(
                              size = 10,
                              angle = 45,
                              hjust = 1,
                              vjust = 1),
                          axis.text.y = element_text(size = 12),
                          legend.text = element_text(size = 10),
                          legend.title = element_text(size = 14),
                          strip.text.x = element_text(size = 10),
                          plot.title = element_text(size = 12))
                plot(p.high)
                rm(temp.data)
                rm(temp.anova)
            }
            results <- NULL
            results <- list(anova.genes.var =  anova.genes.var)
            return(results)
        })
    names(all.aov) <- levels(assay.names)
    # save the results ####
    printColoredMessage(
        message = '-- Save the ANOVA results :',
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
            if (!paste0('gene.', method, '.anova') %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][[paste0('gene.', method, '.anova')]] <-
                    list()
            }
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][[paste0('gene.', method, '.anova')]][[variable]] <-
                log2(all.aov[[x]][['anova.genes.var']][, 'statistic'])
        }
        printColoredMessage(
            message = 'The ANOVA results for indiviaul assay are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)

        ## Plot and save the plot into se.obj@metadata$plot
        if (plot.output == TRUE) {
            printColoredMessage(
                message = '-- Plot the ANOVA results:',
                color = 'magenta',
                verbose = verbose
            )
            printColoredMessage(
                message = 'A boxplot of the ANOVA results are saved to metadata@plot.',
                color = 'blue',
                verbose = verbose)
            se.obj <- plotMetric(
                se.obj,
                assay.names = assay.names,
                metric = paste0('gene.', method, '.anova'),
                variable = variable,
                verbose = verbose)
        }
        printColoredMessage(message = '------------The genesVariableAnova function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
        ## return the results as a list ####
    } else if (save.se.obj == FALSE) {
        printColoredMessage(
            message = 'The ANOVA results for indiviaul assay are saved as list.',
            color = 'blue',
            verbose = verbose)
        all.aov <- lapply(
            levels(assay.names),
            function(x) log2(all.aov[[x]][['anova.genes.var']][, 'statistic']))
        names(all.aov) <- levels(assay.names)
        printColoredMessage(message = '------------The genesVariableAnova function finished.',
                            color = 'white',
                            verbose = verbose)
        return(genes.var.anova = all.aov)
    }
}
