#' is used to compute the correlation between the gene expression (assay)
#' of a SummarizedExperiment class object and a continuous variable (i.e. library size).
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the correlation.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the correlation. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a continuous variable such as
#' library size from colData(se.obj).
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' the log transformation will be selected.
#' @param method A character string indicating which correlation coefficient
#' is to be used for the test: "pearson", "kendall", or "spearman". By default "spearman" will
#' be selected.
#' @param a The significance level used for the confidence intervals in the correlation,
#' by default it is set to 0.05.
#' @param rho The value of the hypothesised correlation to be used in the hypothesis testing,
#' by default it is set to 0.
#' @param save.se.obj Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param boxplot.output Indicates whether to plot the boxplot of the correlation, by default it is set to TRUE.
#' @param plot.top.genes Indicates whether to plot the gene expression of the number of genes
#' from the high or low correlation, by default it is set to FALSE.
#' @param nb.top.genes Defines the number of genes from the high or low correlation to plot,
#' by default is set to 3.
#' @param assess.se.obj Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na TO BE DEFINED.
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#' @param pseudo.count TO BE DEFINED.
#' @param apply.round TO BE DEFINED.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the computed correlation on the continuous variable
#' and if requested the associated plot.
#' @importFrom SummarizedExperiment assays assay colData
#' @importFrom matrixTests row_oneway_equalvar
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom fastDummies dummy_cols
#' @importFrom stats cancor var complete.cases
#' @importFrom stats cor.test p.adjust
#' @importFrom parallel mclapply
#' @importFrom Rfast correls transpose
#' @import ggplot2
#' @export


genesVariableCorrelation<-function(
        se.obj,
        assay.names='All',
        variable,
        apply.log=TRUE,
        method='spearman',
        a = 0.05,
        rho = 0,
        save.se.obj = TRUE,
        boxplot.output=TRUE,
        plot.top.genes = FALSE,
        nb.top.genes = 3,
        assess.se.obj = TRUE,
        remove.na = 'both',
        verbose = TRUE,
        pseudo.count = 1,
        apply.round = TRUE
){
    printColoredMessage(message = '------------The genesVariableCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)
    ### Check the assay names and method input
    if(!method %in% c('pearson', 'spearman')){
        stop('"pearson" and "spearman" are the two supported types for correlations.')
    }
    ### Check se.obj and assay name
    if(assess.se.obj){
        se.obj <- checkSeObj(se.obj = se.obj,
                             assay.names = assay.names,
                             variables = variable,
                             remove.na = remove.na,
                             verbose = verbose)
    }

    ### Check if the variable provided has a variance of 0:
    if(var(se.obj[[variable]]) == 0){
        stop(paste0('The ', variable, ' appears to have no variation.'))
    }

    ## Assays
    if(length(assay.names) == 1 && assay.names=='All'){
        assay.names=as.factor(names(assays(se.obj)))
    }else {
        assay.names=as.factor(unlist(assay.names))
    }

    # Correlation gene expression and continuous variable on all assays
    cor.all<- lapply(
        levels(assay.names),
        function(x){

            ### Log transformation
            printColoredMessage(
                message = '### Data transformation:',
                color = 'magenta',
                verbose = verbose )
            if (apply.log) {
                temp.data <- log2(assay(x = se.obj, i = x) +  pseudo.count)
                printColoredMessage(
                    message = paste0('Performing log2 transformation on the assay ',x),
                    color = 'blue',
                    verbose = verbose )
            } else{
                temp.data <- assay(x = se.obj, i = x)
            }

            ### Correlation between gene expression and continuous variable
            printColoredMessage(message=
                                    paste0(
                                        '### Performing ' ,
                                        method,
                                        ' correlation between individual genes and the ',
                                        variable,
                                        ' variable',
                                        ' in the assay ',
                                        x,
                                        '.'),
                                color = 'magenta',
                                verbose = verbose
            )
            printColoredMessage(message=
                                    paste0(
                                        'Applying the correls function from the Rfast R package. The method = ' ,
                                        method,
                                        ', a = ',
                                        a,
                                        ' and ',
                                        'rho = ',
                                        rho,
                                        ' in the assay ',
                                        x,
                                        '.'),
                                color = 'blue',
                                verbose = verbose
            )
            ## Compute correlation
            corr.genes.var <- correls(
                y = se.obj@colData[, variable],
                x = transpose(temp.data),
                type = method,
                a = a ,
                rho = rho
            )
            row.names(corr.genes.var) <- row.names(se.obj)
            ## Round the correlation obtained to 2 digits
            if(apply.round){
                corr.genes.var <- cbind(
                    round(x = corr.genes.var[,1:4], digits = 2),
                    corr.genes.var[ , 5, drop = FALSE]
                )
            }
            ### Plot top and bottom ranked genes
            if (plot.top.genes) {
                printColoredMessage(message=paste0(
                    '### Plotting top ' ,
                    nb.top.genes,
                    ' highly correlated genes with the ',
                    variable,
                    ' variable', '.'
                ),
                color = 'magenta',
                verbose = verbose
                )
                temp.corr <- corr.genes.var[order(corr.genes.var[, 'correlation'],
                                                  decreasing = TRUE,
                                                  na.last = TRUE) ,]
                ### high positive correlation
                p.pos <- as.data.frame(t(temp.data[row.names(temp.corr)[c(1:nb.top.genes)],]))
                p.pos[ ,'variable'] <- se.obj@colData[, variable]
                p.pos <- p.pos %>% pivot_longer(-variable, names_to = 'genes', values_to = 'expr')
                p.pos <- ggplot(p.pos, aes(x = variable, y = expr)) +
                    geom_point() +
                    ylab(expression(Log[2]~'gene expression')) +
                    xlab(variable) +
                    facet_wrap(~genes) +
                    ggtitle(paste0(nb.top.genes," Top highly positively correlated genes with ",
                                   variable,'\n in the assay ',x))+
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', size = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 10),
                        axis.text.y = element_text(size = 12),
                        legend.text = element_text(size = 10),
                        legend.title = element_text(size = 14),
                        strip.text.x = element_text(size = 10),
                        plot.title = element_text(size = 16)
                    )
                plot(p.pos)

                ### low negative correlation
                temp.corr <- corr.genes.var[order(corr.genes.var[, 'correlation'],
                                                  decreasing = FALSE,
                                                  na.last = TRUE) ,]
                p.neg <- as.data.frame(t(temp.data[row.names(temp.corr)[c(1:nb.top.genes)],]))
                p.neg[ ,'variable'] <- se.obj@colData[, variable]
                p.neg <- p.neg %>% pivot_longer(-variable, names_to = 'genes', values_to = 'expr')
                p.neg <- ggplot(p.neg, aes(x = variable, y = expr)) +
                    geom_point() +
                    ylab(expression(Log[2]~'gene expression')) +
                    xlab(variable) +
                    facet_wrap(~genes) +
                    ggtitle(paste0(
                        'Top highly negatively correlated genes with the variable ',
                        variable,'\n in the assay ',x)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', size = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 10),
                        axis.text.y = element_text(size = 12),
                        legend.text = element_text(size = 10),
                        legend.title = element_text(size = 14),
                        strip.text.x = element_text(size = 10),
                        plot.title = element_text(size = 16)
                    )
                plot(p.neg)
                rm(temp.data)
                rm(temp.corr)

                # results <- list(
                #     corr.genes.var = corr.genes.var,
                #     p.pos=p.pos,
                #     p.neg=p.neg)
            # } else{
                results <- list(
                    corr.genes.var = corr.genes.var)
            }
            return(results)
        })
    names(cor.all) <- levels(assay.names)


    ### Add results to the SummarizedExperiment object
    if(save.se.obj == TRUE){
        printColoredMessage(message= '### Saving the correlation results to the metadata of the SummarizedExperiment object.',
                        color = 'magenta',
                        verbose = verbose)

        for (x in levels(assay.names)){
            ## Check if metadata metric already exist
            if(length(se.obj@metadata)==0 ) {
                se.obj@metadata[['metric']] <- list()
            }
            ## Check if metadata metric already exist for this assay
            if(!'metric' %in% names(se.obj@metadata) ) {
                se.obj@metadata[['metric']] <- list()
            }
            ## Check if metadata metric already exist for this assay
            if(!x %in% names(se.obj@metadata[['metric']]) ) {
                se.obj@metadata[['metric']][[x]] <- list()
            }
            ## Check if metadata metric already exist for this assay and this metric
            if(!paste0('gene.',method,'.corr') %in% names(se.obj@metadata[['metric']][[x]])  ) {
                se.obj@metadata[['metric']][[x]][[paste0('gene.',method,'.corr')]] <- list()
            }
            ## Check if metadata metric already exist for this assay, this metric and this variable
            if(! variable %in% names(se.obj@metadata[['metric']][[x]][[paste0('gene.',method,'.corr')]])  ) {
                se.obj@metadata[['metric']][[x]][[paste0('gene.',method,'.corr')]][[variable]] <- cor.all[[x]][['corr.genes.var']][,'correlation']
            }
        }

        printColoredMessage(message= paste0(
            'The correlation results are saved to metadata@',
            x,
            '$gene.var.corr$',
            variable,
            '.'),
            color = 'blue',
            verbose = verbose)

        ## Plot and save the plot into se.obj@metadata$plot
        if (boxplot.output==TRUE) {
            printColoredMessage(message= '### Plotting and Saving the correlation results to the metadata of the SummarizedExperiment object.',
                                color = 'magenta',
                                verbose = verbose)

            se.obj=plotMetric(se.obj,
                              assay.names =assay.names,
                              metric=paste0('gene.',method,'.corr'),
                              variable=variable,
                              verbose=verbose)
        }
        return(se.obj = se.obj)

    ## Return only the correlation result
    } else if(save.se.obj == FALSE){
        return(gene.corr.var=cor.all)
    }

    printColoredMessage(message = '------------The genesVariableCorrelation function finished.',
                        color = 'white',
                        verbose = verbose)

}



