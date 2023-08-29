#' is used to compute the differential gene expression analysis between the gene expression (assay)
#' of a SummarizedExperiment class object and a categorical variable (i.e. batches) using ANOVA.
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the differential expression ANOVA analysis.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the correlation. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' the log transformation will be selected.
#' @param method A character string indicating which method
#' is to be used for the differential analysis: "aov" or "welch.correction". By default "aov" will
#' be selected.
#' @param save.se.obj Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param boxplot.output Indicates whether to plot the boxplot of the F-test statistics, by default it is set to TRUE.
#' @param plot.top.genes Indicates whether to plot the gene expression of the number of genes
#' from the top listing of anova by default it is set to FALSE.
#' @param nb.top.genes Defines the number of genes from the top or bottom listing of anova to plot,
#' by default is set to 3.
#' @param assess.se.obj Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na TO BE DEFINED.
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#' @param pseudo.count TO BE DEFINED.
#' @param apply.round TO BE DEFINED.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the associated plot or the computed correlation on the continuous variable.

#' @importFrom SummarizedExperiment assay assays
#' @importFrom matrixTests row_oneway_equalvar row_oneway_welch
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom stats anova
#' @import ggplot2
#'
#' @export
genesVariableAnova <- function(se.obj,
                               assay.names='All',
                               variable,
                               apply.log = TRUE,
                               method = 'aov',
                               save.se.obj = TRUE,
                               boxplot.output=TRUE,
                               plot.top.genes = FALSE,
                               nb.top.genes = 3,
                               assess.se.obj = TRUE,
                               remove.na = 'both',
                               verbose = verbose,
                               pseudo.count = 1,
                               apply.round = TRUE
){
    printColoredMessage(message = '------------The genesVariableAnova function starts:',
                        color = 'white',
                        verbose = verbose)
    ### Check the inputs
    if (length(unique(se.obj@colData[, variable])) < 2) {
        stop(paste0('The ', variable,', contains only one variable.'))
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0('The ', variable,', is a numeric, but this should a categorical variable'))
    } else if(!method %in% c('aov', 'welch.correction') ){
        stop('"aov" and "welch.correction" are the two supported types for correlations.')
    }
    ### Assess the
    if(assess.se.obj){
        se.obj <- checkSeObj(se.obj = se.obj,
                             assay.names = assay.names,
                             variables = variable,
                             remove.na = remove.na,
                             verbose = verbose)
    }

    ## Assays
    if(length(assay.names) == 1 && assay.names=='All'){
        assay.names=as.factor(names(assays(se.obj)))
    }else {
        assay.names=as.factor(unlist(assay.names))
    }

    # ANOVA on gene expression and categorical variable on all assays
    anova.all<- lapply(
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
                    message = 'Performing log2 transformation on the data.',
                    color = 'blue',
                    verbose = verbose )
            } else{
                temp.data <- assay(x = se.obj, i = x)
            }

            ### ANOVA on gene expression and categorical variable
            printColoredMessage(message = paste0(
                '### Performing ANOVA between individual genes and the ',
                variable,
                ' variable.'
            ),
            color = 'magenta',
            verbose = verbose)

            if(method == 'aov'){
                printColoredMessage(message = paste0(
                    'Applying the row_oneway_equalvar function from the Rfast R package between individual genes and the ',
                    variable,
                    ' variable.',
                    ' in the assay ',
                    x
                ),
                color = 'blue',
                verbose = verbose)
                anova.genes.var <- row_oneway_equalvar(
                    x = temp.data,
                    g = se.obj@colData[, variable]
                )
                row.names(anova.genes.var) <- row.names(se.obj)
            } else if(method == 'welch.correction'){
                printColoredMessage(message = paste0(
                    'Applying the row_oneway_welch function from the Rfast R package between individual genes and the ',
                    variable,
                    ' variable.',
                    ' in the assay ',
                    x
                ),
                color = 'blue',
                verbose = verbose)
                anova.genes.var <- row_oneway_welch(
                    x = temp.data,
                    g = se.obj@colData[, variable]
                )
                row.names(anova.genes.var) <- row.names(se.obj)
            }

            ## Round the anova statistic obtained to 2 digits
            if(apply.round){
                anova.genes.var <- cbind(
                    round(anova.genes.var[,1:9], digits = 3),
                    anova.genes.var[ , 10, drop = FALSE]
                )
            }
            if (plot.top.genes) {
                temp.anova <- anova.genes.var[ order(anova.genes.var[, 'statistic'],
                                                     decreasing = TRUE,
                                                     na.last = TRUE) , ]
                ### positive correlation
                var<-NULL
                p.high <- as.data.frame(t(temp.data[row.names(temp.anova)[c(1:nb.top.genes)],]))
                p.high <- mutate(p.high , var = se.obj@colData[, variable])
                p.high <- pivot_longer(
                    data = p.high,
                    cols =  -var,
                    names_to = 'genes',
                    values_to = 'expr'
                )
                p.high <- ggplot(p.high, aes(x = var, y = expr)) +
                    geom_boxplot() +
                    ylab(expression(Log[2]~'gene expression')) +
                    xlab(variable) +
                    facet_wrap( ~ genes) +
                    ggtitle(paste0(nb.top.genes," Top affected genes by the variable ", variable, " for ", x)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
                        axis.text.y = element_text(size = 12),
                        legend.text = element_text(size = 10),
                        legend.title = element_text(size = 14),
                        strip.text.x = element_text(size = 10),
                        plot.title = element_text(size = 12)
                    )
                plot(p.high)
                rm(temp.data)
                rm(temp.anova)
                results <- list(
                    anova.genes.var =  anova.genes.var)
            }
            return(results)
        })
    names(anova.all) <- levels(assay.names)

    ### Add results to the SummarizedExperiment object
    if(save.se.obj == TRUE){
        printColoredMessage(message= '### Saving the ANOVA results to the metadata of the SummarizedExperiment object.',
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
            if(!paste0('gene.',method,'.anova') %in% names(se.obj@metadata[['metric']][[x]])  ) {
                se.obj@metadata[['metric']][[x]][[paste0('gene.',method,'.anova')]] <- list()
            }
            ## Check if metadata metric already exist for this assay, this metric and this variable
            if(!variable %in% names(se.obj@metadata[['metric']][[x]][[paste0('gene.',method,'.anova')]])  ) {
                se.obj@metadata[['metric']][[x]][[paste0('gene.',method,'.anova')]][[variable]] <- anova.all[[x]][['anova.genes.var']][,'statistic']
            }
        }
            printColoredMessage(message= paste0(
                'The anova results are saved to metadata@',
                x,
                '$gene.var.anova$',
                variable, '.'),
                color = 'blue',
                verbose = verbose)

            ## Plot and save the plot into se.obj@metadata$plot
            if (boxplot.output==TRUE) {
                printColoredMessage(message= '### Plotting and Saving the anova results to the metadata of the SummarizedExperiment object.',
                                    color = 'magenta',
                                    verbose = verbose)

                se.obj=plotMetric(se.obj,
                                  assay.names =assay.names,
                                  metric=paste0('gene.',method,'.anova'),
                                  variable=variable,
                                  verbose=verbose)
            }
            return(se.obj = se.obj)

            ## Return only the correlation result
        } else if(save.se.obj == FALSE){
            return(gene.anova.var=anova.all[[x]][['anova.genes.var']][,'statistic'])
        }

            printColoredMessage(message = '------------The genesVariableAnova function finished.',
                                color = 'white',
                                verbose = verbose)

}
