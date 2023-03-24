#' is used to compute the correlation between the gene expression (assay)
#' of a SummarizedExperiment class object and a continous variable (i.e. library size).
#'
#' @param se A SummarizedExperiment object that will be used to compute the correlation
#' @param cat_var Vector of a categorical variable such as sample types
#' (i.e. biological subtypes) or batches.
#' @param cat_var_label String or vector of strings of the label of categorical variable(s) such as
#' sample types or batches from colData(se).
#' @param assay_names Optional string or list of strings for selection of the names
#' of the assays of the SummarizedExperiment class object to compute the correlation.
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#' @param plot_output Indicates whether to plot the results
#' @param top.genes.no Defines the number of top genes, by default is set to 3.
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom SummarizedExperiment assays assay
#' @importFrom matrixTests row_oneway_equalvar
#' @import ggplot2
#' @export

anova_gene_exp_contvar<-function(
        se,
        cat_var,
        cat_var_label,
        assay_names=NULL,
        apply.log=FALSE,
        plot_output = TRUE,
        top.genes.no = 3
){
    ### check the inputs
    if( !identical(cat_var,se@colData[, cat_var_label])){
        stop(paste0(
            'There label of the categorical variable ',
            cat_var_label,
            'is different from the categorical variable provided ',
            cat_var,
            'please provide the corresponding label and categorical variable.\n'))
    }

    if( anyNA(se@colData[, cat_var_label]) ){
        stop(paste0(
            'There is/are NA in the ',
            cat_var_label,
            ', please remove them and re-run the function.\n'))
    }
    if( length(unique(se@colData[, cat_var_label])) < 2 ){
        stop(paste0(
            'The ',
            cat_var,
            ', contains only one variable. Please provide at least 2.\n'))
    }
    if( is.numeric(class(se@colData[, cat_var_label]))){
        stop(paste0(
            'The ',
            cat_var,
            ', is numeric, please provide a categorical variable.\n'))
    }
    ## Assays
    if (!is.null(assay_names)){
        normalization=assay_names
    }else{
        normalization=names(assays(se))
    }

    ### ANOVA
    message(paste0(
        'Performing ANOVA between individual genes and the ',
        cat_var_label, ' variable \n')
    )

    anova.all<- lapply(
        normalization,
        function(x){
            ### log transformation
            if(apply.log){
                #message('Performing log + 1 transformation on the data')
                y <- log2(assay(x = se, x) + 1)
            }
            else{
                y <- assay(x = se, x)
            }

            ## Oneway ANOVA on row
            anova.batch.genes <- row_oneway_equalvar(
                x = y,
                g = se@colData[, cat_var_label])

            if(plot_output){
                anova.batch.genes <- anova.batch.genes[
                    order(
                        anova.batch.genes[ , 'statistic'],
                        decreasing = TRUE,
                        na.last = TRUE) , ]
                ### positive correlation
                var<-NULL
                p.high <- as.data.frame(t(y[row.names(anova.batch.genes)[c(1:top.genes.no)], ]))
                p.high <- mutate(p.high , var = se@colData[ , cat_var_label])
                p.high <- pivot_longer(data = p.high, cols =  -var, names_to = 'genes', values_to = 'expr')
                p.high <- ggplot(p.high, aes(x = var, y = expr)) +
                    geom_boxplot() +
                    ylab('Gene expression (log)') +
                    xlab(cat_var_label) +
                    facet_wrap(~genes) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', size = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 10, angle = 45),
                        axis.text.y = element_text(size = 12),
                        legend.text = element_text(size = 10),
                        legend.title = element_text(size = 14),
                        strip.text.x = element_text(size = 10),
                        plot.title = element_text(size = 16)+
                            ggtitle(paste("ANOVA pos. correlation computed on ",cat_var_label,sep=""))
                    )
                ### negative correlation
                p.low <- as.data.frame(t(y[row.names(anova.batch.genes)[c(c(nrow(anova.batch.genes)-c(top.genes.no -1)): nrow(anova.batch.genes))], ]))
                p.low <- mutate(p.low , var = se@colData[ , cat_var_label])
                p.low <- tidyr::pivot_longer(data = p.low, cols = -var, names_to = 'genes', values_to = 'expr')
                p.low <- ggplot(p.low, aes(x = var, y = expr)) +
                    geom_boxplot() +
                    ylab('Gene expression (log)') +
                    xlab(cat_var_label) +
                    facet_wrap(~genes) +
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
                        plot.title = element_text(size = 16)+
                            ggtitle(paste("ANOVA neg. correlationcomputed on ",cat_var_label,sep=""))
                    )
                return(list(
                    anova.batch.genes =anova.batch.genes[row.names(se), ],
                    plot.low = p.low,
                    plot.high = p.high)
                    )

            }
            else{
                return(list(anova.batch.genes = anova.batch.genes[row.names(se), ])#,
                       #plot.low = NULL,
                       #plot.high = NULL)
                       )
            }
    })
    names(anova.all)=normalization
    return(anova.all)
}
