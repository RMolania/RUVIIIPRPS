#' is used to compute the differential gene expression analysis between the gene expression (assay)
#' of a SummarizedExperiment class object and a categorical variable (i.e. batches) using ANOVA.
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
#' @param boxplot_output Indicates whether to plot the boxplot of the F-values of anova.
#' @param ranked_genes_plot_output Indicates whether to plot the gene expression of the number of genes
#' from the top or bottom listing of anova
#' @param nb_ranked_genes Defines the number of genes from the top or bottom listing of anova to plot,
#' by default is set to 3.
#'
#'@return list List containing the associated plot and the anova computed for the categorical variable.

#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom SummarizedExperiment assays assay
#' @importFrom matrixTests row_oneway_equalvar
#' @importFrom kunstomverse geom_boxplot2
#' @importFrom wesanderson wes_palette wes_palette
#' @importFrom stats anova
#' @import ggplot2
#' @export

anova_gene_exp_catvar<-function(
        se,
        cat_var,
        cat_var_label,
        assay_names=NULL,
        apply.log=FALSE,
        boxplot_output = TRUE,
        ranked_genes_plot_output = FALSE,
        nb_ranked_genes = 3
){
    ### Check se and assay names
    if (!class(se)[1] == 'SummarizedExperiment') {
        stop('Please provide a summarized experiment object.\n')
    } else if((!is.null(assay_names))&& (any(assay_names %in% names(assays(se)))=='FALSE')){
        stop('The selected assay is/are not in the assay names of the SummarizedExperiment class object.\n')
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
            cat_var_label,
            ', contains a unique variable. Please provide at least 2.\n'))
    }
    if( is.numeric(class(se@colData[, cat_var_label]))){
        stop(paste0(
            'The ',
            cat_var_label,
            ', is numeric, please provide a categorical variable.\n'))
    }
    ## Assays
    if (!is.null(assay_names)){
        normalization=as.factor(assay_names)
    }else{
        normalization=as.factor(names(assays(se)))
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
                expr <- log2(assay(x = se, x) + 1)
            }
            else{
                expr <- assay(x = se, x)
            }

            ## Oneway ANOVA on row
            anova<- row_oneway_equalvar(
                x = expr,
                g = se@colData[, cat_var_label])

            ## Top and bottom ranked genes
            if(ranked_genes_plot_output){
                ## Bottom ranked genes
                anova_sorted <-anova[order(anova[ , 'statistic'],
                                           decreasing = TRUE,
                                           na.last = TRUE) , ]
                var<-NULL
                p.top <- as.data.frame(t(expr[row.names(anova_sorted)[c(1:nb_ranked_genes)], ]))
                p.top <- mutate(p.top , var = se@colData[ , cat_var_label])
                p.top <- pivot_longer(data = p.top, cols =  -var, names_to = 'genes', values_to = 'expr')
                p.top=ggplot(p.top, aes(x = var, y = expr)) +
                    geom_boxplot() +
                    ylab('Gene expression (log)') +
                    xlab(cat_var_label) +
                    facet_wrap(~genes) +
                    ggtitle(paste(nb_ranked_genes," Top ranked genes from ANOVA computed on ",
                                  cat_var_label," for ",x,sep=""))+
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
                        plot.title = element_text(size = 16)
                    )
                ## Bottom ranked genes
                var<-NULL
                ### Bottom ranked genes
                p.bottom <- as.data.frame(t(expr[row.names(anova_sorted)[c(c(nrow(anova_sorted)-c(nb_ranked_genes -1)): nrow(anova))], ]))
                p.bottom <- mutate(p.bottom , var = se@colData[ , cat_var_label])
                p.bottom <- tidyr::pivot_longer(data = p.bottom, cols = -var, names_to = 'genes', values_to = 'expr')
                p.bottom<- ggplot(p.bottom, aes(x = var, y = expr)) +
                    geom_boxplot() +
                    ylab('Gene expression (log)') +
                    xlab(cat_var_label) +
                    facet_wrap(~genes) +
                    ggtitle(paste(nb_ranked_genes," Bottom ranked genes from ANOVA computed on ",
                                  cat_var_label," for ",x,sep=""))+
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
                results <- list(
                    anova = anova,
                    plot.top=p.top,
                    plot.bottom= p.bottom)
            } else{
                results <- list(
                    anova = anova)
            }

            return(results)

        })
    names(anova.all) <- normalization

        ### Boxplot of the F-test association between the variable and gene expression
        if(boxplot_output){
            ftest.all <- lapply(
                normalization,
                function(x){
                    as.numeric(anova.all[[x]]$anova$statistic)
                })
            names(ftest.all) <- normalization
            datasets<-fval<-everything<-NULL
            ## length of assays
            assays_nb=length(normalization)
            ftest.all = as.data.frame(ftest.all) %>% pivot_longer( everything(),
                        names_to = 'datasets',values_to = 'fval') %>%
                        mutate(datasets = factor(datasets,levels=normalization))
            # color
            dataSets.colors <- wes_palette(
                n = length(normalization),
                name = "GrandBudapest1")[seq(1:length(normalization))]
            boxplot_p=ggplot(ftest.all, aes(x = datasets, y = log2(fval),fill = datasets)) +
                geom_boxplot2() +
                ylab(expression(Log[2]~'F statistics')) +
                xlab('') +
                scale_fill_manual(values = dataSets.colors, guide = 'none') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black'),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 8),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14))+
            ggtitle(paste(" F-test distribution from ANOVA between the gene expression and ",cat_var_label,sep=""))
            }

        if(boxplot_output==T ){
            return(list(anova =anova.all,
                    boxplot_ftest=boxplot_p))
        } else {
            return(list(anova = anova.all))
        }

}
