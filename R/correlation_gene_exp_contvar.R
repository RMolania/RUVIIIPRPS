#' is used to compute the correlation between the gene expression (assay)
#' of a SummarizedExperiment class object and a continuous variable (i.e. library size).
#'
#' @param se A SummarizedExperiment object that will be used to compute the correlation
#' @param cont_var Vector of a continous variable such as sample types
#' (i.e. biological subtypes) or batches.
#' @param cont_var_label String or vector of strings of the label of continous variable(s) such as
#' sample types or batches from colData(se).
#' @param method A character string indicating which correlation coefficient
#' is to be used for the test: "pearson", "kendall", or "spearman". By default 'spearman will
#' be selected.
#' @param assay_names Optional string or list of strings for selection of the names
#' of the assays of the SummarizedExperiment class object to compute the correlation.
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#' @param boxplot_output Indicates whether to plot the boxplot of the correlation.
#' @param ranked_genes_plot_output Indicates whether to plot the gene expression of the number of genes
#' from the high or low correlation.
#' @param a The significance level used for the confidence intervals in the correlation,
#' by default it is set to 0.05.
#' @param rho The value of the hypothesised correlation to be used in the hypothesis testing,
#' by default it is set to 0.
#' @param nb_ranked_genes Defines the number of genes from the top or bottom listing of anova to plot,
#' by default is set to 3.
#'
#' @return list List containing the associated plot and the computed correlation on the continuous variable.
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr rename mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom SummarizedExperiment assays assay
#' @importFrom stats cor.test p.adjust
#' @importFrom parallel mclapply
#' @importFrom Rfast correls transpose
#' @import ggplot2
#' @export

correlation_gene_exp_contvar<-function(
        se,
        cont_var,
        cont_var_label,
        assay_names=NULL,
        method='spearman',
        apply.log=FALSE,
        boxplot_output = TRUE,
        ranked_genes_plot_output = FALSE,
        a = 0.05,
        rho = 0,
        nb_ranked_genes = 3
){
    ### check the inputs
    # if( !identical(cont_var,as.vector(se@colData[, cont_var_label]))){
    #     stop(paste0(
    #         'The label of the continous variable ',
    #         cont_var_label,
    #         ' is different from the continous variable provided,
    #         please provide the corresponding label and continous variable.\n'))
    # }
    if (!is.null(assay_names)){
        normalization=assay_names
    }else{
        normalization=names(assays(se))
    }
    # Correlation gene expression and continous variable
    cor.all<- lapply(
        normalization,
        function(x){
            correlation_gene_exp_contvar_single_assay <- function(expr.data,
                                                 apply.log,
                                                 cont_var_label,
                                                 method,
                                                 a,
                                                 rho){
                                                 #n.cores) {
                if(apply.log==FALSE){
                    expr.data <- expr.data
                }
                else{
                    expr.data <- log2(expr.data + 1)
                }
                corr.genes.var <- correls(
                    y = se@colData[ , cont_var_label],
                    x = transpose(expr.data),
                    type = method,
                    a = a ,
                    rho = rho
                )
                corr.genes.var

                # rho <- mclapply(
                #     1:nrow(expr.data),
                #     function(x){
                #         round(cor.test(
                #             x = expr.data[x, ],
                #             y = cont_var,
                #             method = method)[[4]], 6)},
                #     mc.cores = n.cores
                # )
                # pval <- mclapply(
                #     1:nrow(expr.data),
                #     function(x){
                #         cor.test(
                #             x = expr.data[x, ],
                #             y = cont_var,
                #             method = method)[[3]]},
                #     mc.cores = n.cores)
                #
                # results <- data.frame(
                #     genes = row.names(expr.data),
                #     rho = unlist(rho),
                #     pvalue = unlist(pval),
                #     adj.pvalue = p.adjust(unlist(pval), 'BH')
                # )
                # return(results)
            }


            data <- as.matrix(assay(se, x))
            cor <- correlation_gene_exp_contvar_single_assay(
                expr.data = data,
                apply.log=apply.log,
                cont_var_label=cont_var_label,
                method=method,
                a,
                rho)
                #n.cores = n.cores)
            row.names(cor) <- row.names(se)
            cor

            ## Top and bottom ranked genes
            if(ranked_genes_plot_output){
                cor.sorted <- cor[order(cor[ , 'correlation'],
                        decreasing = TRUE,
                        na.last = TRUE) , ]
                ### high positive correlation scatterplot
                p.high <- as.data.frame(t(data[row.names(cor.sorted)[c(1:nb_ranked_genes)], ]))
                p.high[ , 'cont_var'] <- se@colData[ , cont_var_label]
                p.high <- p.high %>% pivot_longer(-cont_var, names_to = 'genes', values_to = 'expr')
                p.high <- ggplot(p.high, aes(x = cont_var, y = expr)) +
                    geom_point() +
                    ylab('Gene expression (log)') +
                    xlab(cont_var_label) +
                    facet_wrap(~genes) +
                    ggtitle(paste(nb_ranked_genes," Top ranked genes from correlation computed on ",
                                  cont_var_label," for ",x,sep=""))+
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

                ### low negative correlation
                p.low <- as.data.frame(t(data[row.names(cor.sorted)[c(c(nrow(cor.sorted)-c(nb_ranked_genes -1)): nrow(cor.sorted))], ]))
                p.low[ , 'cont_var'] <- se@colData[ , cont_var_label]
                p.low <- p.low %>% pivot_longer(-cont_var, names_to = 'genes', values_to = 'expr')
                p.low <- ggplot(p.low, aes(x = cont_var, y = expr)) +
                    geom_point() +
                    ylab('Gene expression (log)') +
                    xlab(cont_var_label) +
                    facet_wrap(~genes) +
                    ggtitle(paste(nb_ranked_genes," Bottom ranked genes from correlation computed on ",
                                  cont_var_label," for ",x,sep=""))+
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
                    cor =  cor,
                    plot.high=p.high,
                    plot.low= p.low)
            } else{
                results <- list(
                    cor = cor)

            }
        })
    names(cor.all) <- normalization
    cor.all.coeff <- lapply(
        normalization,
        function(x){
            cor.all[[x]]$cor[,'correlation']
        })
    names(cor.all.coeff) <- normalization
    everything<-datasets<-corr.coeff<-NULL
    cor.all.coeff <- as.data.frame(cor.all.coeff)
    cor.all.coeff= cor.all.coeff %>% pivot_longer(
        everything(),
        names_to = 'datasets',
        values_to = 'corr.coeff') %>% mutate(datasets = factor(
            datasets))
    ### Boxplot of the F-test association between the variable and gene expression
    if(boxplot_output){
        dataSets.colors <- wes_palette(
            n = length(normalization),
            name = "GrandBudapest1")[c(1,2,4,3)]
        p=ggplot(cor.all.coeff, aes(x = datasets, y = corr.coeff, fill = datasets)) +
            geom_boxplot() +
            ylab(paste(method,"correlation",sep=" ")) +
            xlab('') +
            geom_hline(yintercept=0)+
            scale_fill_manual(values = dataSets.colors, guide = 'none') +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', size = 1),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12))+
            ggtitle(paste("Correlation between the gene expression and ",cont_var_label,sep=""))
    }
    if(boxplot_output==T ){
        return(list(cor=cor.all,
                    plot=p,
                    cont_var_label=cont_var_label))
    } else {
        return(list(cor=cor.all,
                    cont_var_label=cont_var_label))
    }
}

