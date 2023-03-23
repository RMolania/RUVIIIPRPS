#' is used to compute the correlation between the gene expression (assay)
#' of a SummarizedExperiment class object and a continous variable (i.e. library size).
#'
#' @param se A SummarizedExperiment object that will be used to compute the correlation
#' @param cont_var Vector of a categorical variable such as sample types
#' (i.e. biological subtypes) or batches.
#' @param cont_var_label String or vector of strings of the label of categorical variable(s) such as
#' sample types or batches from colData(se).
#' @param method A character string indicating which correlation coefficient
#' is to be used for the test: "pearson", "kendall", or "spearman". By default 'spearman will
#' be selected.
#' @param assay_names Optional string or list of strings for selection of the names
#' of the assays of the SummarizedExperiment class object to compute the correlation.
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#' @param n.cores is the number of cpus used for mclapply parallelization.
#'
#' @return list List containing the association plot and the computed correlation
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr rename mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom SummarizedExperiment assays assay
#' @importFrom stats cor.test p.adjust
#' @importFrom parallel mclapply
#' @import ggplot2
#' @export

correlation_gene_exp_contvar<-function(
        se,
        cont_var,
        cont_var_label,
        assay_names=NULL,
        method='spearman',
        apply.log=FALSE,
        n.cores=5
){
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
                                                 cont_var,
                                                 method,
                                                 n.cores) {
                if(apply.log==FALSE){
                    expr.data <- expr.data
                }
                else{
                    expr.data <- log2(expr.data + 1)
                }
                rho <- mclapply(
                    1:nrow(expr.data),
                    function(x){
                        round(cor.test(
                            x = expr.data[x, ],
                            y = cont_var,
                            method = method)[[4]], 6)},
                    mc.cores = n.cores
                )
                pval <- mclapply(
                    1:nrow(expr.data),
                    function(x){
                        cor.test(
                            x = expr.data[x, ],
                            y = cont_var,
                            method = method)[[3]]},
                    mc.cores = n.cores)

                results <- data.frame(
                    genes = row.names(expr.data),
                    rho = unlist(rho),
                    pvalue = unlist(pval),
                    adj.pvalue = p.adjust(unlist(pval), 'BH')
                )
                return(results)
            }

            data <- as.matrix(assay(se, x))
            cor <- correlation_gene_exp_contvar_single_assay(
                expr.data = data,
                apply.log=apply.log,
                cont_var=cont_var,
                method=method,
                n.cores = n.cores)
            cor
        })
    names(cor.all) <- normalization
    cor.all.coeff <- lapply(
        normalization,
        function(x){
            cor.all[[x]]$rho
        })
    names(cor.all.coeff) <- normalization
    everything<-datasets<-corr.coeff<-NULL
    cor.all.coeff <- as.data.frame(cor.all.coeff)
    cor.all.coeff= cor.all.coeff %>% pivot_longer(
        everything(),
        names_to = 'datasets',
        values_to = 'corr.coeff') %>% mutate(datasets = factor(
            datasets))
    ### Plot
    # color
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


    return(list(plot=p,corr.coeff=cor.all.coeff,cont_var_label=cont_var_label))
}

