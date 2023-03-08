#' is used to compute the Spearman correlation between the gene expression of all assays
#' and a continous variable (i.e. library size)
#'
#' @param sce the dataset that will be used for this analysis
#' @param cont_var The continuous variable that will be used to compute to correlation
#' @param apply.log Indicates whether to apply a log-transformation to the data
#' @param n.cores is the number of cpus used for mclapply parallelization
#'
#' @return list List containing the association plot and the computed correlation
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr rename mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom SummarizedExperiment assay
#' @import ggplot2
#' @export

correlation_gene_exp_contvar_all_assays<-function(
        sce,
        cont_var,
        apply.log=FALSE,
        n.cores=5
){
    normalization=names(assays(sce))
    # Correlation gene expression and continous variable
    cor.all<- lapply(
        normalization,
        function(x){
            data <- as.matrix(assay(sce, x))
            cor <- correlation_gene_exp_contvar_single_assay(
                expr.data = data,
                apply.log=apply.log,
                cont_var=cont_var,
                method='spearman',
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
        n = 4,
        name = "GrandBudapest1")[c(1,2,4,3)]
    p=ggplot(cor.all.coeff, aes(x = datasets, y = corr.coeff, fill = datasets)) +
        geom_boxplot() +
        ylab("Spearman correlation") +
        xlab('') +
        geom_hline(yintercept=0)+
        scale_fill_manual(values = dataSets.colors, guide = 'none') +
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', size = 1),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12))

    return(list(plot=p,corr.coeff=cor.all.coeff))
}

