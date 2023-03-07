#' is used to compute the Spearman correlation between the data and library size
#' on all assay
#'
#' @param sce the dataset that will be used for this analysis
#' @param corr_var The continuous variable that will be used to compute to correlation
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
        corr_var,
        apply.log=FALSE,
        n.cores=5
){
    normalization=names(assays(sce))
    # Correlation gene expression and lib size
    cor.ls <- lapply(
        normalization,
        function(x){
            data <- as.matrix(assay(sce, x))
            cor <- correlation_gene_exp_contvar_single_assay(
                expr.data = data,
                apply.log=apply.log,
                variable=corr_var,
                method='spearman',
                n.cores = n.cores)
            cor
        })
    names(cor.ls) <- normalization
    cor.ls.coeff <- lapply(
        normalization,
        function(x){
            cor.ls[[x]]$rho
        })
    names(cor.ls.coeff) <- normalization
    everything<-datasets<-corr.coeff<-NULL
    cor.ls.coeff <- as.data.frame(cor.ls.coeff)
    cor.ls.coeff= cor.ls.coeff %>% pivot_longer(
        everything(),
        names_to = 'datasets',
        values_to = 'corr.coeff') %>% mutate(datasets = factor(
            datasets))
    ### Plot
    # color
    dataSets.colors <- wes_palette(
        n = 4,
        name = "GrandBudapest1")[c(1,2,4,3)]
    p=ggplot(cor.ls.coeff, aes(x = datasets, y = corr.coeff, fill = datasets)) +
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

    return(list(plot=p,corr.ls.coeff=cor.ls.coeff))
}

