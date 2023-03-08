
#' is used to compute the differential gene expression analysis on all the assays
#' given a categorical variable of two groups of the samples using the Wilcoxon test
#'
#' @param sce the dataset that will be used for this analysis
#' @param cat_var The categorical variable of two groups that will be used
#' to computed the differential expression analysis (i.e. high vs low library size)
#' @param apply.log Indicates whether to apply a log-transformation to the data
#' @param n.cores is the number of cpus used for mclapply parallelization
#'
#' @return list List containing the association plot and the computed regression
#' @importFrom stats lm
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr rename mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom SummarizedExperiment assay
#' @import ggplot2
#' @export

da_analysis_wilcoxon_gene_exp_catvar_all_assays<-function(
        sce,
        cat_var,
        apply.log=FALSE,
        n.cores=5
){
normalization=names(assays(sce))
# Wilcoxon test
de <- lapply(
    normalization,
    function(x){
        data <- assay(sce, x)
        de <- da_analysis_wilcoxon_gene_exp_catvar_single_assay(
            expr.data = data,
            apply.log,
            cat_var,
            n.cores = n.cores)
        de
    })
names(de) <- normalization
pval.de <- lapply(
    normalization,
    function(x){
        de[[x]]$pvalue
    })
names(pval.de) <- normalization
everything<-datasets<-p.val<-NULL
pval.de <- as.data.frame(pval.de)
pval.de= pval.de %>% pivot_longer(
        everything(),
        names_to = 'datasets',
        values_to = 'p.val') %>% mutate(datasets = factor(
        datasets))
### Plot
p=ggplot(pval.de, aes(p.val)) +
    geom_histogram(binwidth = .1) +
    scale_x_continuous(breaks = c(seq(0, 1, .5))) +
    xlab('p_values') + ylab('Frequency') +
    facet_wrap( ~ datasets, ncol = 4) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black', size = 1),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 16))

return(list(plot=p,pval.de=pval.de))
}

