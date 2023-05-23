#' is used to compute the adjusted rank index (ARI) from the first PC of a SummarizedExperiment class
#' object given a categorical variable.
#'
#' It can be used to assess how a group of biological samples
#' are distributed across batches (i.e. example subtypes vs batch), and how batches ### update
#'
#'
#' @param pca PCA components of a SummarizedExperiment variable.
#' @param cat_var Vector of a categorical variable such as sample types
#' (i.e. biological subtypes) or batches.
#' @param cat_var_label String or vector of strings of the label of categorical variable(s) such as
#' sample types or batches from colData(se).
#' @param assay_names Optional string or list of strings for selection of the names
#' of the assays of the SummarizedExperiment class object to compute the ARI.
#' @param plot Optional output of a plot, default set to FALSE.
#' @param nPCs is the number of PCs used to measure the distance, default is set to 3.
#'
#' @return list List containing the association plot and the computed ari and the categorical variable.
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr rename mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom stats dist
#' @importFrom mclust mclustBIC Mclust adjustedRandIndex
#' @import ggplot2
#' @export
#'

compute_ari <-function(
        pca,
        cat_var,
        cat_var_label,
        assay_names=NULL,
        plot=FALSE,
        nPCs=3
){
    if (!is.null(assay_names)){
        normalization=as.factor(assay_names)
    }else{
        normalization=as.factor(names(pca))
    }
    # ARI on all assays
    ari <- lapply(
        normalization,
        function(x){
            ari_catvar_single_assay <- function(pca,
                                                cat_var,
                                                nPCs){
                BIC <- mclustBIC(data = pca)
                mod <- Mclust(data = pca, x = BIC)
                ari=adjustedRandIndex(
                    mod$classification,
                    cat_var)
                return(ari)
            }
            coef=ari_catvar_single_assay(pca[[x]]$sing.val$u,
                                        cat_var,
                                        nPCs)
            coef
        })
    names(ari) <- normalization
    everything<-datasets<-silh.coeff<-NULL
    pcs.ari <- as.data.frame(ari)
    pcs.ari =  pcs.ari %>% pivot_longer(
        everything(),
        names_to = 'datasets',
        values_to = 'ari') %>% mutate(datasets = factor(
            datasets),levels=normalization)

    ### Plot
    if (isTRUE(plot)){
    dataSets.colors <- wes_palette(
        n = length(normalization),
        name = "GrandBudapest1")[seq(1:length(normalization))]
    p=ggplot(pcs.ari , aes(x = datasets, y = ari, fill = datasets)) +
        geom_point() +
        ylab("ARI") +
        xlab('') +
        scale_fill_manual(values = dataSets.colors, guide = 'none')+
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', size = 1),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12))+
        ggtitle(paste("ARI computed on ",cat_var_label,sep=""))
    return(list(plot=p,ari=pcs.ari,cat_var_label=cat_var_label))
    }else{
        return(list(ari=pcs.ari,cat_var_label=cat_var_label))}
}

