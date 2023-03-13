#' is used to compute the adjusted rank index (ARI) from the first PC of a SummarizedExperiment class
#' object given a categorical variable.
#'
#'
#' @param pca PCA components of a SummarizedExperiment variable.
#' @param cat_var Vector of a categorical variable such as sample types or batches.
#' @param assay_names Optional string or list of strings for selection of the names
#' of the assays of the SummarizedExperiment class object to compute the ARI.
#' @param plot Optional output of a plot, default set to FALSE.
#' @param nPCs is the number of PCs used to measure the distance, default is set to 3.
#'
#' @return list List containing the association plot and the computed ari
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
        assay_names=NULL,
        plot=FALSE,
        nPCs=3
){
    if (!is.null(assay_names)){
        normalization=assay_names
    }else{
        normalization=names(pca)
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
            datasets))

    ### Plot
    if (isTRUE(plot)){
    dataSets.colors <- wes_palette(
        n = length(normalization),
        name = "GrandBudapest1")[c(1,2,4,3)]
    p=ggplot(pcs.ari , aes(x = datasets, y = ari, fill = datasets)) +
        geom_col() +
        ylab("ARI") +
        xlab('') +
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', size = 1),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12))
    return(list(plot=p,ari=pcs.ari))
    }else{
        return(ari=pcs.ari)}
}

