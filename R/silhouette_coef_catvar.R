#' is used to compute the Silhouette coefficient from the PC of all assays
#' given a categorical variable
#'
#'
#' @param pca PCs of the dataset that will be used
#' @param cat_var is a categorical variable such as sample types or batches
#' @param assay_names Optional selection of names of the assays to compute the PCA
#' @param nPCs is the number of PCs used to measure the distance
#'
#' @return list List containing the association plot and the computed silhouette
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr rename mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom stats dist
#' @importFrom cluster silhouette
#' @import ggplot2
#' @export

silhouette_coef_catvar<-function(
        pca,
        assay_names=NULL,
        cat_var,
        nPCs=3
){
    if (!is.null(assay_names)){
        normalization=assay_names
    }else{
        normalization=names(pca)
    }
    # Silhouette coefficients on all assays
    silCoef <- lapply(
        normalization,
        function(x){
            silhouette_coef_catvar_single_assay <- function(pca,
                                                            cat_var,
                                                            nPCs){
                d.matrix <- as.matrix(dist(pca[, seq_len(nPCs)]))
                avg=summary(silhouette(
                    as.numeric(as.factor(cat_var)),
                    d.matrix))$avg.width
                return(avg)
            }
            sil=silhouette_coef_catvar_single_assay(
                pca[[x]]$sing.val$u,
                cat_var,
                nPCs)
            sil
        })
    names(silCoef) <- normalization
    everything<-datasets<-silh.coeff<-NULL
    pcs.silCoef <- as.data.frame(silCoef)
    pcs.silCoef =  pcs.silCoef %>% pivot_longer(
        everything(),
        names_to = 'datasets',
        values_to = 'silh.coeff') %>% mutate(datasets = factor(
            datasets))

    ### Plot
    # color
    dataSets.colors <- wes_palette(
        n = length(normalization),
        name = "GrandBudapest1")[c(1,2,4,3)]
    p=ggplot(pcs.silCoef , aes(x = datasets, y = silh.coeff, fill = datasets)) +
        geom_col() +
        ylab("Silhouette coefficient") +
        xlab('') +
        scale_fill_manual(values = dataSets.colors, guide = 'none')+
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', size = 1),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12))

    return(list(plot=p,silh.coeff=pcs.silCoef))
}

