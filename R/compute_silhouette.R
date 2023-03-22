#' is used to compute the Silhouette coefficient from the first PC of a SummarizedExperiment class
#' object given a categorical variable.
#'
#'
#' @param pca PCA components of a SummarizedExperiment variable.
#' @param cat_var Vector of a categorical variable such as sample types
#' (i.e. biological subtypes) or batches.
#' @param cat_var_label String or vector of strings of the label of categorical variable(s) such as
#' sample types or batches from colData(se).
#' @param assay_names Optional string or list of strings for selection of the names
#' of the assays of the SummarizedExperiment class object to compute the PCA.
#' @param plot Optional output of a plot, default set to FALSE.
#' @param nPCs is the number of PCs used to measure the distance, by default it is set to 3.
#'
#' @return list List containing the association plot and the computed silhouette
#' and the cat_var_label
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr rename mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom stats dist
#' @importFrom cluster silhouette
#' @import ggplot2
#' @export

compute_silhouette<-function(
        pca,
        cat_var,
        cat_var_label,
        assay_names=NULL,
        plot=FALSE,
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
    if (isTRUE(plot)){
        # color
        dataSets.colors <- wes_palette(
            n = length(normalization),
            name = "GrandBudapest1")[c(1,2,4,3)]
        names(dataSets.colors)=normalization
        p=ggplot(pcs.silCoef , aes(x = datasets, y = silh.coeff, fill = datasets)) +
            geom_point(aes(colour=datasets)) +
            ylab("Silhouette coefficient") +
            xlab('') +
            scale_color_manual(values = dataSets.colors) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', size = 1),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12))+
            ggtitle(paste("Silhouette coefficient computed on ",cat_var_label,sep=""))
        return(list(plot=p,silh.coeff=pcs.silCoef,cat_var_label=cat_var_label))
    }else{
        return(silh.coeff=pcs.silCoef)}
}

