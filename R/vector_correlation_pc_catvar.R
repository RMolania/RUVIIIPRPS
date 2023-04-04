#' is used to compute the vector correlation between the first cumulative PCs of the gene expression (assay)
#' of a SummarizedExperiment class object and a categorical variable (i.e. batch).
#'
#' @param pca PCA components of a SummarizedExperiment variable that will be used in the plot.
#' @param cat_var Vector of a categorical variable such as sample types
#' (i.e. biological subtypes) or batches.
#' @param cat_var_label String or vector of strings of the label of categorical variable(s) such as
#' sample types or batches from colData(se).
#' @param assay_names Optional string or list of strings for selection of the names
#' of the assays of the SummarizedExperiment class object to compute the correlation.
#' @param plot_output Indicates whether to plot the results
#' @param nb_pca_comp The number of components of the PCA used to compute the regression.
#'
#' @return list List containing the associated plot and correlation for the categorical variable.
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom SummarizedExperiment assays assay
#' @importFrom matrixTests row_oneway_equalvar
#' @importFrom fastDummies dummy_cols
#' @importFrom wesanderson wes_palette
#' @importFrom stats cancor
#' @import ggplot2
#' @export

vector_correlation_pc_catvar<-function(
        pca,
        cat_var,
        cat_var_label,
        assay_names=NULL,
        plot_output = TRUE,
        nb_pca_comp=10
){

    ## Assays
    if (!is.null(assay_names)){
        normalization=assay_names
    }else{
        normalization=names(pca)
    }

    ## Correlation
    message(paste0(
        'Performing vector correlation between PCs and the ',
        cat_var_label, ' variable \n')
    )

    ### across all samples
    catvar.dummies <- dummy_cols(cat_var)
    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]

    cca.all<- lapply(
        normalization,
        function(x){
            pcs <- pca[[x]]$sing.val$u
            sapply(
                1:nb_pca_comp,
                function(y){
                    ## Canonical correlations
                    cca.fcch <- cancor(
                        x = pcs[, 1:y, drop = FALSE],
                        y = catvar.dummies)
                    1 - prod(1 - cca.fcch$cor^2)
                })
        })
    names(cca.all)=normalization


    ### Plot the association between the variable and the PC using the computed regression
    pcs<-datasets<-pcs<-cca.coef<-NULL
    pcs.cca=as.data.frame(cca.all)
    pcs.cca = pcs.cca %>% mutate(pcs=c(1:nb_pca_comp)) %>% pivot_longer( -pcs,
                                            names_to = 'datasets',values_to = 'cca.coef') %>%
        mutate(datasets = factor(
            datasets))
    ## length of assays
    assays_nb=length(normalization)
    # color
    dataSets.colors <- wes_palette(
        n = assays_nb,
        name = "GrandBudapest1")[c(1,2,4,3)]
    names(dataSets.colors) <- normalization
    p=ggplot(pcs.cca, aes(x = pcs, y = cca.coef, group = datasets)) +
        geom_line(aes(color = datasets), size = 1) +
        geom_point(aes(color = datasets), size = 3) +
        xlab('PCs') + ylab ("Vector correlation") +
        scale_color_manual(
            values = c(dataSets.colors),
            name = 'Datasets',
            labels = normalization) +
        scale_x_continuous(
            breaks = (1:nb_pca_comp),
            labels = c('PC1', paste0('PC1:', 2:nb_pca_comp)) ) +
        scale_y_continuous(
            breaks = scales::pretty_breaks(n = 5),
            limits = c(0,1)) +
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', size = 1),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 12, angle = 35, hjust = 1),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 14))+
        ggtitle(paste("Vector correlation between ",cat_var_label," and the first cumulative PCs",sep=""))

    return(list(plot=p,cca=pcs.cca,cat_var_label=cat_var_label))
}
