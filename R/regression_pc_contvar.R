
#' is used to compute the linear regression of a continuous variable and the first cumulative PCs
#' of the gene expression (assay) of a SummarizedExperiment class object.
#'
#'
#' @param pca PCA components of a SummarizedExperiment variable that will be used in the plot.
#' @param cont_var Vector of a categorical variable such as sample types
#' (i.e. biological subtypes) or batches.
#' @param cont_var_label String or vector of strings of the label of categorical variable(s) such as
#' sample types or batches from colData(se).
#' @param assay_names Optional selection of names of the assays to compute the PCA.
#' @param nb_pca_comp The number of components of the PCA used to compute the regression.
#'
#' @return list List containing the associated plot and the computed regression for the continuous variable.
#' @importFrom stats lm
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr rename mutate
#' @importFrom tidyr pivot_longer %>%
#' @import ggplot2
#' @export


regression_pc_contvar<-function(
    pca,
    cont_var,
    cont_var_label,
    assay_names=NULL,
    nb_pca_comp=10
){

    ## Assays
    if (!is.null(assay_names)){
        normalization=as.factor(unlist(assay_names))
    }else{
        normalization=as.factor(names(pca))
    }
    ### Compute the regression
    lreg.pcs<- lapply(
        normalization,
        function(x){
            pcs <- pca[[x]]$sing.val$u
            rSquared <- sapply(
                1:nb_pca_comp,
                function(y) {
                    lm.ls <- summary(lm(
                        cont_var ~ pcs[, 1:y])
                    )$r.squared
                })
        })
    names(lreg.pcs) <- normalization

    ### Plot the association between the variable and the PC using the computed regression
    pcs<-datasets<-r.sq<-NULL
    pcs.lnreg=cbind(as.data.frame(lreg.pcs),pcs=(1:10))
    ## length of assays
    assays_nb=length(normalization)
    pcs.lnreg = pcs.lnreg %>% pivot_longer( -(assays_nb+1),
    names_to = 'datasets',values_to = 'r.sq') %>%
        mutate(datasets = factor(
            datasets,levels=normalization))
    # color
    dataSets.colors <- wes_palette(
        n = assays_nb,
        name = "GrandBudapest1")[seq(1:assays_nb)]
    names(dataSets.colors) <- normalization
    p=ggplot(pcs.lnreg, aes(x = pcs, y = r.sq, group = datasets)) +
        geom_line(aes(color = datasets), size = 1) +
        geom_point(aes(color = datasets), size = 3) +
        xlab('PCs') + ylab (expression("R"^"2")) +
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
    ggtitle(paste("Regression of ",cont_var_label," and the first cumulative PCs",sep=""))

    return(list(plot=p,reg=pcs.lnreg,cont_var_label=cont_var_label))
}
