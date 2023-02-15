
#' is used to compute the linear regression of a variable to the PCA of the data
#'
#'
#' @param pca PCs of the dataset that will be used in the plot
#' @param normalization All the available assays for the data (i.e. normalizations methods)
#' @param regression_var The regression variable that will be computed to the PCAof the data
#' @param variable.name The label of the variable that will be used on the PCA plot
#' @param nb_pca_comp The number of components of the PCA used to compute the regression
#'
#' @return list List containing the association plot and the computed regression
#' @importFrom stats lm
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr rename mutate
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export


regression_pc<-function(
    pca,
    normalization,
    regression_var,
    nb_pca_comp=10
){
    ### Compute the regression
    lreg.pcs<- lapply(
        normalization,
        function(x){
            pcs <- pca[[x]]$sing.val$u
            rSquared <- sapply(
                1:nb_pca_comp,
                function(y) {
                    lm.ls <- summary(lm(
                        regression_var ~ pcs[, 1:y])
                    )$r.squared
                })
        })
    names(lreg.pcs) <- normalization

    ### Plot the association between the variable and the PC using the computed regression
    pcs.lnreg <- as.data.frame(lreg.pcs) %>%
        rename(
            'Raw counts' =  normalization[1],
            FPKM =  normalization[2],
            FPKM.UQ =  normalization[3],
            'RUV-PRPS' =  normalization[4]) %>%
        mutate(pcs = c(1:10)) %>%
        pivot_longer(
            -pcs,
            names_to = 'datasets',
            values_to = 'r.sq') %>%
        mutate(
            datasets = factor(
                datasets,
                levels = c(
                    'Raw counts',
                    'FPKM',
                    'FPKM.UQ',
                    'RUV-PRPS')))
    # color
    dataSets.colors <- wes_palette(
        n = 4,
        name = "GrandBudapest1")[c(1,2,4,3)]
    names(dataSets.colors) <- c(
        'FPKM',
        'Quantile',
        'UQ',
        'RUV-PRPS')
    p=ggplot(pcs.lnreg, aes(x = pcs, y = r.sq, group = datasets)) +
        geom_line(aes(color = datasets), size = .5) +
        geom_point(aes(color = datasets), size = 2) +
        xlab('PCs') + ylab (expression("R"^"2")) +
        scale_color_manual(
            values = c(dataSets.colors3),
            name = 'Datasets',
            labels = c('Raw counts', 'FPKM','FPKM.UQ', 'RUV-PRPS')) +
        scale_x_continuous(
            breaks = (1:10),
            labels = c('PC1', paste0('PC1:', 2:10)) ) +
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
            legend.title = element_text(size = 14))

    return(list(plot=p,reg=pcs.lnreg))
}
