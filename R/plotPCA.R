#' is used to plot principal components.

#' @description
#' is used to plot the pairwise plots of the first principal components of the assays of a SummarizedExperiment object.

#' @param se.obj A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s) of the assay(s) of the
#' SummarizedExperiment class object to plot the PCA. By default all the assays of the SummarizedExperiment class object
#' will be selected.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param fast.pca Logical. Indicates whether to calculate a specific number of PCs instead of the full range to speed up
#' the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to be calculated for the fast pca process, by default is set to 10.
#' @param plot.type Symbol.
#' @param point.color The color of the variable that will be used on the PCA plot
#' @param point.size geom_point aesthetics
#' @param stroke.color geom_point aesthetics
#' @param stroke.size geom_point aesthetics
#' @param alpha.point geom_point aesthetics
#' @param alpha.density geom_point aesthetics
#' @param ncol.plot is the argument of gtable for the layout specifying ncol, by default it is set to 4.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return plot PCA plot of the data colored by one variable

#' @importFrom ggpubr ggarrange theme_pubr
#' @importFrom patchwork plot_spacer plot_layout
#' @importFrom tidyr pivot_longer
#' @importFrom utils combn
#' @import ggplot2
#' @export

plotPCA <- function(
        se.obj,
        assay.names = 'All',
        variable,
        fast.pca = TRUE,
        nb.pcs = 3,
        plot.type = 'scatter',
        point.color = NULL,
        point.size = 1.5,
        stroke.color = 'gray30',
        stroke.size = .2,
        alpha.point = .5,
        alpha.density = .5,
        ncol.plot = 4,
        assess.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The plotPCA function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (is.null(assay.names)){
        stop('The "assay.names" cannot be empty.')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0('The ', variable, 'must be a categorical variable'))
    } else if (is.null(nb.pcs)){
        stop('The "nb.pcs" cannot be empty.')
    } else if (!plot.type %in% c('scatter', 'boxplot')){
        stop('The "plot.type" must be one of "scatter" or "boxplot"')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'All') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- as.factor(unlist(assay.names))

    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = 'none',
            verbose = verbose)
    }
    # pca plots ####
    if(plot.type == 'scatter'){
        printColoredMessage(message = '-- Plot PCA for individual assays:',
                            color = 'magenta',
                            verbose = verbose)
        p.pca <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('Obtain the first ', nb.pcs, ' PCs of ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose)
                if (fast.pca) {
                    if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]]))
                        stop('To plot the PCA, the fast PCA must be computed first on the assay ', x, ' .')
                    pca.data <- se.obj@metadata[['metric']][[x]][['fastPCA']]$svd$u
                    pc.var <- se.obj@metadata[['metric']][[x]][['fastPCA']]$percentage.variation
                } else {
                    if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]]))
                        stop('To plot the PCA, the PCA must be computed first on the assay ', x, ' .')
                    pca.data <- se.obj@metadata[['metric']][[x]][['PCA']]$svd$u[, 1:nb.pcs]
                    pc.var <- se.obj@metadata[['metric']][[x]][['PCA']]$percentage.variation
                }
                if(ncol(pca.data) < nb.pcs){
                    printColoredMessage(
                        message = paste0('The number of PCs of the assay', x, 'are ', ncol(pca.data), '.'),
                        color = 'blue',
                        verbose = verbose)
                    stop(paste0('The number of PCs of the assay ', x, ' are less than', nb.pcs, '.',
                                'Re-run the computePCA function with nb.pcs = ', nb.pcs, '.'))
                }
                pca.data <- pca.data[ , seq_len(nb.pcs)]
                pair.pcs <- combn(ncol(pca.data), 2)
                var <- colData(se.obj)[, variable]
                printColoredMessage(
                    message = paste0('Create all possible pairwise scatter plots of the first ', nb.pcs, ' PCs.'),
                    color = 'blue',
                    verbose = verbose)
                p.per.data <- lapply(
                    1:ncol(pair.pcs),
                    function(i) {
                        plot1 <- ggplot(mapping = aes(x = pca.data[, pair.pcs[1, i]], y = pca.data[, pair.pcs[2, i]])) +
                            geom_point(
                                aes(fill = se.obj@colData[, variable]),
                                color = stroke.color,
                                pch = 21,
                                stroke = stroke.size,
                                size = point.size,
                                alpha = alpha.point) +
                            scale_x_continuous(
                                name = paste0('PC', pair.pcs[1, i], ' (', pc.var[pair.pcs[2, i]], '%)'),
                                breaks = scales::pretty_breaks(n = 5)) +
                            scale_y_continuous(
                                name = paste0('PC', pair.pcs[2, i], ' (', pc.var[pair.pcs[1, i]], '%)'),
                                breaks = scales::pretty_breaks(n = 5)) +
                            ggtitle(x) +
                            theme_pubr() +
                            theme(
                                legend.position = "none",
                                legend.background = element_blank(),
                                legend.text = element_text(size = 12),
                                legend.title = element_text(size = 14),
                                legend.key = element_blank(),
                                axis.text.x = element_text(size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.x = element_text(size = 12),
                                axis.title.y = element_text(size = 12)) +
                            scale_fill_discrete(name = variable) +
                            guides(fill = guide_legend(override.aes = list(size = 3, shape = 21)))
                        dense.x <-
                            ggplot(mapping = aes(x = pca.data[, pair.pcs[1, i]], fill = se.obj@colData[, variable])) +
                            geom_density(alpha = 0.4) +
                            theme_void() +
                            theme(
                                legend.position = "none",
                                legend.text = element_text(size = 12),
                                legend.title = element_text(size = 14)) +
                            scale_fill_discrete(name = variable) +
                            guides(fill = guide_legend(override.aes = list(size = 3)))

                        dense.y <- ggplot(mapping = aes(x = pca.data[, pair.pcs[2, i]], fill = se.obj@colData[, variable])) +
                            geom_density(alpha = 0.4) +
                            theme_void() +
                            theme(
                                legend.position = "none",
                                legend.text = element_text(size = 12),
                                legend.title = element_text(size = 14)) +
                            coord_flip() +
                            scale_fill_discrete(name = variable) +
                            guides(fill = guide_legend(override.aes = list(size = 3)))

                        dense.x +
                            plot_spacer() +
                            plot1 + dense.y +
                            plot_layout(
                                ncol = 2,
                                nrow = 2,
                                widths = c(4, 1),
                                heights = c(1, 4))
                    })
                p.per.data
            })
        names(p.pca) <- levels(assay.names)
        # prepare plot ####
        printColoredMessage(message = '-- Put all the plots togather:',
                            color = 'magenta',
                            verbose = verbose)
        p <- p.pca[[levels(assay.names)[1]]]
        if (length(assay.names) > 1) {
            for (n in 2:length(assay.names))
                p = c(p, p.pca[[levels(assay.names)[n]]])
        }
        plot <- ggarrange(
            plotlist = p,
            common.legend = TRUE,
            legend = "bottom",
            nrow = length(levels(assay.names)),
            ncol = ncol(combn(nb.pcs, 2))
        )
    } else{
        printColoredMessage(message = '-- Plot PCA for individual assays:',
                            color = 'magenta',
                            verbose = verbose)
        p.pca <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('Obtain the first ', nb.pcs, ' PCs of ', x, ' data.'),
                    color = 'blue',
                    verbose = verbose)
                if (fast.pca) {
                    if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]]))
                        stop('To plot the PCA, the fast PCA must be computed first on the assay ', x, ' .')
                    pca.data <- se.obj@metadata[['metric']][[x]][['fastPCA']]$svd$u
                    pc.var <- se.obj@metadata[['metric']][[x]][['fastPCA']]$percentage.variation
                } else {
                    if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]]))
                        stop('To plot the PCA, the PCA must be computed first on the assay ', x, ' .')
                    pca.data <- se.obj@metadata[['metric']][[x]][['PCA']]$svd$u[, 1:nb.pcs]
                    pc.var <- se.obj@metadata[['metric']][[x]][['PCA']]$percentage.variation
                }
                if(ncol(pca.data) < nb.pcs){
                    printColoredMessage(
                        message = paste0('The number of PCs of the assay', x, 'are ', ncol(pca.data), '.'),
                        color = 'blue',
                        verbose = verbose)
                    stop(paste0('The number of PCs of the assay ', x, ' are less than', nb.pcs, '.',
                                'Re-run the computePCA function with nb.pcs = ', nb.pcs, '.'))
                }
                var <- NULL
                pca.data <- as.data.frame(pca.data[ , seq_len(nb.pcs)])
                colnames(pca.data) <- paste0('PC', seq_len(nb.pcs), ' (',pc.var[seq_len(nb.pcs)], '%)')
                pca.data$var <- colData(se.obj)[, variable]
                printColoredMessage(
                    message = paste0('Create boxplots of the first ', nb.pcs, ' PCs.'),
                    color = 'blue',
                    verbose = verbose)
                pca.data <- tidyr::pivot_longer(
                    data = pca.data,
                    -var,
                    names_to = 'PCs',
                    values_to = 'PC')
                pca.data <- as.data.frame(pca.data)
                plot.p <- ggplot(pca.data, aes(x = var, y =  PC, fill = var)) +
                    geom_boxplot() +
                    facet_wrap(~PCs, scales = 'free') +
                    xlab(variable) +
                    ylab('PC') +
                    ggtitle(x) +
                    # scale_fill_manual(values = studies.color, name = 'Studies') +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', size = 1),
                        axis.title.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        plot.title = element_text(size = 18),
                        axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1),
                        axis.text.y = element_text(size = 10),
                        strip.text.x = element_text(size = 12),
                        legend.position = 'none')
                plot.p
            })
        names(p.pca) <- levels(assay.names)
        # prepare plot ####
        printColoredMessage(message = '-- Put all the plots togather:',
                            color = 'magenta',
                            verbose = verbose)
        plot <- ggarrange(
            plotlist = p.pca,
            nrow = length(levels(assay.names)),
            ncol = 1 )
    }
    printColoredMessage(message = '------------The plotPCA function finished.',
                        color = 'white',
                        verbose = verbose)
    return(plot)
}
