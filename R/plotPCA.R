#' is used to plot the pairwise plots of the first principal components of the assays of a SummarizedExperiment object.
#'
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to plot the PCA. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param fast.pca Logical. Indicates whether to calculate a specific number of PCs instead of the full range
#' to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to be calculated for the fast pca process, by default is set to 10.
#' @param ncol.plot is the argument of gtable for the layout specifying ncol, by default it is set to 4.
#' @param color The color of the variable that will be used on the PCA plot
#' @param strokeSize geom_point aesthetics
#' @param pointSize geom_point aesthetics
#' @param strokeColor geom_point aesthetics
#' @param alpha geom_point aesthetics
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return plot PCA plot of the data colored by one variable
#' @importFrom ggpubr ggarrange as_ggplot
#' @importFrom patchwork plot_spacer plot_layout
#' @import ggplot2
#' @export


plotPCA <- function(
        se.obj,
        assay.names = 'All',
        variable,
        color = NULL,
        fast.pca = TRUE,
        nb.pcs = 3,
        ncol.plot = 4,
        strokeSize = .2,
        pointSize = 1.5,
        strokeColor = 'gray30',
        alpha = .5,
        assess.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The plotPCA function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (is.null(nb.pcs)) {
        stop('To plot the PCA, the number of PCs (nb.pcs) must be specified.')
    } else if (is.null(assay.names)) {
        stop('Please provide at least an assay name.')
    } else if (is.null(variable)) {
        stop('Please provide a categorical variable.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0(
            'The ',
            variable,
            ', is a numeric, but this should a categorical variable'
        ))
    }

    # find assays ####
    if (length(assay.names) == 1 && assay.names == 'All') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else {
        assay.names <- as.factor(unlist(assay.names))
    }
    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = 'none',
            verbose = verbose
        )
    }
    # pca analysis ####
    p.pca <- lapply(
        levels(assay.names),
        function(x) {
            if (fast.pca) {
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]]) )
                    stop('To plot the PCA, the fast PCA must be computed first on the assay ', x, ' .')
                pca.data <- se.obj@metadata[['metric']][[x]][['fastPCA']]$svd$u[, 1:nb.pcs]
                pc.var <- se.obj@metadata[['metric']][[x]][['fastPCA']]$percentage.variation
            } else {
                if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To plot the PCA, the PCA must be computed first on the assay ', x, ' .')
                pca.data <- se.obj@metadata[['metric']][[x]][['PCA']]$svd$u[, 1:nb.pcs]
                pc.var <- se.obj@metadata[['metric']][[x]][['PCA']]$percentage.variation
            }
            pair.pcs <- combn(ncol(pca.data), 2)
            var <- colData(se.obj)[ , variable]
            p.per.data <- lapply(
                1:ncol(pair.pcs),
                function(i){
                    plot1 <- ggplot(
                        mapping = aes(
                            x = pca.data[, pair.pcs[1, i]],
                            y = pca.data[, pair.pcs[2, i]] ) ) +
                        geom_point(aes(fill = se.obj@colData[, variable]),
                                   color = strokeColor,
                                   pch = 21,
                                   stroke = strokeSize,
                                   size = pointSize,
                                   alpha = alpha) +
                        scale_x_continuous(
                            name = paste0('PC', pair.pcs[1, i], ' (', pc.var[pair.pcs[2, i]], '%)'),
                            breaks = scales::pretty_breaks(n = 5)) +
                        scale_y_continuous(
                            name = paste0('PC', pair.pcs[2, i], ' (', pc.var[pair.pcs[1, i]], '%)'),
                            breaks = scales::pretty_breaks(n = 5)) +
                        ggtitle(x) +
                        theme_pubr() +
                        theme(legend.position = "none",
                              legend.background = element_blank(),
                              legend.text = element_text(size = 12),
                              legend.title = element_text(size = 14),
                              legend.key = element_blank(),
                              axis.text.x = element_text(size = 10),
                              axis.text.y = element_text(size = 10),
                              axis.title.x = element_text(size = 12),
                              axis.title.y = element_text(size = 12)) +
                        scale_fill_discrete(name = variable) +
                        guides(fill = guide_legend(override.aes = list(size = 3, shape = 21) ) )

                    dens1 <- ggplot(mapping = aes(x = pca.data[, pair.pcs[1, i]], fill = se.obj@colData[, variable])) +
                        geom_density(alpha = 0.4) +
                        theme_void() +
                        theme(legend.position = "none",
                              legend.text = element_text(size = 12),
                              legend.title = element_text(size = 14)) +
                        scale_fill_discrete(name = variable) +
                        guides(fill = guide_legend(override.aes = list(size = 3)))

                    dens2 <- ggplot(mapping = aes(x = pca.data[, pair.pcs[2, i]], fill = se.obj@colData[, variable])) +
                        geom_density(alpha = 0.4) +
                        theme_void() +
                        theme(legend.position = "none",
                              legend.text = element_text(size = 12),
                              legend.title = element_text(size = 14)) +
                        coord_flip() +
                        scale_fill_discrete(name = variable) +
                        guides(fill = guide_legend(override.aes = list(size = 3)))

                    dens1 + plot_spacer() + plot1 + dens2 +
                        plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

                })
            p.per.data
        })
    names(p.pca) <- levels(assay.names)
    ## Prepare plot
    p <- p.pca[[levels(assay.names)[1]]]
    if (length(assay.names) > 1) {
        for (n in 2:length(assay.names)) {
            p = c(p, p.pca[[levels(assay.names)[n]]])
        }
    }
    plot <- ggarrange(
        plotlist = p,
        common.legend = TRUE,
        legend = "bottom",
        nrow = length(levels(assay.names)),
        ncol = nb.pcs)
    printColoredMessage(message = '------------The plotPCA function starts:',
                        color = 'white',
                        verbose = verbose)
    return(plot = as_ggplot(plot))
}
