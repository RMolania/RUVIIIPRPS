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
#' @importFrom ggpubr get_legend as_ggplot
#' @importFrom grid textGrob
#' @importFrom cowplot axis_canvas ggdraw insert_xaxis_grob insert_yaxis_grob
#' @import ggplot2 scales
#' @importFrom gridExtra grid.arrange
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
        stop(paste0('The ', variable,', is a numeric, but this should a categorical variable'))
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
            verbose = verbose)
    }

    if (fast.pca) {
        ppca <- lapply(
            levels(assay.names),
            function(x) {
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]])) {
                    stop('To plot the PCA, the fast PCA must be computed first on the assay ', x, ' .')
                }
                pca.data  <- se.obj@metadata[['metric']][[x]][['fastPCA']]$sing.val$u[,1:nb.pcs]
                pc.var <- se.obj@metadata[['metric']][[x]][['fastPCA']]$sing.val$
                pair.pcs <- combn(ncol(pca.data), 2)
                pList <- list()
                for(i in 1:ncol(pair.pcs)){
                    if(i == 1){
                        p <- ggplot(mapping = aes(x = pca.data[,pair.pcs[1,i]], y = pca.data[,pair.pcs[2,i]], fill = var)) +
                            geom_point(aes(fill = se.obj@colData[, variable]), pch = 21, color = strokeColor, stroke = strokeSize, size = pointSize, alpha = alpha) +
                            scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            xlab(paste0('PC', x, ' (', pc.var[x], '%)')) +
                            ylab(paste0('PC', y, ' (', pc.var[y], '%)')) +
                            theme(legend.position = "right",
                                  panel.background = element_blank(),
                                  axis.line = element_line(colour = "black", size = 1.1),
                                  legend.background = element_blank(),
                                  legend.text = element_text(size = 12),
                                  legend.title = element_text(size = 14),
                                  legend.key = element_blank(),
                                  axis.text.x = element_text(size = 10),
                                  axis.text.y = element_text(size = 10),
                                  axis.title.x = element_text(size = 14),
                                  axis.title.y = element_text(size = 14),
                                  aspect.ratio = 1) +
                            guides(fill = guide_legend(override.aes = list(size = 4)))
                        if (!is.null(color))
                            p <- p + scale_fill_manual(name = variable, values = color)
                        le <- get_legend(p)
                    } else{
                        p <- ggplot(mapping = aes(x = pcs[,pair.pcs[1,i]], y = pcs[,pair.pcs[2,i]], fill = var)) +
                            geom_point(aes(fill = se.obj@colData[, variable]), pch = 21, color = strokeColor, stroke = strokeSize, size = pointSize,alpha = alpha) +
                            scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            xlab(paste0('PC', x, ' (', pc.var[x], '%)')) +
                            ylab(paste0('PC', y, ' (', pc.var[y], '%)')) +
                            theme(panel.background = element_blank(),
                                  axis.line = element_line(colour = "black", size = 1.1),
                                  legend.position = "none",
                                  axis.text.x = element_text(size = 10),
                                  axis.text.y = element_text(size = 10),
                                  axis.title.x = element_text(size = 14),
                                  axis.title.y = element_text(size = 14),
                                  aspect.ratio = 1)
                        if (!is.null(color))
                            p <- p + scale_fill_manual(name = variable, values = color)
                    }
                    p <- p + theme(legend.position = "none")
                    xdens <- axis_canvas(p, axis = "x")+
                        geom_density( mapping = aes(x = pcs[,x], fill = se.obj@colData[, variable]), alpha = 0.7, size = 0.2) +
                        theme(legend.position = "none")
                    if (!is.null(color))
                        xdens <- xdens + scale_fill_manual(name = variable, values = color)
                    ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
                        geom_density(mapping = aes(x = pcs[,y], fill = se.obj@colData[, variable]), alpha = 0.7, size = 0.2) +
                        theme(legend.position = "none")+
                        coord_flip()
                    if (!is.null(color))
                        ydens <- ydens + scale_fill_manual(name = variable, values = color)
                    p1 <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
                    p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
                    pList[[i]] <- ggdraw(p2)
                }
                pList[[i+1]] <- le
                return(pList)

            })
        names(ppca) <- levels(assay.names)

    } else {
        ppca <- lapply(
            levels(assay.names),
            function(x) {
                if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To plot the PCA, the PCA must be computed first on the assay ', x, ' .')
                pca.data  <- se.obj@metadata[['metric']][[x]][['PCA']]$sing.val$u[,1:nb.pcs]
                pc.var <- se.obj@metadata[['metric']][[x]][['PCA']]$sing.val$
                    pair.pcs <- combn(ncol(pca.data), 2)
                pList <- list()
                for(i in 1:ncol(pair.pcs)){
                    if(i == 1){
                        p <- ggplot(mapping = aes(x = pca.data[,pair.pcs[1,i]], y = pca.data[,pair.pcs[2,i]], fill = var)) +
                            geom_point(aes(fill = se.obj@colData[, variable]), pch = 21, color = strokeColor, stroke = strokeSize, size = pointSize, alpha = alpha) +
                            scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            xlab(paste0('PC', x, ' (', pc.var[x], '%)')) +
                            ylab(paste0('PC', y, ' (', pc.var[y], '%)')) +
                            theme(legend.position = "right",
                                  panel.background = element_blank(),
                                  axis.line = element_line(colour = "black", size = 1.1),
                                  legend.background = element_blank(),
                                  legend.text = element_text(size = 12),
                                  legend.title = element_text(size = 14),
                                  legend.key = element_blank(),
                                  axis.text.x = element_text(size = 10),
                                  axis.text.y = element_text(size = 10),
                                  axis.title.x = element_text(size = 14),
                                  axis.title.y = element_text(size = 14),
                                  aspect.ratio = 1) +
                            guides(fill = guide_legend(override.aes = list(size = 4)))
                        if (!is.null(color))
                            p <- p + scale_fill_manual(name = variable, values = color)
                        le <- get_legend(p)
                    } else{
                        p <- ggplot(mapping = aes(x = pcs[,pair.pcs[1,i]], y = pcs[,pair.pcs[2,i]], fill = var)) +
                            geom_point(aes(fill = se.obj@colData[, variable]), pch = 21, color = strokeColor, stroke = strokeSize, size = pointSize,alpha = alpha) +
                            scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            xlab(paste0('PC', x, ' (', pc.var[x], '%)')) +
                            ylab(paste0('PC', y, ' (', pc.var[y], '%)')) +
                            theme(panel.background = element_blank(),
                                  axis.line = element_line(colour = "black", size = 1.1),
                                  legend.position = "none",
                                  axis.text.x = element_text(size = 10),
                                  axis.text.y = element_text(size = 10),
                                  axis.title.x = element_text(size = 14),
                                  axis.title.y = element_text(size = 14),
                                  aspect.ratio = 1)
                        if (!is.null(color))
                            p <- p + scale_fill_manual(name = variable, values = color)
                    }
                    p <- p + theme(legend.position = "none")
                    xdens <- axis_canvas(p, axis = "x")+
                        geom_density( mapping = aes(x = pcs[,x], fill = se.obj@colData[, variable]), alpha = 0.7, size = 0.2) +
                        theme(legend.position = "none")
                    if (!is.null(color))
                        xdens <- xdens + scale_fill_manual(name = variable, values = color)
                    ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
                        geom_density(mapping = aes(x = pcs[,y], fill = se.obj@colData[, variable]), alpha = 0.7, size = 0.2) +
                        theme(legend.position = "none")+
                        coord_flip()
                    if (!is.null(color))
                        ydens <- ydens + scale_fill_manual(name = variable, values = color)
                    p1 <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
                    p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
                    pList[[i]] <- ggdraw(p2)
                }
                pList[[i+1]] <- le
                return(pList)
            })
        names(ppca) <- levels(assay.names)
    }
    ## Prepare plot
    p = ppca[[levels(assay.names)[1]]]
    if (length(assay.names) > 1) {
        for (n in 2:length(assay.names)) {
            p = c(p, ppca[[levels(assay.names)[n]]])
        }
    }
    if (fast.pca) {
        title = paste0("FastPCA ordered as ", paste(levels(assay.names), collapse = ","))
    } else{
        title = paste0("PCA ordered as ", paste(levels(assay.names), collapse = ","))
    }
    plot = do.call(grid.arrange,
                   c(p,
                     ncol = ncol.plot,
                     top = title))

    return(plot = as_ggplot(plot))
}
