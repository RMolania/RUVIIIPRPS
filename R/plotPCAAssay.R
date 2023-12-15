#' is used to plot the pairwise plots of the first PCA of the gene expression
#' of a single assay of the SummarizedExperiment class object.
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the differential expression ANOVA analysis.
#' @param pca_x The PCA of a SummarizedExperiment object that will be used to plot.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param fast.pca Logical. Indicates whether to calculate a specific number of PCs instead of the full range
#' to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to be calculated for the fast pca process, by default is set to 10.
#' @param color The color of the variable that will be used on the PCA plot
#' @param strokeSize geom_point aesthetics
#' @param pointSize geom_point aesthetics
#' @param strokeColor geom_point aesthetics
#' @param alpha geom_point aesthetics
#'
#' @return plot PCA plot of the data colored by one variable
#' @importFrom ggpubr get_legend as_ggplot
#' @importFrom grid textGrob
#' @importFrom cowplot axis_canvas ggdraw insert_xaxis_grob insert_yaxis_grob
#' @import ggplot2 scales
#' @importFrom utils combn
#' @importFrom gridExtra grid.arrange


plotPCAassay <- function(
        se.obj,
        pca_x,
        variable,
        nb.pcs = 3,
        fast.pca = TRUE,
        color = NULL,
        strokeSize = .2,
        pointSize = 1.5,
        strokeColor = 'gray30',
        alpha = .5
) {
    pcs <- pca_x$sing.val$u[,1:nb.pcs]
    pc.var <- pca_x$variation
    pc.var.nb.pcs <- length(pca_x$variation)
    pair.pcs <- combn(ncol(pcs), 2)
    pList <- list()
    var <- se.obj@colData[, variable]
    var.label <- paste0(variable)

    for(i in 1:ncol(pair.pcs)){
        if(i == 1){
            p <- ggplot(mapping = aes(x = pcs[,pair.pcs[1,i]], y = pcs[,pair.pcs[2,i]], fill = var)) +
                geom_point(aes(fill = var), pch = 21, color = strokeColor, stroke = strokeSize, size = pointSize, alpha = alpha) +
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
                p <- p + scale_fill_manual(name = var.label, values = color)
            le <- get_legend(p)
        } else{
            p <- ggplot(mapping = aes(x = pcs[,pair.pcs[1,i]], y = pcs[,pair.pcs[2,i]], fill = var)) +
                geom_point(aes(fill = var), pch = 21, color = strokeColor, stroke = strokeSize, size = pointSize,alpha = alpha) +
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
                    p <- p + scale_fill_manual(name = var.label, values = color)
        }
        p <- p + theme(legend.position = "none")
        xdens <- axis_canvas(p, axis = "x")+
            geom_density( mapping = aes(x = pcs[,x],fill = var), alpha = 0.7, size = 0.2) +
            theme(legend.position = "none")
            if (!is.null(color))
                xdens <- xdens + scale_fill_manual(name = var.label,values = color)
        ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
            geom_density(mapping = aes(x = pcs[,y], fill = var), alpha = 0.7, size = 0.2) +
            theme(legend.position = "none")+
            coord_flip()
            if (!is.null(color))
                ydens <- ydens + scale_fill_manual(name = var.label,values = color)
        p1 <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
        p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
        pList[[i]] <- ggdraw(p2)
    }
    pList[[i+1]] <- le
    return(pList)
}
