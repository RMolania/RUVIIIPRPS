#' is used to plot the pairwise plots of the first PCA of the gene expression (assay)
#' of a SummarizedExperiment class object.
#'
#'
#' @param pca PCA components of a SummarizedExperiment variable that will be used in the plot.
#' @param assay_names Optional selection of names of the assays to compute the PCA.
#' @param cat_var Vector of a categorical variable such as sample types
#' (i.e. biological subtypes) or batches.
#' @param cat_var_label String or vector of strings of the label of categorical variable(s) such as
#' sample types or batches from colData(se).
#' @param nPCs is the number of PCs used to measure the distance, by default it is set to 3.
#' @param ncol_plot is the argument of gtable for the layout specifying ncol, by default it is set to 4.
#' @param color The color of the variable that will be used on the PCA plot
#' @param strokeSize geom_point aesthetics
#' @param pointSize geom_point aesthetics
#' @param strokeColor geom_point aesthetics
#' @param alpha geom_point aesthetics
#'
#' @return p PCA plot of the data colored by one variable
#' @importFrom ggpubr get_legend as_ggplot
#' @importFrom grid textGrob
#' @importFrom cowplot axis_canvas ggdraw insert_xaxis_grob insert_yaxis_grob
#' @import ggplot2 scales
#' @importFrom gridExtra grid.arrange
#' @export

plot_pca=function(
        pca,
        assay_names=NULL,
        cat_var,
        cat_var_label,
        nPCs=3,
        ncol_plot=4,
        color,
        strokeSize = .2,
        pointSize = 1.5,
        strokeColor = 'gray30',
        alpha = .5
){
    if (!is.null(assay_names)){
        normalizations=assay_names
    }else{
         normalizations=names(pca)
    }
    ppca <- lapply(
        normalizations,
        function(x){
            plot_pca_single_assay<-function(pca_x,
                                            cat_var,
                                            cat_var_label,
                                            color,
                                            strokeSize = strokeSize,
                                            pointSize = pointSize,
                                            strokeColor = strokeColor,
                                            alpha = alpha){
                pcs = pca_x$sing.val$u[,1:nPCs]
                pc.var = pca_x$var
                pair.pcs <- utils::combn(ncol(pcs), 2)
                pList <- list()
                for(i in 1:ncol(pair.pcs)){
                    if(i == 1){
                        x <- pair.pcs[1,i]
                        y <- pair.pcs[2,i]
                        p <- ggplot(mapping = aes(
                            x = pcs[,x],
                            y = pcs[,y],
                            fill = cat_var)) +
                            xlab(paste0('PC', x, ' (', pc.var[x], '%)')) +
                            ylab(paste0('PC', y, ' (', pc.var[y], '%)')) +
                            geom_point(
                                aes(fill = cat_var),
                                pch = 21,
                                color = strokeColor,
                                stroke = strokeSize,
                                size = pointSize,
                                alpha = alpha) +
                            scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            theme(
                                legend.position = "right",
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
                                aspect.ratio=1) +
                            guides(fill = guide_legend(override.aes = list(size = 4))) +
                            scale_fill_manual(name = cat_var_label, values = color)

                        le <- get_legend(p)
                    }else{
                        x <- pair.pcs[1,i]
                        y <- pair.pcs[2,i]
                        p <- ggplot(mapping = aes(
                            x = pcs[,x],
                            y = pcs[,y],
                            fill = cat_var)) +
                            xlab(paste0('PC', x, ' (',pc.var[x],  '%)')) +
                            ylab(paste0('PC', y, ' (',pc.var[y], '%)')) +
                            geom_point(
                                aes(fill = cat_var),
                                pch = 21,
                                color = strokeColor,
                                stroke = strokeSize,
                                size = pointSize,
                                alpha = alpha) +
                            scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
                            theme(
                                panel.background = element_blank(),
                                axis.line = element_line(colour = "black", size = 1.1),
                                legend.position = "none",
                                axis.text.x = element_text(size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.x = element_text(size = 14),
                                axis.title.y = element_text(size = 14),
                                aspect.ratio=1) +
                            scale_fill_manual(values = color, name = cat_var_label)
                    }
                    p <- p + theme(legend.position = "none")
                    xdens <- axis_canvas(p, axis = "x")+
                        geom_density(
                            mapping = aes(
                                x = pcs[,x],
                                fill = cat_var),
                            alpha = 0.7,
                            size = 0.2
                        ) +
                        theme(legend.position = "none") +
                        scale_fill_manual(values = color)

                    ydens <- axis_canvas(
                        p,
                        axis = "y",
                        coord_flip = TRUE) +
                        geom_density(
                            mapping = aes(
                                x = pcs[,y],
                                fill = cat_var),
                            alpha = 0.7,
                            size = 0.2) +
                        theme(legend.position = "none") +
                        scale_fill_manual(name = cat_var_label, values = color) +
                        coord_flip()

                    p1 <- insert_xaxis_grob(
                        p,
                        xdens,
                        grid::unit(.2, "null"),
                        position = "top"
                    )
                    p2 <- insert_yaxis_grob(
                        p1,
                        ydens,
                        grid::unit(.2, "null"),
                        position = "right"
                    )
                    pList[[i]] <- ggdraw(p2)
                }
                pList[[i+1]] <- le
                return(pList)
            }
            pca_x <- pca[[x]]
            p1=plot_pca_single_assay(pca_x,cat_var,cat_var_label,color,strokeSize,pointSize,strokeColor,alpha)
            p1
        })
        names(ppca) <- normalizations

        ## Prepare plot
        p=ppca[[1]]
        if (length(normalizations)>1){
            for (n in 2:length(normalizations)){
                p=c(p,ppca[[n]])
            }
        }
        plot=do.call(grid.arrange,
            c(p,
              ncol = ncol_plot,
           top=paste0("PCA ordered as ", paste(normalizations, collapse = ","))))

        return(plot=as_ggplot(plot))
}

