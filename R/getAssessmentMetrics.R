#' Create all possible assessment metrics for the variables.

#' @author Ramyar Molania

#' @description
#' This functions provides the names of all possible assessment metrics for the given variable(s). The list will be used
#' in the 'assessVariation' and 'assessNormalization' functions.

#' @param se.obj A summarized experiment object.
#' @param variables Symbols. A symbol and a vector of symbols indicating the columns names of variables in the samples
#' annotation in the SummarizedExperiment object. The 'variables' can be categorical and continuous.
#' @param output.file Symbol. A name for the output file.
#' @param plot.output Logical. Whether to print the plot or not.
#'
#' @return A list of all possible assessment metrics for the variables.

#' @importFrom SummarizedExperiment colData
#' @importFrom tibble tibble
#' @importFrom igraph vertex_attr layout_as_tree graph_from_data_frame
#' @importFrom ggpubr ggarrange


getAssessmentMetrics <- function(
        se.obj,
        variables,
        output.file,
        plot.output = TRUE) {
    categorical.var <- continuous.var <- NULL
    if (!is.null(variables)) {
        var.class <- sapply(
            variables,
            function(x)
                class(colData(se.obj)[[x]]))
        categorical.var <- names(var.class[var.class %in% c('character', 'factor')])
        continuous.var <- names(var.class[var.class %in% c('numeric', 'integer')])
    }

    cat.char <- lapply(
        categorical.var,
        function(x){
            all.char <- gregexpr(pattern ='_', text = x)[[1]]
            all.char[1:length(all.char)]
        })
    categorical.var <- gsub('_', '.', categorical.var)
    names(cat.char) <- categorical.var
    cont.char <- lapply(
        continuous.var,
        function(x){
            all.char <- gregexpr(pattern ='_', text = x)[[1]]
            all.char[1:length(all.char)]
        } )
    continuous.var <- gsub('_', '.', continuous.var)
    names(cont.char) <- continuous.var
    all.var.char <- c(cat.char, cont.char)

    # general plots #####
    general.plot <- data.frame(
        Variables = 'General',
        Factors = 'General',
        PlotTypes = 'RLEplot',
        Metrics = 'RLE')

    # possible metrics for continuous variables #####
    if(length(continuous.var) > 0){
        metrics.for.cont.var <- c(
            'rleMedians||scatterPlot||RLE',
            'rleIqr||scatterPlot||RLE',
            'pcs||scatterPlot||PCA',
            'pcs||lineDotPlot||PcaReg',
            'geneCorr||boxPlot||GeneVarCorr')
        metrics.for.cont.var <- expand.grid(
            continuous.var,
            metrics.for.cont.var)
        colnames(metrics.for.cont.var) <- c('Variables', 'MetricsPlots')
        metrics.for.cont.var <- metrics.for.cont.var[order(metrics.for.cont.var$Variables), ]
        metrics.for.cont.var$Metrics <- unlist(lapply(
            metrics.for.cont.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][3]))
        metrics.for.cont.var$PlotTypes <- unlist(lapply(
            metrics.for.cont.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][2]))
        metrics.for.cont.var$Factors <- unlist(lapply(
            metrics.for.cont.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][1]))
        metrics.for.cont.var.table <- metrics.for.cont.var[,c(1,5,4,3)]
        metrics.for.cont.var.list <- paste0(
            metrics.for.cont.var$Variables,
            '_',
            metrics.for.cont.var$MetricsPlots)

    } else {
        metrics.for.cont.var.table <- NULL
        metrics.for.cont.var.list <- NULL
    }

    # possible metrics for continuous variables #####
    if(length(categorical.var) > 0){
        metrics.for.cat.var <- c(
            'rle||coloredRLEplot||RLE',
            'rleMedians||boxPlot||RLE',
            'rleIqr||boxPlot||RLE',
            'pcs||boxPlot||PCA',
            'pcs||scatterPlot||PCA',
            'pcs||lineDotPlot||PcaVecCorr',
            'ariCoeff||barPlot||ARI',
            'silhouetteCoeff||barPlot||Silhouette',
            'geneAnov||boxPlot||GeneVarAov',
            'pvaluse||pvalHist||DGE')
        metrics.for.cat.var <- expand.grid(
            categorical.var,
            metrics.for.cat.var)
        colnames(metrics.for.cat.var) <- c('Variables', 'MetricsPlots')
        metrics.for.cat.var <- metrics.for.cat.var[order(metrics.for.cat.var$Variables), ]
        metrics.for.cat.var$Metrics <- unlist(lapply(
            metrics.for.cat.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][3]))
        metrics.for.cat.var$PlotTypes <- unlist(lapply(
            metrics.for.cat.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][2]))
        metrics.for.cat.var$Factors <- unlist(lapply(
            metrics.for.cat.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][1]))
        metrics.for.cat.var.table <- metrics.for.cat.var[,c(1,5,4,3)]
        metrics.for.cat.var.list <- paste0(
            metrics.for.cat.var$Variables,
            '_',
            metrics.for.cat.var$MetricsPlots)

    } else {
        metrics.for.cat.var.table <- NULL
        metrics.for.cat.var.list <- NULL
        }

    # combined ari and silhouette
    if(length(categorical.var) > 1){
        metrics.for.two.cat.var <- c(
            'ariCoeff||combinedPlot||ARI',
            'silhouetteCoeff||combinedPlot||Silhouette')
        categorical.var <- names(var.class[var.class %in% c('character', 'factor')])
        two.cat.vars <- apply(
            combn(x = categorical.var, m = 2),
            2,
            function(x) paste0(x, collapse = '&'))

        metrics.for.two.cat.var <- expand.grid(
            two.cat.vars,
            metrics.for.two.cat.var)
        colnames(metrics.for.two.cat.var) <- c('Variables', 'MetricsPlots')
        metrics.for.two.cat.var <- metrics.for.two.cat.var[order(metrics.for.two.cat.var$Variables), ]
        metrics.for.two.cat.var$Metrics <- unlist(lapply(
            metrics.for.two.cat.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][3]))
        metrics.for.two.cat.var$PlotTypes <- unlist(lapply(
            metrics.for.two.cat.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][2]))
        metrics.for.two.cat.var$Factors <- unlist(lapply(
            metrics.for.two.cat.var$MetricsPlots,
            function(x) strsplit(x = as.character(x), split = '\\|\\|')[[1]][1]))
        metrics.for.two.cat.var.table <- metrics.for.two.cat.var[,c(1,5,4,3)]
        metrics.for.two.cat.var.list <- paste0(
            metrics.for.two.cat.var$Variables,
            '_',
            metrics.for.two.cat.var$MetricsPlots)
    } else {
        metrics.for.two.cat.var.table <- NULL
        metrics.for.two.cat.var.list <- NULL
    }

    # final metric files ####
    final.metrics.list <- c(
        metrics.for.cont.var.list,
        metrics.for.cat.var.list,
        metrics.for.two.cat.var.list,
        'GeneralRLEplot')
    final.metrics.table <- as.data.frame(rbind(
        metrics.for.cont.var.table,
        metrics.for.cat.var.table,
        metrics.for.two.cat.var.table))
    final.metrics.table$Variables <- as.character(final.metrics.table$Variables)
    rownames(final.metrics.table) <- c(1:nrow(final.metrics.table))
    all.var.char
    for(i in names(all.var.char)){
        if(all.var.char[[i]][1] != -1){
            cur.name <- names(all.var.char[i])
            for(j in 1:length(all.var.char[[i]]) ){
                pos.to.rep <- all.var.char[[i]][j]
                init.name <- paste0(
                    substring(cur.name, 1,  pos.to.rep - 1),
                    '_',
                    substring(cur.name, pos.to.rep + 1))
                cur.name <- init.name
                }
            final.metrics.table$Variables[final.metrics.table$Variables == names(all.var.char[i])] <- cur.name
        }
    }
    final.metrics.table.toplot <- final.metrics.table
    final.metrics.table <- rbind(final.metrics.table , general.plot)

    # plot ####
    plot.metrics <- lapply(
        variables,
        function(x){
            sub.final.metrics.table.toplot <- final.metrics.table.toplot[final.metrics.table.toplot$Variables == x, ]
            metrics.tree <- tibble::tibble(
                from = sub.final.metrics.table.toplot$Variables,
                to = paste(
                    paste0('T: ', sub.final.metrics.table.toplot$Metrics),
                    paste0('V: ', sub.final.metrics.table.toplot$Factors),
                    paste0('P: ', sub.final.metrics.table.toplot$PlotTypes),
                    sep = '\n'))
            g <- igraph::graph_from_data_frame(metrics.tree, directed = TRUE)
            coords <- igraph::layout_as_tree(g)
            colnames(coords) <- c("x", "y")
            output.df <- tibble::as_tibble(coords) %>%
                mutate(
                    step = igraph::vertex_attr(g, "name"),
                    label = gsub("\\d+$", "", step),
                    x = x * -1,
                    type = factor(1))
            plot.nodes = output.df %>%
                mutate(
                    xmin = x - 0.45,
                    xmax = x + 0.45,
                    ymin = y - 0.1,
                    ymax = y + 0.1)
            plot.edges <- metrics.tree %>%
                dplyr::mutate(id = dplyr::row_number()) %>%
                pivot_longer(
                    cols = c("from", "to"),
                    names_to = "s_e",
                    values_to = "step") %>%
                dplyr::left_join(plot.nodes, by = "step") %>%
                dplyr::select(-c(label, type, y, xmin, xmax)) %>%
                dplyr::mutate(y = ifelse(s_e == "from", ymin, ymax)) %>%
                dplyr::select(-c(ymin, ymax))
            plot.nodes$xmin[1] <- 0
            plot.nodes$ymin[1] <- 0
            plot.nodes$xmax[1] <- 0
            plot.nodes$ymax[1] <- 0
            p <- ggplot() + geom_rect(
                data = plot.nodes,
                mapping = aes(
                    xmin = xmin,
                    ymin = ymin,
                    xmax = xmax,
                    ymax = ymax,
                    fill = type,
                    colour = type),
                alpha = 0.5,
                color = 'grey',
                fill = 'white')
            p <- p + geom_text(
                data = plot.nodes[1,],
                mapping = aes(x = x, y = y, label = label),
                color = "black", size = 10
            )
            p <- p + geom_text(
                data = plot.nodes[-1,],
                mapping = aes(x = x, y = y, label = label),
                color = "black",angle = 30
            )
            p <- p + geom_path(
                data = plot.edges,
                mapping = aes(x = x, y = y, group = id),
                colour = "#585c45",
                arrow = arrow(length = unit(0.3, "cm"), type = "closed")
            )
            p <- p + theme(
                panel.background = element_blank(),
                axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = 'none')
        })

    pdf(paste0(output.file, '.pdf'), width = 8*length(variables), height = 14)
    plot.caption <- expression(atop(
        scriptstyle("T: test | V: variable | P: plot type"))
        )
    all.plots <- ggarrange(plotlist = plot.metrics, common.legend = TRUE)
    all.plots <- ggpubr::annotate_figure(
        p = all.plots,
        bottom = ggpubr::text_grob(label = plot.caption, size = 20))
    print(all.plots)
    dev.off()
    if(isTRUE(plot.output)) print(all.plots)
    return(all.metrics = list(
        final.metrics.list = unlist(final.metrics.list),
        final.metrics.table = final.metrics.table))
}
