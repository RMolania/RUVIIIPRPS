#' is used to find repeating factors ro characters in a vector.
#'
#' @param se.obj A vector of factors or characters.
#' @param variables Numeric. Indicates the minimum repeat of individual factors ro characters in the vector.
#'
#'
#' @return A vec of factors ro characters that are repeated at least "n.repeat" times.

getAssessmentMetrics <- function(
        se.obj,
        variables) {
    categorical.var <- continuous.var <- NULL
    if (!is.null(variables)) {
        var.class <- sapply(
            variables,
            function(x)
                class(colData(se.obj)[[x]]))
        categorical.var <- names(var.class[var.class %in% c('character', 'factor')])
        continuous.var <- names(var.class[var.class %in% c('numeric', 'integer')])
    }

    # general plots #####
    general.plot <- 'General||boxPlot||RLE'

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
            metrics.for.cont.var
        )
        metrics.for.cont.var <- metrics.for.cont.var[order(metrics.for.cont.var$Var1), ]
        metrics.for.cont.var <- unname(unlist(apply(
            metrics.for.cont.var,
            1,
            function(x) paste0(x, collapse = '_'))))

    } else metrics.for.cont.var <- NULL

    # possible metrics for continuous variables #####
    if(length(categorical.var) > 0){
        metrics.for.cat.var <- c(
            '||coloredRLEplot||RLE',
            'rleMedians||boxPlot||RLE',
            'rleIqr||boxPlot||RLE',
            'pcs||boxPlot||PCA',
            'pcs||scatterPlot||PCA',
            'pcs||lineDotPlot||PcaVecCorr',
            '||barPlot||ARI',
            '||barPlot||Silhouette',
            'geneAnov||boxPlot||GeneVarAov',
            '||pvalHist||DGE'
        )
        metrics.for.cat.var <- expand.grid(
            categorical.var,
            metrics.for.cat.var
            )
        metrics.for.cat.var <- metrics.for.cat.var[order(metrics.for.cat.var$Var1) , ]
        metrics.for.cat.var <- unname(unlist(apply(
            metrics.for.cat.var,
            1,
            function(x) paste0(x, collapse = '_'))))
        metrics.for.cat.var <- sub("_\\|\\|", "||", metrics.for.cat.var)
    } else metrics.for.cat.var <- NULL

    # combined ari and silhouette
    if(length(categorical.var) > 1){
        metrics.for.two.cat.var <- c(
            '||combinedPlot||ARI',
            '||combinedPlot||Silhouette')
        two.cat.vars <- apply(
            combn(x = categorical.var, m = 2),
            2,
            function(x) paste0(x, collapse = '_'))

        metrics.for.two.cat.var <- expand.grid(
            two.cat.vars,
            metrics.for.two.cat.var)
        metrics.for.two.cat.var <- unlist(apply(
            metrics.for.two.cat.var,
            1,
            function(x) paste0(x, collapse = '')))
    } else metrics.for.two.cat.var <- NULL
    final.metrics <- c(
        general.plot,
        metrics.for.cont.var,
        metrics.for.cat.var,
        metrics.for.two.cat.var)

    return(unlist(final.metrics))
}
