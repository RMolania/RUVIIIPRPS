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
    final.metrics.table <- rbind(final.metrics.table , general.plot)
    return(all.metrics = list(
        final.metrics.list = unlist(final.metrics.list),
        final.metrics.table = final.metrics.table))
}
