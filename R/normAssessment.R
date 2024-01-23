#' is used to assess the performance of the normalisation of a SummarizedExperiment class object.
#'
#' Several assessment will be performed:
#' For each categorical variable:
#' - PCA plot of the categorical variable.
#' - Silhouette and ARI computed on the categorical variable.
#' - Differential analysis based ANOVA between the gene expression and the categorical variable.
#' - Vector correlation between the first cumulative PCs of the gene expression and the categorical variable.
#' For each continous variable:
#' - Linear regression between the first cumulative PC and continuous variable.
#' - Correlation between gene expression and continuous variable.
#'
#' It will output the following plots:
#' - PCA plot of each categorical variable.
#' - Boxplot of the F-test distribution from ANOVA between the gene expression and each categorical variable.
#' - Vector correlation between the first cumulative PCs of the gene expression and each categorical variable.
#' - Combined Silhouette plot of the combined pair of all categorical variables.
#' - Linear regression between the first cumulative PC and continuous variable.
#' - Boxplot of the correlation between gene expression and continuous variable.
#' - It will also output the RLE plot distribution.
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variables String or vector of strings of the label of continuous or categorical variable(s)
#' such as samples types, batch or library size from colData(se).
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#' @param output.file Path and name of the output file to save the assessments plots in a pdf format.
#' @param fast.pca Logical. Indicates whether to calculate a specific number of PCs instead of the full range
#' to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to be calculated for the fast pca process, by default is set to 10.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param metrics Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param metrics.to.exclude TTTT
#' @param center.pca TTT
#' @param scale.pca TTTT
#' @param bsparam TTTT
#' @param remove.na TTTT
#' @param save.se.obj TTTT
#' @param silhouette.dist.measure TTTT
#' @param plot.ncol TTTT
#' @param rle.outputs.to.return TTT
#' @param ylim.rle.plot TTT
#' @param rle.median.points.size TTT
#' @param rle.median.points.color TTT
#' @param rle.iqr.width TTT
#' @param rle.geom.hline.color TTT
#' @param rle.plot.ncol TTT
#' @param compute.nb.pcs TTT
#'
#' @return  SummarizedExperiment A SummarizedExperiment object containing all the assessments plots and metrics.
#' If specified it will generate a pdf containing the assessments plots and metrics used for the assessment.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom kunstomverse geom_boxplot2
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom gridExtra grid.arrange
#' @importFrom SummarizedExperiment assays colData
#' @export

## remove n.core and merge all variable together
normAssessment <- function(
        se.obj,
        assay.names = 'all',
        variables,
        metrics = 'all',
        metrics.to.exclude = 'NULL',
        rle.outputs.to.return = 'all',
        ylim.rle.plot = c(-2,2),
        rle.median.points.size = 1,
        rle.median.points.color = "red",
        rle.iqr.width = 2,
        rle.geom.hline.color = "cyan",
        rle.plot.ncol = 1,
        fast.pca = TRUE,
        compute.nb.pcs = 10,
        center.pca = TRUE,
        scale.pca = FALSE,
        bsparam = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        plot.ncol = 1,
        silhouette.dist.measure = NULL,
        assess.se.obj = FALSE,
        remove.na = 'none',
        save.se.obj = TRUE,
        output.file = NULL,
        verbose = TRUE
){
    printColoredMessage(message = '------------The normAssessment function starts:',
                        color = 'white',
                        verbose = verbose)
    # check the inputs of PCA ####
    if (length(assay.names) == 1 && assay.names != 'All') {
        if (!assay.names %in% names(assays(se.obj)))
            stop('The assay name cannot be found in the SummarizedExperiment object.')
    }
    if (length(assay.names) > 1) {
        if (!assay.names %in% names(assays(se.obj)))
            stop('The assay names cannot be found in the SummarizedExperiment object.')
    }
    if (fast.pca & is.null(nb.pcs)) {
        stop('To perform fast PCA, the number of PCs (nb.pcs) must specified.')
    } else if (fast.pca & nb.pcs == 0) {
        stop('To perform fast PCA, the number of PCs (nb.pcs) must specified.')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- as.factor(unlist(assay.names))

    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = NULL,
            remove.na = 'measurements',
            verbose = verbose
        )
    }

    # find categorical and continuous variables ####
    categorical.var <- continuous.var <- NULL
    if (!is.null(variables)) {
        var.class <- sapply(variables,
                            function(x)
                                class(colData(se.obj)[[x]]))
        categorical.var <- names(var.class[var.class %in% c('character', 'factor')])
        continuous.var <- names(var.class[var.class %in% c('numeric', 'integer')])
    }
    # all possible metrics for each variable #####
    all.metrics <- getAssessmentMetrics(
        se.obj = se.obj,
        variables = variables)

    # metrics #####
    metrics.to.compute <- unlist(lapply(
        all.metrics,
        function(x){
            le <- length(strsplit(x , '\\|\\|')[[1]])
            strsplit(x , '\\|\\|')[[1]][le]}))
    plots.to.generate <- unlist(lapply(
        all.metrics,
        function(x) strsplit(x , '\\|\\|')[[1]][2]))

    ## RLE #####
    ### compute rle #####
    if('RLE' %in% metrics.to.compute){
        se.obj <- computeRLE(
            se.obj = se.obj,
            assay.names = assay.names,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            outputs.to.return = rle.outputs.to.return,
            assess.se.obj = assess.se.obj,
            remove.na = remove.na,
            save.se.obj = TRUE,
            verbose = verbose)}

    ### plot general rle #####
    if('General||boxPlot||RLE' %in% all.metrics ){
        se.obj <- plotRLE(
            se.obj = se.obj,
            assay.names = assay.names,
            variable = NULL,
            ylim.rle.plot = ylim.rle.plot,
            median.points.size = 1,
            median.points.color = "red",
            iqr.width = 2,
            geom.hline.color = "cyan",
            plot.ncol = 1,
            plot.output = FALSE,
            save.se.obj = TRUE,
            verbose = TRUE)
    }

    ### plot colored rle #####
    if('coloredRLEplot' %in% plots.to.generate){
        to.plot <- all.metrics[plots.to.generate == 'coloredRLEplot']
        vars <- unlist(lapply(
            to.plot,
            function(x) strsplit(x = x, split = '\\|\\|')[[1]][1]))
        for(i in vars){
            se.obj <- plotRLE(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = i,
                ylim.rle.plot = c(-2, 2),
                median.points.size = .5,
                median.points.color = "black",
                iqr.width = 2,
                geom.hline.color = "cyan",
                plot.ncol = 1,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = TRUE)
        }
    }
    ### plot rle medinas\iqr with variable #####
    rle.var.plots <- c(
        all.metrics[plots.to.generate == 'scatterPlot' & metrics.to.compute == 'RLE'],
        all.metrics[plots.to.generate == 'boxPlot' & metrics.to.compute == 'RLE']
    )
    rle.var.plots <- rle.var.plots[!rle.var.plots %in% "General||boxPlot||RLE"]
    if(length(rle.var.plots) > 0){
        vars <- unlist(lapply(
            rle.var.plots,
            function(x) {
                x <- strsplit(x = x, split = '\\|\\|')[[1]][1]
                strsplit(x = x, split = '_')[[1]][1]
            }))
        for(i in vars){
            se.obj <- plotRleVariable(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = i,
                rle.data.type = 'both',
                ylim.rle.med.plot = NULL,
                ylim.rle.iqr.plot = NULL,
                points.size = 1,
                plot.ncol = 1,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }

    # PCA ####
    ## compute pca ####
    if('PCA' %in% metrics.to.compute){
        se.obj <- computePCA(
            se.obj = se.obj,
            assay.names = assay.names,
            fast.pca = fast.pca,
            nb.pcs = nb.pcs,
            scale = scale.pca,
            center = center.pca,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            bsparam = NULL,
            assess.se.obj = assess.se.obj,
            remove.na = remove.na,
            save.se.obj = save.se.obj,
            verbose = verbose)
    }
    ## plot pca ####
    pca.var.plots <- c(
        all.metrics[plots.to.generate == 'scatterPlot' & metrics.to.compute == 'PCA'],
        all.metrics[plots.to.generate == 'boxPlot' & metrics.to.compute == 'PCA']
    )
    if(length(pca.var.plots) > 0){
        for(i in pca.var.plots){
            var.meric <- strsplit(i, '\\|\\|')
            var <- strsplit(var.meric[[1]], split = '_')[[1]][1]
            if(var.meric[[1]][2] == 'scatterPlot'){
                plot.type.pca <- 'scatter'
            } else plot.type.pca <- 'boxplot'
                se.obj <- plotPCA(
                    se.obj = se.obj,
                    assay.names = assay.names,
                    variable = var,
                    plot.type = plot.type.pca,
                    fast.pca = fast.pca,
                    nb.pcs = 3,
                    save.se.obj = TRUE,
                    verbose = TRUE)
            }
        }

    # Vector correlation ####
    ## compute and plot vector correlation ####
    if('PcaVecCorr' %in% metrics.to.compute){
        var.metric <- all.metrics[metrics.to.compute == 'PcaVecCorr']
        for(i in var.metric){
            var <- strsplit(x = i, split = '\\|\\|')[[1]][1]
            var <- strsplit(x = var, split = '_')[[1]][1]
            se.obj <- computePCVariableCorrelation(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = var,
                fast.pca = fast.pca,
                nb.pcs = nb.pcs,
                save.se.obj = save.se.obj,
                verbose = verbose)
            se.obj <- plotPCVariableCorrelation(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = var,
                fast.pca = fast.pca,
                nb.pcs = nb.pcs,
                plot.output = FALSE,
                save.se.obj = save.se.obj,
                verbose = verbose)
        }
    }

    # Linear regression ####
    ## compute and plot linear regression ####
    if('PcaReg' %in% metrics.to.compute){
        var.metric <- all.metrics[metrics.to.compute == 'PcaReg']
        for(i in var.metric){
            var <- strsplit(x = i, split = '\\|\\|')[[1]][1]
            var <- strsplit(x = var, split = '_')[[1]][1]
            se.obj <- computePCVariableRegression(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = var,
                fast.pca = fast.pca,
                nb.pcs = nb.pcs,
                save.se.obj = save.se.obj,
                verbose = verbose)
            se.obj <- plotPCVariableRegression(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = var,
                fast.pca = fast.pca,
                nb.pcs = nb.pcs,
                plot.output = FALSE,
                save.se.obj = save.se.obj,
                verbose = verbose)
        }
    }

    # Silhouette coefficient ####
    ## compute and plot adjusted Silhouette coefficient ####
    if('Silhouette' %in% metrics.to.compute){
        var.metric <- all.metrics[metrics.to.compute == 'Silhouette']
        vars <- unique(unlist(lapply(
            var.metric,
            function(x) {
                vars <- strsplit(x = x, split = '\\|\\|')[[1]][1]
                vars <- unlist(strsplit(x = vars, '_'))})))
        for(i in vars){
            se.obj <- computeSilhouette(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = i,
                dist.measure = silhouette.dist.measure,
                fast.pca = fast.pca,
                nb.pcs = nb.pcs,
                save.se.obj = save.se.obj,
                verbose = verbose)
        }
        for(i in var.metric){
            plot.type <- strsplit(x = i, split = '\\|\\|')[[1]][2]
            if(plot.type == 'barPlot'){
                var <- strsplit(x = i, split = '\\|\\|')[[1]][1]
                plot.type <- 'single.plot'
            } else{
                var <- strsplit(x = i, split = '\\|\\|')[[1]][1]
                var <- unlist(strsplit(x = var, split = '_'))
                plot.type <- 'combined.plot'
            }
            se.obj <- plotSilhouette(
                se.obj = se.obj,
                assay.names = assay.names,
                variables = var,
                plot.type = plot.type,
                silhouette.method = 'sil.euclidian',
                plot.output = FALSE,
                save.se.obj = save.se.obj,
                verbose = verbose)
        }
    }

    # ARI ####
    ## compute and plot adjusted rand index ####
    if('ARI' %in% metrics.to.compute){
        var.metric <- all.metrics[metrics.to.compute == 'ARI']
        vars <- unique(unlist(lapply(
            var.metric,
            function(x) {
                vars <- strsplit(x = x, split = '\\|\\|')[[1]][1]
                vars <- unlist(strsplit(x = vars, '_'))})))
        for(i in vars){
            se.obj <- computeARI(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = i,
                clustering.method = 'hclust',
                hclust.method = 'complete',
                hclust.dist.measure = 'euclidian',
                fast.pca = fast.pca,
                nb.pcs = nb.pcs,
                save.se.obj = save.se.obj,
                verbose = verbose)
        }
        for(i in var.metric){
            plot.type <- strsplit(x = i, split = '\\|\\|')[[1]][2]
            if(plot.type == 'barPlot'){
                var <- strsplit(x = i, split = '\\|\\|')[[1]][1]
                plot.type <- 'single.plot'
            } else{
                var <- strsplit(x = i, split = '\\|\\|')[[1]][1]
                var <- unlist(strsplit(x = var, split = '_'))
                plot.type <- 'combined.plot'
            }
            se.obj <- plotARI(
                se.obj = se.obj,
                assay.names = assay.names,
                variables = var,
                plot.type = plot.type,
                ari.method = 'hclust.complete.euclidian',
                plot.output = FALSE,
                save.se.obj = save.se.obj,
                verbose = verbose)
        }
    }
    # Gene variable correlation ####
    ## compute and plot gene variable correlation ####
    if('GeneVarAov' %in% metrics.to.compute){
        var.metric <- all.metrics[metrics.to.compute == 'GeneVarCorr']
        vars <- unique(unlist(lapply(
            var.metric,
            function(x) {
                vars <- strsplit(x = x, split = '\\|\\|')[[1]][1]
                vars <- unlist(strsplit(x = vars, '_')[[1]][1])})))
        for(i in vars){
            se.obj <- computeGenesVariableCorrelation(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = i,
                method = 'spearman',
                a = 0.05,
                rho = 0,
                plot.top.genes = FALSE,
                nb.top.genes = NULL,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                apply.round = TRUE,
                assess.se.obj = assess.se.obj,
                remove.na = remove.na,
                save.se.obj = save.se.obj)
            se.obj <- plotGenesVariableCorrelation(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = i,
                correlation.method = 'gene.spearman.corr',
                plot.output = FALSE,
                save.se.obj = save.se.obj,
                verbose = verbose)
        }
    }

    # Gene variable ANOVA ####
    ## compute and plot gene variable ANOVA ####
    if('GeneVarAov' %in% metrics.to.compute){
        var.metric <- all.metrics[metrics.to.compute == 'GeneVarAov']
        vars <- unique(unlist(lapply(
            var.metric,
            function(x) {
                vars <- strsplit(x = x, split = '\\|\\|')[[1]][1]
                vars <- unlist(strsplit(x = vars, '_')[[1]][1])})))
        for(i in vars){
            se.obj <- computeGenesVariableAnova(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = i,
                method = 'aov',
                plot.output = FALSE,
                plot.top.genes = FALSE,
                nb.top.genes = NULL,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                apply.round = TRUE,
                assess.se.obj = assess.se.obj,
                remove.na = remove.na,
                save.se.obj = save.se.obj)
            se.obj <- plotGenesVariableAnova(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = i,
                anova.method = "genes.aov.anova",
                plot.output = FALSE,
                save.se.obj = save.se.obj,
                verbose = verbose)
        }
    }

    # DGE ####
    ## compute and plot gene variable ANOVA ####
    if('DGE' %in% metrics.to.compute){
        var.metric <- all.metrics[metrics.to.compute == 'DGE']
        vars <- unique(unlist(lapply(
            var.metric,
            function(x) strsplit(x = x, split = '\\|\\|')[[1]][1])))
        for(i in vars){
            se.obj <- computeDGE(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = i,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                assess.se.obj = assess.se.obj,
                remove.na = remove.na,
                save.se.obj = save.se.obj,
                verbose = verbose)
            se.obj <- plotDGE(
                se.obj = se.obj,
                assay.names = assay.names,
                variable = i,
                plot.ncol = plot.ncol,
                plot.output = FALSE,
                save.se.obj = save.se.obj,
                verbose = verbose)
        }
    }

    # save all plots ####
    # all.vars <- unique(unlist(lapply(
    #     all.metrics,
    #     function(x){
    #         vars <- strsplit(x , '\\|\\|')[[1]][1]
    #         strsplit(vars , '_')[[1]][1]})))
    #
    # pdf_file <- "output.pdf"
    # pdf(pdf_file)
    # plot.new()
    # text(.5, .5, "Hello, this is some text in the PDF.", font = 2, cex = 1.5)
    # if("General||boxPlot||RLE" %in% all.metrics){
    #     print(se.obj@metadata$plot$RLE$GeneralRLE)
    # }
    # for(i in all.vars){
    #     if(colData(se.obj)[[i]] %in% c('numeric', 'integer')){
    #         if(i %in% se.obj@metadata$plot$PCA$fastPCA)
    #     }
    # }

    dev.off()
    printColoredMessage(message = '------------The normAssessment function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}
