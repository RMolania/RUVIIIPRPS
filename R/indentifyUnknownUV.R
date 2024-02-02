#' find potential unknown sources of unwanted variation in RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function uses three different approaches: 'rle', 'pca' and 'sample.scoring' to find potential unknown sources of
#' unwanted variation in RNA-seq data. In the rle approach, a clustering method will be
#' applied on either medians or IQR or both of the rle data.In the absent on unwanted variation, we should not find any
#' clustering. In the pca approach, first, a principal component analysis on either a set of negative control genes or
#' all genes will be applied and then a clustering method will be used to clusters the first principal components to
#' find potential unknown sources of unwanted variation. In the sample.scoring approach, first,  all samples will be
#' scored against a set of signature genes related to unwanted variation, then  clustering method will be applied
#' on the scoring to find potential unknown sources of unwanted variation.

#' @references
#' Gandolfo L. C. & Speed, T. P., RLE plots: visualizing unwanted variation in high dimensional data. PLoS ONE, 2018.\
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. A symbol for the selection of the name of the assay in the SummarizedExperiment object to
#' indentify possible unknown sources of unwanted variation.
#' @param approach Symbol. A symbol indicating the approach to be employed for identifying unknown sources of unwanted
#' variation. This should be one of 'rle', 'pca', or 'sample.scoring'.In the rle approach, a clustering method is applied
#' to either medians or IQR or both of the rle data. In the absence of unwanted variation, no distinguishable clustering
#' should occur. In the pca approach, a principal component analysis is initially performed on either a set of negative
#' control genes specified in 'ncg' or all genes. Subsequently, a clustering method is employed to cluster the first
#' principal components and identify potential unknown sources of unwanted variation. In the 'sample.scoring' approach, all
#' samples are initially scored against a set of signature genes related to unwanted variation specified in 'uv.gene.sets'.
#' Then, a clustering method is applied to the scores to identify potential unknown sources of unwanted variation.
#' @param regress.out.bio.variables Symbol. A symbol or a vector of symbols specifying the column names of biological
#' variables in the sample annotation of the SummarizedExperiment object. These variables can be either categorical or
#' continuous variables. The default is 'NULL'.
#' @param regress.out.bio.gene.sets List. A list of biological gene signatures. Individual gene sets will be used to
#' score samples and then each scores will be regressed out from the data before identifying sources of unwanted variation.
#' The default is 'NULL'.
#' @param uv.gene.sets List. A list of gene signatures related any possible unwanted variables. Individual gene sets will
# 'be used to score samples and then each scores will be clustered to identify possible unknown batches. The default is NULL.
#' @param ncg Vector. Specifies a set of negative control genes to identify sources of unwanted variation. The default
#' is 'NULL'. If not 'NULL', the RLE or PCA and sample scoring will be conducted solely based on the negative control genes.
#' @param clustering.methods Symbol. Specifies the clustering method to be utilized for identifying potential unknown batches.
#' This should be one of the following: 'kmeans', 'cut', 'quantile', 'nbClust'. The default is 'nbClust'.
#' @param nbClust.diss dissimilarity matrix to be used. By default, diss=NULL, but if it is replaced by a dissimilarity
#' matrix, distance should be "NULL".
#' @param nbClust.distance the distance measure to be used to compute the dissimilarity matrix. This must be one of:
#' "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski" or "NULL".
#' By default, distance="euclidean". If the distance is "NULL", the dissimilarity matrix (diss) should be given by the
#' user. If distance is not "NULL", the dissimilarity matrix should be "NULL".
#' @param nbClust.min.nc minimal number of clusters, between 1 and (number of objects - 1)
#' @param nbClust.max.nc Numeric. Indicates maximal number of clusters, between 2 and (number of objects - 1),
#' greater or equal to min.nc. By default, max.nc = 15. Refer to the NbClsut function for more details.
#' @param nbClust.method Symbol. indicates the cluster method to be used.
#' This should be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid",
#' "kmeans". Refer to the NbClsut function for more details.
#' @param nbClust.index The index to be calculated. This should be one of : 'kl', 'ch', 'hartigan', 'ccc', 'scott',
#' 'marriot', 'trcovw', 'tracew', 'friedman', 'rubin', 'cindex', 'db', 'silhouette', 'duda', 'pseudot2', 'beale',
#' 'ratkowsky', 'ball', 'ptbiserial', 'gap', 'frey', 'mcclain', 'gamma', 'gplus', 'tau', 'dunn', 'hubert', 'sdindex',
#' 'dindex', 'sdbw', 'all' (all indices except GAP, Gamma, Gplus and Tau),'alllong' (all indices with Gap, Gamma, Gplus
#' and Tau included).
#' @param nbClust.alphaBeale significance value for Beale's index.
#' @param max.samples.per.batch Numeric. Indicates the maximum number of samples per cluster when the clustering.methods
#' is nbClust. The default is 0.1 (10%) of total samples in the SummarizedExperiment object.
#' @param nb.clusters Numeric. Specifies the number of clusters to be identified when "clustering.methods" is set to 'kmeans',
# 'cut', and 'quantile'. The default is 3.
#' @param rle.comp A symbol. Specifies which properties, either 'median' or 'iqr' or 'both', of the RLE data should be used for
#' clustering when the approach is set to 'RLE'. The default is 'median'.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to 'TRUE'.
#' The data must be in log transformation before comouting RLE or PCA.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param nb.pcs Numeric. A value determining the number of first principal components to be selected for clustering
#' when the approach is set to 'PCA'. The default is 2.
#' @param center Logical. Indicates whether to center the data or not before performing PCA. The default is 'TRUE'.
#' @param scale Logical. Indicates whether to scale the data or not before performing PCA. The default is 'FASLE'.
#' @param svd.bsparam A BiocParallelParam object specifying how parallelization should be performed. The default is bsparam().
#' We refer to the 'runSVD' function from the BiocSingular R package for more details.
#' @param assess.se.obj Logical. Whether to assess the SummarizedExperiment object or not. If 'TRUE', the function
#' 'checkSeobj' will be applied. The default is 'TRUE'.
#' @param remove.na A symbol. Indicates whether to remove NA or missing values from either the 'assays',
#' 'sample.annotation', 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be
#' excluded. If 'sample.annotation' is selected, the samples that contains NA or missing values for any 'bio.variables'
#' and 'uv.variables' will be excluded. By default, it is set to both'.
#' @param save.se.obj Logical. Indicates whether to save the results in the metadata of the SummarizedExperiment class.
#' Th default is 'TRUE', the results will be save to the 'metadata$UV$Unknown'.
#' @param plot.out TTTT
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom singscore rankGenes simpleScore
#' @importFrom BiocSingular bsparam runSVD
#' @importFrom NbClust NbClust
#' @importFrom stats as.formula
#' @export

identifyUnknownUV <- function(
        se.obj,
        assay.name,
        approach = 'rle',
        rle.comp = 'median',
        regress.out.bio.variables = NULL,
        regress.out.bio.gene.sets = NULL,
        uv.gene.sets = NULL,
        ncg = NULL,
        clustering.methods = 'nbClust',
        nbClust.diss = NULL,
        nbClust.distance = "euclidean",
        nbClust.min.nc = 2,
        nbClust.max.nc = 5,
        nbClust.method = 'kmeans',
        nbClust.index = 'silhouette',
        nbClust.alphaBeale = 0.1,
        max.samples.per.batch = .1,
        nb.clusters = 3,
        apply.log = TRUE,
        pseudo.count = 1,
        nb.pcs = 2,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        plot.out = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The indentifyUnknownUV function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if (is.null(assay.name)){
        stop('The "assay.name" cannot be empty.')
    } else if (length(assay.name) > 1) {
        stop('The "assay.name" must be a single assay name.')
    }
    if(is.null(approach)){
        stop('The "approach" cannot be empty.')
    } else if (length(approach) > 1){
        stop('The approach must be one of the "rle", "pca" or "sample.scoring".')
    } else if (!approach %in% c('rle', 'pca', 'sample.scoring') ) {
        stop('The approach must be one of the "rle", "pca" or "sample.scoring".')
    }
    if (!is.null(regress.out.bio.variables)) {
        if (!regress.out.bio.variables %in% colnames(colData(se.obj)))
            stop('The "regress.out.bio.variables" are not found in the SummarizedExperiment object.')
    }
    if(!is.null(regress.out.bio.gene.sets)){
        lapply(
            regress.out.bio.gene.sets,
            function(x){
                if(sum(regress.out.bio.gene.sets %in% row.names(se.obj)) == 0)
                    stop('The "regress.out.bio.gene.sets" are not found in the SummarizedExperiment object.')
            })
    }
    if(!is.null(uv.gene.sets)){
        lapply(
            names(uv.gene.sets),
            function(x){
                if(sum(uv.gene.sets[[x]] %in% row.names(se.obj)) < 2)
                    stop('The "uv.gene.sets" are not found in the SummarizedExperiment object.')
            })
    }
    if (!is.null(ncg)) {
        if (sum(ncg %in% row.names(se.obj)) == 0)
            stop('The "ncg" genes are found in the SummarizedExperiment object.')
    }
    if (length(clustering.methods) > 1) {
        stop('A single method should be provided for the clustering.methods.')
    }
    if(clustering.methods == 'nbClust'){
        if(is.null(nbClust.min.nc)){
            stop('The "nbClust.min.nc" must be specified, when the clustering.methods is nbClust.')
        } else if (nbClust.min.nc < 0 | nbClust.min.nc == 1){
            stop('The "nbClust.min.nc" must be equal or more than 2, when the clustering.methods is nbClust.')
        } else if(is.null(nbClust.max.nc)){
            stop('The "nbClust.max.nc" must be specified, when the clustering.methods is nbClust.')
        } else if (nbClust.max.nc < 0 | nbClust.max.nc < 2){
            stop('The "nbClust.max.nc" must be equal or more than 2, when the clustering.methods is nbClust.')
        } else if(!nbClust.method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans")){
            stop('The "nbClust.method" must be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans."')
        } else if (!nbClust.index %in% c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db",
                                         "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain",
                                         "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw",
                                         "all", "alllong")){
            stop('The nbClust.index must of one of: "kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db",
                                         "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain",
                                         "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw",
                                         "all", "alllong"')
        } else if (max.samples.per.batch == 0 | max.samples.per.batch < 0  | max.samples.per.batch >= 1){
            stop('The value of max.samples.per.batch must be between 0<max.samples.per.batch<1.')
        }
    }

    if (!clustering.methods %in% c('kmeans', 'cut', 'quantile', 'nbClust')) {
        stop('The clustering.methods should be one of kmeans, cut, quantile or nbClust.')
    }
    if(clustering.methods %in% c('kmeans', 'cut', 'quantile')){
        if(is.null(nb.clusters))
            stop('The "nb.clusters" must be specified when the "clustering.methods" is kmeans, cut or quantile.')
    }
    if(approach == 'rle'){
        if(!rle.comp %in% c('median', 'iqr', 'both'))
            stop('The "rle.comp" should be one of "median", "iqr" or "both".')
    }
    if(approach == 'pca'){
        if(nb.pcs == 0 | is.null(max.samples.per.batch))
            stop('The value of nb.pcs should be more than 0 when the arroach is equal to pca.')
    }
    if(apply.log){
        if (pseudo.count < 0)
            stop('The valuse of pseudo.count cannot be negative.')
    }
    if(approach == 'pca' & nb.pcs > 1 & clustering.methods %in% c('cut', 'quantile')){
        stop(paste0('The nb.pcs should be 1 to use the ', clustering.methods, ' method for clustering.'))
    }
    if(is.null(regress.out.bio.variables) & remove.na == 'both'){
        stop('The "remove.na" cannot be set to "both" when the "regress.out.bio.variables = NULL".')
    } else if (is.null(regress.out.bio.variables) & remove.na == 'sample.annotation'){
        stop('The "remove.na" cannot be set to "sample.annotation" when the "regress.out.bio.variables = NULL".')
    }

    # check the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(regress.out.bio.variables),
            remove.na = remove.na,
            verbose = verbose)
    }
    # data transformation ####
    printColoredMessage(message = '-- Data transformation:',
        color = 'magenta',
        verbose = verbose)
    if (isTRUE(apply.log) & !is.null(pseudo.count)) {
        printColoredMessage(
            message = paste0('Apply log2 + ', pseudo.count,  ' (pseudo.count) on the ', assay.name, ' data.'),
            color = 'blue',
            verbose = verbose)
        temp.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else if (isTRUE(apply.log) & is.null(pseudo.count)){
        printColoredMessage(
            message = paste0('Apply log2 on the ', assay.name, ' data.'),
            color = 'blue',
            verbose = verbose)
        temp.data <- log2(assay(x = se.obj, i = assay.name))
    } else {
        printColoredMessage(
            message = paste0('The ', assay.name, ' data will be used without log transformation.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = 'Please note, the assay should be in log scale before computing RLE or PCA.',
            color = 'red',
            verbose = verbose)
        temp.data <- assay(x = se.obj, i = assay.name)
    }

    # regress out bio variables and bio gene sets ####
    if(!is.null(regress.out.bio.variables) | !is.null(regress.out.bio.gene.sets)){
        printColoredMessage(
            message = '-- Regress out "regress.out.bio.variables" and "regress.out.bio.gene.sets" from the data:',
            color = 'magenta',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'We do not recommend regressing out any biological variation',
                'if they may be largely associated with the unwanted variation.'),
            color = 'red',
            verbose = verbose)

        # sample scoring for regress.out.bio.gene.sets ####
        if(!is.null(regress.out.bio.gene.sets)){
            printColoredMessage(
                message = '-Calculate sample scores for individual gene sets of the "regress.out.bio.gene.sets".',
                color = 'blue',
                verbose = verbose)
            ranked.data <- rankGenes(temp.data)
            regress.out.bio.gene.sets <- sapply(
                regress.out.bio.gene.sets,
                function(x) singscore::simpleScore(
                    rankData = ranked.data,
                    upSet = x)$TotalScore)
            rm(ranked.data)
        }
        if(!is.null(regress.out.bio.variables) & is.null(regress.out.bio.gene.sets)){
            ## regress out regress.out.bio.variables ####
            printColoredMessage(
                message = paste0(
                    'The ',
                    paste0(regress.out.bio.variables, collapse = ' & '),
                    ' variables will be regressed out from the data,',
                    ' please make sure your data is log transformed.'),
                color = 'blue',
                verbose = verbose)
            temp.data <- t(temp.data)
            lm.formula <- paste('se.obj', regress.out.bio.variables, sep = '$')
            adjusted.data <- lm(as.formula(paste('temp.data', paste0(lm.formula, collapse = '+') , sep = '~')))
            temp.data <- t(adjusted.data$residuals)
            colnames(temp.data) <- colnames(se.obj)
            row.names(temp.data) <- row.names(se.obj)
        } else if(is.null(regress.out.bio.variables) & !is.null(regress.out.bio.gene.sets)){
            ## regress out regress.out.bio.gene.sets ####
            printColoredMessage(
                message = paste0(
                    'The sample scores of individual gene list of ',
                    'regress.out.bio.gene.sets',
                    ' will be regressed out from the data,',
                    ' please make sure your data is log transformed.'),
                color = 'blue',
                verbose = verbose)
            temp.data <- t(temp.data)
            adjusted.data <- lm(temp.data~ regress.out.bio.gene.sets)
            temp.data <- t(adjusted.data$residuals)
            colnames(temp.data) <- colnames(se.obj)
            row.names(temp.data) <- row.names(se.obj)
        } else if (!is.null(regress.out.bio.variables) & !is.null(regress.out.bio.gene.sets)){
            printColoredMessage(
                message = paste0(
                    'The sample scores of individual gene list of ',
                    'regress.out.bio.gene.sets and the ',
                    paste0(regress.out.bio.variables, collapse = ' & '),
                    ' variables will be regressed out from the data,',
                    ' please make sure your data is log transformed.'),
                color = 'blue',
                verbose = verbose)
            all.variables <- as.data.frame(cbind(
                regress.out.bio.gene.sets,
                as.data.frame(colData(se.obj)[, regress.out.bio.variables, drop = FALSE])))
            lm.formula <- paste('all.variables', colnames(all.variables), sep = '$')
            temp.data <- t(temp.data)
            adjusted.data <- lm(as.formula(paste('temp.data', paste0(lm.formula, collapse = '+') , sep = '~')))
            temp.data <- t(adjusted.data$residuals)
            colnames(temp.data) <- colnames(se.obj)
            row.names(temp.data) <- row.names(se.obj)
        }
    }
    # select data input for clustering ####
    printColoredMessage( message = '-- Select input data for clustering:',
        color = 'magenta',
        verbose = verbose)
    if(approach == 'pca' & is.null(ncg)){
        printColoredMessage(
            message = paste0(
                '- Apply PCA on the data and use the first ',
                nb.pcs,
                ' PCs as an input for clustering.'),
            color = 'blue',
            verbose = verbose)
        sv.dec <- runSVD(
            x = t(temp.data),
            k = nb.pcs,
            BSPARAM = svd.bsparam,
            center = center,
            scale = scale)
        input.data = sv.dec$u
        if(clustering.methods == 'nbClust'){
            input.data.name <- paste0(approach, 'onAllGenes_', nbClust.method, 'Clustering')
        } else input.data.name <- paste0(approach, 'onAllGenes_', clustering.methods, 'Clustering')
    } else if (approach == 'pca' & !is.null(ncg)){
        printColoredMessage(
            message = paste0(
                '- Apply PCA on the data using the "ncg" gene only, and use the first ',
                nb.pcs,
                ' PCs as an input for clustering.'),
            color = 'blue',
            verbose = verbose)
        sv.dec <- runSVD(
            x = t(temp.data[ncg , ]),
            k = nb.pcs,
            BSPARAM = svd.bsparam,
            center = center,
            scale = scale)
        input.data = sv.dec$u
        if(clustering.methods == 'nbClust'){
            input.data.name <- paste0(approach, 'onNCG_', nbClust.method, 'Clustering')
        } else input.data.name <- paste0(approach, 'onNCG_', clustering.methods, 'Clustering')
    } else if (approach == 'rle'){
        if(is.null(ncg)){
            printColoredMessage(
                message = paste0('-Apply RLE on the data.'),
                color = 'blue',
                verbose = verbose)
            rle.data <- temp.data - rowMedians(temp.data)
            if(clustering.methods == 'nbClust'){
                input.data.name <- paste0(approach, rle.comp, 'onAllGenes_', nbClust.method, 'Clustering')
            } else input.data.name <- paste0(approach, rle.comp,'onAllGenes_', clustering.methods, 'Clustering')
        } else if (!is.null(ncg)){
            printColoredMessage(
                message = paste0('-Apply RLE on the data using only "ncg" genes.'),
                color = 'blue',
                verbose = verbose)
            rle.data <- temp.data[ncg , ] - rowMedians(temp.data[ncg , ])
            if(clustering.methods == 'nbClust'){
                input.data.name <- paste0(approach, 'onNCG_', nbClust.method, 'Clustering')
            } else input.data.name <- paste0(approach, 'onNCG_', clustering.methods, 'Clustering')
        }
        if(rle.comp == 'median'){
            input.data <- colMedians(rle.data)
        } else if (rle.comp == 'iqr'){
            input.data <- colIQRs(rle.data)
        } else if (rle.comp == 'both'){
            input.data <- list(
                rle.med = colMedians(rle.data),
                rle.iqr = colIQRs(rle.data))
        }
    } else if(approach == 'sample.scoring'){
        if(is.null(uv.gene.sets)){
            all.uv.gene.sets <- list(ncg = ncg)
        } else if (!is.null(uv.gene.sets) & !is.null(ncg)){
            all.uv.gene.sets <- uv.gene.sets[['ncg']] <- ncg
        } else if (!is.null(uv.gene.sets) & is.null(ncg)){
            all.uv.gene.sets <- uv.gene.sets
        }
        printColoredMessage(
            message = 'Calculate sample scores for individual gene sets of the uv.gene.sets and use the scores as input for clustering.',
            color = 'blue',
            verbose = verbose)
        ranked.data <- rankGenes(temp.data)
        input.data <- lapply(
            names(all.uv.gene.sets),
            function(x) singscore::simpleScore(rankData = ranked.data, upSet = all.uv.gene.sets[[x]])$TotalScore)
        names(input.data) <- names(all.uv.gene.sets)
        rm(ranked.data)
        gc()
        if(clustering.methods == 'nbClust'){
            input.data.name <- paste0(approach, '_', nbClust.method, 'Clustering')
        } else input.data.name <- paste0(approach, '_', clustering.methods, 'Clustering')
    }
    # clustering ####
    printColoredMessage(message = '- Cluster the data',
        color = 'magenta',
        verbose = verbose)
    if(clustering.methods == 'kmeans'){
        printColoredMessage(
            message = paste0('- Apply kmeans with centers = ', nb.clusters, ' on the data'),
            color = 'blue',
            verbose = verbose)
        ## kmeans ####
        set.seed(3344)
        if(is.list(input.data)){
            uv.sources <- lapply(
                names(input.data),
                function(x){
                    groups <- kmeans(x = input.data[[x]], centers = nb.clusters, iter.max = 10000)$cluster
                    paste0(input.data.name, '_', x, '_batch' , groups)
                })
            names(uv.sources) <- names(input.data)
        } else {
            groups <- kmeans(x = input.data, centers = nb.clusters, iter.max = 10000)$cluster
            uv.sources <- paste0(input.data.name, '_batch' , groups)
            }
    } else if (clustering.methods == 'cut'){
        printColoredMessage(
            message = paste0('- Apply the cut method with breaks = ', nb.clusters, ' on the data'),
            color = 'blue',
            verbose = verbose)
        ## cut ####
        if(is.list(input.data)){
            uv.sources <- lapply(
                names(input.data),
                function(x){
                    groups <- as.numeric(cut(x = input.data[[x]], breaks = nb.clusters, include.lowest = TRUE))
                    paste0(input.data.name, '_', x, '_batch' , groups)
                })
            names(uv.sources) <- names(input.data)
        } else {
            groups <- as.numeric(cut(x = input.data, breaks = nb.clusters, include.lowest = TRUE))
            uv.sources <- paste0(input.data.name, '_batch' , groups)
            }
    } else if(clustering.methods == 'quantile'){
        printColoredMessage(
            message = paste0('- Apply the quantile method with probs = ',
                             paste0(round(seq(0, 1, 1/nb.clusters), digits = 2)), ' on the data'),
            color = 'blue',
            verbose = verbose)
        ## quantile ####
        if(is.list(input.data)){
            uv.sources <- lapply(
                names(input.data),
                function(x){
                    quantiles <- quantile(x = input.data[[x]], probs = seq(0, 1, 1 / nb.clusters))
                    groups <- as.numeric(cut(x = input.data[[x]], breaks = quantiles, include.lowest = TRUE))
                    paste0(input.data.name, '_', x, '_batch' , groups)
                })
            names(uv.sources) <- names(input.data)
        } else {
            quantiles <- quantile(x = input.data, probs = seq(0, 1, 1 / nb.clusters))
            groups <- as.numeric(cut(x = input.data, breaks = quantiles, include.lowest = TRUE))
            uv.sources <- paste0(input.data.name, '_batch' , groups)
        }
    } else if(clustering.methods == 'nbClust'){
        printColoredMessage(
            message = '- Apply the nbClust method on the data.',
            color = 'blue',
            verbose = verbose)
        ## nbClust ####
        if(!is.list(input.data)){
            initial.clusters <- NbClust(
                data = input.data,
                diss = nbClust.diss,
                distance = nbClust.distance,
                min.nc = nbClust.min.nc,
                max.nc = nbClust.max.nc,
                method = nbClust.method,
                index = nbClust.index,
                alphaBeale = nbClust.alphaBeale
            )
            batch.samples <- data.frame(
                id = colnames(se.obj),
                batch = initial.clusters$Best.partition
            )
            selected.clusters <- findRepeatingPatterns(
                vec = batch.samples$batch,
                n.repeat = round(max.samples.per.batch * ncol(se.obj), digits = 0)
            )
            while(length(selected.clusters) > 0){
                more.clusters <- lapply(
                    selected.clusters,
                    function(x){
                        index <- batch.samples$batch == x
                        if(is.matrix(input.data)){
                            sub.input.data <- input.data[index , ]
                        } else sub.input.data <- input.data[index]
                        sub.clusters <- NbClust(
                            data = sub.input.data,
                            diss = nbClust.diss,
                            distance = nbClust.distance,
                            min.nc = nbClust.min.nc,
                            max.nc = nbClust.max.nc,
                            method = nbClust.method,
                            index = nbClust.index,
                            alphaBeale = nbClust.alphaBeale)
                        data.frame(
                            id = batch.samples$id[index],
                            batch = paste0(x, sub.clusters$Best.partition))
                    })
                more.clusters <- do.call(rbind, more.clusters)
                batch.samples$batch[match(more.clusters$id, batch.samples$id)] <- more.clusters$batch
                selected.clusters <- findRepeatingPatterns(
                    vec = batch.samples$batch,
                    n.repeat = round(max.samples.per.batch * ncol(se.obj), digits = 0))
            }
            uv.sources <- paste0(input.data.name, '_nbClust_batch', batch.samples$batch)
        } else {
            uv.sources <- lapply(
                names(input.data),
                function(x){
                    initial.clusters <- NbClust(
                        data = input.data[[x]],
                        diss = nbClust.diss,
                        distance = nbClust.distance,
                        min.nc = nbClust.min.nc,
                        max.nc = nbClust.max.nc,
                        method = nbClust.method,
                        index = nbClust.index,
                        alphaBeale = nbClust.alphaBeale
                    )
                    batch.samples <- data.frame(
                        id = colnames(se.obj),
                        batch = initial.clusters$Best.partition)
                    selected.clusters <- findRepeatingPatterns(
                        vec = batch.samples$batch,
                        n.repeat = round(max.samples.per.batch * ncol(se.obj), digits = 0))
                    while(length(selected.clusters) > 0){
                        more.clusters <- lapply(
                            selected.clusters,
                            function(y){
                                index <- batch.samples$batch == y
                                if(is.matrix(input.data)){
                                    sub.input.data <- input.data[[x]][index , ]
                                } else sub.input.data[[x]][index]
                                sub.clusters <- NbClust(
                                    data = sub.input.data,
                                    diss = nbClust.diss,
                                    distance = nbClust.distance,
                                    min.nc = nbClust.min.nc,
                                    max.nc = nbClust.max.nc,
                                    method = nbClust.method,
                                    index = nbClust.index,
                                    alphaBeale = nbClust.alphaBeale
                                )
                                data.frame(
                                    id = batch.samples$id[index],
                                    batch = paste0(y, sub.clusters$Best.partition))
                            })
                        more.clusters <- do.call(rbind, more.clusters)
                        all.equal(more.clusters$id, batch.samples$id)
                        batch.samples$batch[match( more.clusters$id, batch.samples$id )] <- more.clusters$batch
                        selected.clusters <- findRepeatingPatterns(
                            vec = batch.samples$batch,
                            n.repeat = round(max.samples.per.batch * ncol(se.obj), digits = 0))
                    }
                    return(paste0(x, '_nbclust_batch', batch.samples$batch))
                })
            names(uv.sources) <- names(input.data)
        }
    }
    # number of possible batches ####
    printColoredMessage(
        message = paste0(length(unique(uv.sources)),' potential batches are found in the ', assay.name, ' data.'),
        color = 'blue',
        verbose = verbose)
    # saving the data results ####
    printColoredMessage(message = '- Save the the results:',
                        color = 'magenta',
                        verbose = verbose)
    if(save.se.obj == TRUE){
        if (!'SourcesOfUV' %in%  names(se.obj@metadata)) {
            se.obj@metadata[['SourcesOfUV']] <- list()
        }
        if (!assay.name %in%  names(se.obj@metadata[['SourcesOfUV']])){
            se.obj@metadata[['SourcesOfUV']][['assay.name']] <- list()
        }
        if (!'Unknown' %in%  names(se.obj@metadata[['SourcesOfUV']][[assay.name]])){
            se.obj@metadata[['SourcesOfUV']][['assay.name']][['Unknown']] <- list()
        }
        if (!input.data.name %in%  names(se.obj@metadata[['SourcesOfUV']][[assay.name]][['Unknown']])){
            se.obj@metadata[['SourcesOfUV']][[assay.name]][['Unknown']][[input.data.name]] <- list()
        }
        se.obj@metadata[['SourcesOfUV']][[assay.name]][['Unknown']][[input.data.name]][['batches']] <- uv.sources
        se.obj@metadata[['SourcesOfUV']][[assay.name]][['Unknown']][[input.data.name]][['input.data']] <- input.data
        printColoredMessage(
            message = 'The potentail unknow sources of variation are saved to the metadata of the SummarizedExperiment object',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The indentifyUnknownUV function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj)
    } else {
        printColoredMessage(
            message = 'The results are outputed as list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The indentifyUnknownUV function finished.',
            color = 'white',
            verbose = verbose)
        return(list(batches = uv.sources, input.data = input.data))
    }
}




