#' Identify potential unknown sources of unwanted variation in RNA-seq data.

#' @author Ramyar Molania

#' @description
#' Identifies sources of unwanted variation in RNA-seq data when none are known.

#' @details
#' This function uses three different approaches: 'rle', 'pca' and 'sample.scoring' to find potential sources of unwanted
#' variation in RNA-seq data when none are known. In the 'rle' approach, a clustering method specified by
#' 'clustering.methods' will be applied on either the RLE medians or IQRs or both separately. In the absence of unwanted
#' variation,there should be no clearly distinguishable clusters. In the 'pca' approach, first, a principal component analysis on
#' either a set of negative control genes or all genes will be applied and then a clustering method will be used to cluster the
#' first principal components to find unknown sources of unwanted variation. In the 'sample.scoring' approach,
#' first, individual samples will be scored against a set of gene set(s) e.g. housekeeping gene, whose visible variation
#' indicates the existence of unwanted variation s, then clustering method will be applied on the scoring to find potential
#' unknown sources of unwanted variation.

#' @references
#' Gandolfo L. C. & Speed, T. P., RLE plots: visualizing unwanted variation in high dimensional data. PLoS ONE, 2018.\
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. A symbol for the selection of the name of the assay in the SummarizedExperiment object to be
#' used to identify possible sources of unwanted variation.
#' @param approach Symbol. A symbol indicating the approach to be employed for identifying unknown sources of unwanted
#' variation. This should be one of 'rle', 'pca', or 'sample.scoring'.In the rle approach, a clustering method is applied
#' to either medians or IQR or both of the rle data. In the absence of unwanted variation, no distinguishable clustering
#' should occur. In the 'pca' approach, a principal component analysis is initially performed on either a set of negative
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
#' The data must be in log transformation before computing RLE or PCA.
#' @param pseudo.count Numeric. A value serving as a pseudo count to be added to all measurements in the assay(s) before
#' applying log-transformation. This helps prevent -Inf values for measurements equal to 0. The default is 1.
#' @param nb.pcs Numeric. A value determining the number of first principal components to be selected for clustering
#' when the approach is set to 'PCA'. The default is 2.
#' @param center Logical. Indicates whether to center the data or not before performing PCA. The default is 'TRUE'.
#' @param scale Logical. Indicates whether to scale the data or not before performing PCA. The default is 'FASLE'.
#' @param svd.bsparam A BiocParallelParam object specifying how parallelization should be performed. We refer to the
#' BiocParallelParam R package for more details. The default is bsparam().
#' We refer to the 'runSVD' function from the BiocSingular R package for more details.
#' @param remove.current.estimates Symbol. Specifies whether to remove the current estimates of the unknown batch in the
#' SummarizedExperiment object or not. The default is set to 'TRUE'.
#' @param output.name Symbol. A symbol specifies the name of the output file. If is 'NULL', the functions creates a name.
#' @param assess.se.obj Logical. Whether to assess the SummarizedExperiment object or not. If 'TRUE', the function
#' 'checkSeobj' will be applied. The default is 'TRUE'.
#' @param remove.na A symbol. Indicates whether to remove NA or missing values from either the 'assays', sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be
#' excluded. If 'sample.annotation' is selected, the samples that contains NA or missing values for any 'bio.variables'
#' and 'uv.variables' will be excluded. By default, it is set to both'.
#' @param save.se.obj Logical. Indicates whether to save the results in the metadata of the SummarizedExperiment class.
#' Th default is 'TRUE', the results will be save to the 'metadata$UV$Unknown'.
#' @param plot.output Logical. When set to 'TRUE', the function generates a plot of the input data for clustering, with
#' colors representing the identified groups. The default value is 'TRUE
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom singscore rankGenes simpleScore
#' @importFrom BiocSingular bsparam runSVD
#' @importFrom NbClust NbClust
#' @importFrom stats as.formula
#' @importFrom GGally ggpairs
#' @importFrom RColorBrewer brewer.pal.info
#' @import ggplot2
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
        remove.current.estimates = FALSE,
        output.name = NULL,
        assess.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        plot.output = TRUE,
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
    if(is.logical(output.name)){
        stop('The "output.name" should be eitehr NULL or a character.')
    }

    # remove current estimates for the assay ####
    printColoredMessage(
        message = paste0('The current estimated unknown batches:'),
        color = 'magenta',
        verbose = verbose)
    if(isTRUE(remove.current.estimates)){
        if (!'UknownUV' %in%  names(se.obj@metadata)) {
            printColoredMessage(
                message = paste0('There is not any estimated unknown batches for the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
        } else  if (assay.name %in% names(se.obj@metadata[['UknownUV']])) {
            printColoredMessage(
                message = paste0('The current estimated unknown batches for the  ', assay.name, ' data is removed.'),
                color = 'blue',
                verbose = verbose)
            se.obj@metadata[['UknownUV']][[assay.name]] <- list()
        } else{
            printColoredMessage(
                message = paste0('There is not any estimated unknown batches for the  ', assay.name, ' data.'),
                color = 'blue',
                verbose = verbose)
        }
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
        colnames(input.data) <- c(paste0('PC', 1:ncol(input.data)))
        if(clustering.methods == 'nbClust'){
            input.data.name <- paste0(approach, '|AllGenes_nbClust.', nbClust.method, 'Clustering')
        } else input.data.name <- paste0(approach, '|AllGenes|', clustering.methods, 'Clustering')
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
            input.data.name <- paste0(approach, '|NCG|nbClust.', nbClust.method, 'Clustering')
        } else input.data.name <- paste0(approach, '|NCG|', clustering.methods, 'Clustering')
    } else if (approach == 'rle'){
        if(is.null(ncg)){
            printColoredMessage(
                message = paste0('-Apply RLE on the data.'),
                color = 'blue',
                verbose = verbose)
            rle.data <- temp.data - rowMedians(temp.data)
            if(clustering.methods == 'nbClust'){
                input.data.name <- paste0(approach, '.',rle.comp, '|AllGenes|nbClust.', nbClust.method, 'Clustering')
            } else input.data.name <- paste0(approach, '.', rle.comp,'|AllGenes|', clustering.methods, 'Clustering')
        } else if (!is.null(ncg)){
            printColoredMessage(
                message = paste0('-Apply RLE on the data using only "ncg" genes.'),
                color = 'blue',
                verbose = verbose)
            rle.data <- temp.data[ncg , ] - rowMedians(temp.data[ncg , ])
            if(clustering.methods == 'nbClust'){
                input.data.name <- paste0(approach, '.', rle.comp, '|NCG|nbClust.', nbClust.method, 'Clustering')
            } else input.data.name <- paste0(approach, '.', rle.comp,'|NCG|', clustering.methods, 'Clustering')
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
            input.data.name <- paste0(approach, '|nbClust.', nbClust.method, 'Clustering')
        } else input.data.name <- paste0(approach, '|', clustering.methods, 'Clustering')
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
                    paste0('Batch' , groups)
                })
            names(uv.sources) <- names(input.data)
        } else {
            groups <- kmeans(x = input.data, centers = nb.clusters, iter.max = 10000)$cluster
            uv.sources <- paste0('Batch' , groups)
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
                    paste0('Batch' , groups)
                })
            names(uv.sources) <- names(input.data)
        } else {
            groups <- as.numeric(cut(x = input.data, breaks = nb.clusters, include.lowest = TRUE))
            uv.sources <- paste0('Batch' , groups)
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
                    paste0('Batch' , groups)
                })
            names(uv.sources) <- names(input.data)
        } else {
            quantiles <- quantile(x = input.data, probs = seq(0, 1, 1 / nb.clusters))
            groups <- as.numeric(cut(x = input.data, breaks = quantiles, include.lowest = TRUE))
            uv.sources <- paste0('Batch' , groups)
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
            uv.sources <- paste0('Batch', batch.samples$batch)
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
                    return(paste0('Batch', batch.samples$batch))
                })
            names(uv.sources) <- names(input.data)
        }
    }
    # number of possible batches ####
    printColoredMessage(
        message = paste0(
            length(unique(uv.sources)),
            ' potential batches are found in the ',
            assay.name,
            ' data.'
        ),
        color = 'blue',
        verbose = verbose
    )

    # plot outputs ####
    currentCols <-  c(
        RColorBrewer::brewer.pal(8, "Dark2")[-5],
        RColorBrewer::brewer.pal(10, "Paired"),
        RColorBrewer::brewer.pal(12, "Set3"),
        RColorBrewer::brewer.pal(9, "Blues")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Oranges")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Greens")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Purples")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Reds")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Greys")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "BuGn")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "PuRd")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "BuPu")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "YlGn")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(10, "Paired")
    )
    colors.selected <- currentCols[1:length(unique(uv.sources))]
    if (!is.matrix(input.data)) {
        data.to.plot <- data.frame(input.data = input.data,
                                   # samples = c(1:ncol(se.obj)),
                                   batches = factor(
                                       x = paste0('Batch', as.numeric(as.factor(uv.sources))),
                                       levels = paste0('Batch', sort(unique(
                                           as.numeric(as.factor(uv.sources))
                                       )))
                                   ))
        data.to.plot <-
            data.to.plot[order(data.to.plot$batches) ,]
        data.to.plot$samples <- c(1:ncol(se.obj))
        p <-
            ggplot(data = data.to.plot, aes(x = samples, y = input.data, color = batches)) +
            geom_point() +
            ggtitle('Possible sources of batches') +
            scale_color_manual(values = colors.selected, name = 'Batch') +
            theme(
                panel.background = element_blank(),
                legend.key = element_blank(),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 12),
                axis.title.y = element_text(size = 12),
                axis.text.x = element_text(size = 9),
                axis.text.y = element_text(size = 9)
            ) +
            guides(colour = guide_legend(override.aes = list(size = 5)))
        if (isTRUE(plot.output))
            print(p)

    } else{
        if (ncol(input.data) == 1) {
            data.to.plot <- as.data.frame(input.data)
            data.to.plot$batches <- factor(x = paste0('Batch', as.numeric(as.factor(uv.sources))),
                                           levels = paste0('Batch', sort(unique(
                                               as.numeric(as.factor(uv.sources))
                                           ))))
            data.to.plot$samples <- c(1:ncol(se.obj))
            data.to.plot <-
                data.to.plot[order(data.to.plot$batches) ,]
            p <-
                ggplot(data = data.to.plot, aes(
                    x = samples,
                    y = input.data,
                    color = batches
                )) +
                geom_point() +
                ggtitle('Possible sources of batches') +
                scale_color_manual(values = colors.selected, name = 'Batch') +
                theme(
                    panel.background = element_blank(),
                    legend.key = element_blank(),
                    legend.text = element_text(size = 12),
                    legend.title = element_text(size = 14),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 12),
                    axis.title.y = element_text(size = 12),
                    axis.text.x = element_text(size = 9),
                    axis.text.y = element_text(size = 9)
                ) +
                guides(colour = guide_legend(override.aes = list(size = 5)))
            if (isTRUE(plot.output))
                print(p)

        } else{
            data.to.plot <- as.data.frame(input.data)
            data.to.plot$batches <- factor(x = paste0('Batch', as.numeric(as.factor(uv.sources))),
                                           levels = paste0('Batch', sort(unique(
                                               as.numeric(as.factor(uv.sources))
                                           ))))
            p <- GGally::ggpairs(
                data = data.to.plot[, 1:(ncol(data.to.plot) - 1)],
                mapping = ggplot2::aes(colour = data.to.plot[, ncol(data.to.plot)]),
                upper = NULL
            ) +
                scale_color_manual(values = colors.selected)
            if (isTRUE(plot.output))
                print(p)
        }
    }

    # out put name ####
    if(is.null(output.name)){
        input.data.name <- paste0(length(unique(uv.sources)), 'batches|', input.data.name)
    } else input.data.name <- output.name


    # saving the data results ####
    printColoredMessage(message = '- Save the the results:',
                        color = 'magenta',
                        verbose = verbose)
    if(save.se.obj == TRUE){
        if (!'UknownUV' %in%  names(se.obj@metadata)) {
            se.obj@metadata[['UknownUV']] <- list()
        }
        if (!assay.name %in%  names(se.obj@metadata[['UknownUV']])){
            se.obj@metadata[['UknownUV']][[assay.name]] <- list()
        }
        if (!input.data.name %in%  names(se.obj@metadata[['UknownUV']][[assay.name]])){
            se.obj@metadata[['UknownUV']][[assay.name]][[input.data.name]] <- list()
        }
        se.obj@metadata[['UknownUV']][[assay.name]][[input.data.name]][['batches']] <- uv.sources
        se.obj@metadata[['UknownUV']][[assay.name]][[input.data.name]][['input.data']] <- input.data
        se.obj@metadata[['UknownUV']][[assay.name]][[input.data.name]][['plot']] <- p
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
        return(list(
            batches = paste0('Batch', as.numeric(as.factor(uv.sources))),
            input.data = input.data,
            plot = p ))
    }
}




