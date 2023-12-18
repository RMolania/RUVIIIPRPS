#' is used to find potential unknown sources of unwanted variation..
#'
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name A symbol indicates an assay name in SummarizedExperiment object.
#' @param approach A symbol indicates which approach should be used to identify unknown sources of unwanted variation. RLE:
#' @param regress.out.bio.variables One or more symbol indicate the names of columns in the sample annotation (colData) that will be regress out
#' from the data before finding sources of unwanted variation.
#' @param regress.out.bio.gene.sets List. A list of biological gene signatures. Individual gene set will be used to score samples and
#' then each scores will be regressed out from the data before finding negative control genes.
#' @param uv.gene.sets List. A list of unwanted variation gene signatures. Individual gene set will be used to score samples and
#' then each scores will be clustered to find possible unknown batches. The default is null.
#' @param ncg Logical or factor. Indicates a set of negative control genes to find sources of unwanted variation. If is not NULL,
#' the RLE or PCA and sample scoring will be based on only negative control genes. The default is null.
#' @param clustering.methods A symbol. Indicates winch clustering methods: 'kmeans', 'cut', 'quantile', 'nbClust' should be used to find
#' possible unknown batches. The default is nbClust.
#' @param nbClust.diss dissimilarity matrix to be used. By default, diss=NULL, but if it is replaced by a dissimilarity matrix, distance should be "NULL".
#' @param nbClust.distance the distance measure to be used to compute the dissimilarity matrix. This must be one of: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski" or "NULL".
#' By default, distance="euclidean". If the distance is "NULL", the dissimilarity matrix (diss) should be given by the user. If distance is not "NULL", the dissimilarity matrix should be "NULL".
#' @param nbClust.min.nc minimal number of clusters, between 1 and (number of objects - 1)
#' @param nbClust.max.nc minimal number of clusters, between 1 and (number of objects - 1)
#' @param nbClust.method maximal number of clusters, between 2 and (number of objects - 1), greater or equal to min.nc. By default, max.nc=15.
#' @param nbClust.index The index to be calculated. This should be one of : "kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette",
#' "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw", "all" (all indices except GAP, Gamma, Gplus and Tau),
#' "alllong" (all indices with Gap, Gamma, Gplus and Tau included).
#' @param nbClust.alphaBeale significance value for Beale's index.
#' @param max.samples.per.batch Numeric. Indicates the maximum number of samples per cluster when the clustering.methods is nbClust.
#' The default is .1 (10%) of total samples in the SummarizedExperiment object.
#' @param nb.clusters Numeric. Indicates how many clusters should be found when clustering.methods kmeans, cut and quantile.
#' @param rle.comp A symbol. Indicates which properties: 'median' or 'IQR' or 'both' of the RLE data should be used for clustering when the approach is RLE.
#' The default is median.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param nb.pcs Numeric. A value to select the first principal components for clustering when the approach is set to 'PCA'.
#' @param center Logical. Indicates whether to center the data or not before performing PCA. The default is TRUE.
#' @param scale Logical. Indicates whether to scale the data or not before performing PCA. The default is FASLE.
#' @param BSPARAM TBBBB
#' @param remove.na A symbol. Indicates whether to remove NA or missing values from either the 'assays', 'sample.annotation', 'both' or 'none'.
#' If 'assays' is selected, the genes that contains NA or missing values will be excluded. If 'sample.annotation' is selected, the
#' samples that contains NA or missing values for any 'bio.variables' and 'uv.variables' will be excluded. By default, it is set to
#' 'both'.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' By default it is set to TRUE.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the
#' execution of the functions, by default it is set to TRUE.

#' @importFrom singscore rankGenes simpleScore
#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocSingular bsparam runSVD
#' @importFrom NbClust NbClust
#' @export

## correlate with known others uv
indentifyUnknownUV <- function(
        se.obj,
        assay.name,
        approach = 'RLE',
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
        nb.clusters = NULL,
        rle.comp = 'median',
        apply.log = TRUE,
        pseudo.count = 1,
        nb.pcs = 2,
        center = TRUE,
        scale = FALSE,
        BSPARAM = bsparam(),
        remove.na = 'none',
        assess.se.obj = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The indentifyUnknownUV function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if (length(assay.name) > 1) {
        stop('Please provide only one assay name.')
    }
    if (!is.null(regress.out.bio.variables)) {
        if (!regress.out.bio.variables %in% colnames(colData(se.obj))) {
            stop('The regress out bio variables are not found in the SummarizedExperiment object.')
        }
    }
    if(!is.null(regress.out.bio.gene.sets)){
        lapply(
            regress.out.bio.gene.sets,
            function(x){
                if(sum(regress.out.bio.gene.sets %in% row.names(se.obj)) == 0){
                    stop('The regress.out.bio.gene.sets are not found in the SummarizedExperiment object.')
                }
            })
    }
    if(!is.null(uv.gene.sets)){
        lapply(
            names(uv.gene.sets),
            function(x){
                if(sum(uv.gene.sets[[x]] %in% row.names(se.obj)) < 2){
                    stop('The uv.gene.sets are not found in the SummarizedExperiment object.')
                }
            })
    }
    if (!approach %in% c('RLE', 'PCA', 'sample.scoring') ) {
        stop('The approach should be one of the RLE,PCA or sample.scoring.')
    }
    if (!is.null(ncg)) {
        if (sum(ncg %in% row.names(se.obj)) == 0) {
            stop('The ncg genes are found in the SummarizedExperiment object.')
        }
    }
    if (length(clustering.methods) > 1) {
        stop('A single method should be provided for the clustering.methods.')
    }
    if (!clustering.methods %in% c('kmeans', 'cut', 'quantile', 'nbClust')) {
        stop('The clustering.methods should be one of kmeans, cut, quantile or nbClust.')
    }
    if (pseudo.count < 0) {
        stop('The valuse of pseudo.count cannot be negative.')
    }
    if(approach == 'PCA' & nb.pcs > 1 & clustering.methods %in% c('cut', 'quantile')){
        stop(paste0('The nb.pcs should be 1 to use the ', clustering.methods, ' method for clustering.'))
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
    # data transformation_ log 2####
    printColoredMessage(
        message = '### Data transformation:',
        color = 'magenta',
        verbose = verbose)
    if (apply.log) {
        printColoredMessage(
            message = paste0('Applying log2 transformation of the ', assay.name, ' + ', pseudo.count, '.'),
            color = 'blue',
            verbose = verbose)
        temp.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else {
        printColoredMessage(
            message = paste0('It seems the ', assay.name, ' data is already log transformed.'),
            color = 'blue',
            verbose = verbose)
        temp.data <- assay(x = se.obj, i = assay.name)
    }
    # sample scoring for regress.out.bio.gene.sets ####
    if(!is.null(regress.out.bio.gene.sets)){
        printColoredMessage(
            message = 'Calculate sample scores for individual gene sets of the regress.out.bio.gene.sets',
            color = 'blue',
            verbose = verbose)
        ranked.data <- rankGenes(temp.data)
        regress.out.bio.gene.sets <- sapply(
            regress.out.bio.gene.sets,
            function(x) singscore::simpleScore(rankData = ranked.data, upSet = x)$TotalScore)
        rm(ranked.data)
        gc()
    }
    # regress out bio variables and bio gene sets ####
    printColoredMessage(
        message = '### Regress out variables and gene sets from the data:',
        color = 'magenta',
        verbose = verbose)
    printColoredMessage(
        message = paste0(
            'We do not recommend regressing out any biological variation',
            'if they may be largely associated with the unwanted variation.'),
        color = 'red',
        verbose = verbose)
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
    # select data input for clustering ####
    printColoredMessage(
        message = '### Selecting input data for clustering:',
        color = 'magenta',
        verbose = verbose)
    if(approach == 'PCA' & is.null(ncg)){
        printColoredMessage(
            message = paste0(
                'Applying PCA on the transformed data and use the first ',
                nb.pcs,
                ' PCs as an input for clustering.'),
            color = 'blue',
            verbose = verbose)
        sv.dec <- runSVD(
            x = t(temp.data),
            k = nb.pcs,
            BSPARAM = BSPARAM,
            center = center,
            scale = scale)
        input.data = sv.dec$u
        if(clustering.methods == 'nbClust'){
            input.data.name <- paste0(approach, 'onAllGenes_', nbClust.method, 'Clustering')
        } else input.data.name <- paste0(approach, 'onAllGenes_', clustering.methods, 'Clustering')
    } else if (approach == 'PCA' & !is.null(ncg)){
        printColoredMessage(
            message = paste0(
                'Applying PCA on the transformed data using the ncg only, and use the first ',
                nb.pcs,
                ' PCs as an input for clustering.'),
            color = 'blue',
            verbose = verbose)
        sv.dec <- runSVD(
            x = t(temp.data[ncg , ]),
            k = nb.pcs,
            BSPARAM = BSPARAM,
            center = center,
            scale = scale)
        input.data = sv.dec$u
        if(clustering.methods == 'nbClust'){
            input.data.name <- paste0(approach, 'onNCG_', nbClust.method, 'Clustering')
        } else input.data.name <- paste0(approach, 'onNCG_', clustering.methods, 'Clustering')
    } else if (approach == 'RLE'){
        if(is.null(ncg)){
            printColoredMessage(
                message = paste0('Applying RLE on the transformed data.'),
                color = 'blue',
                verbose = verbose)
            rle.data <- temp.data - rowMedians(temp.data)
            if(clustering.methods == 'nbClust'){
                input.data.name <- paste0(approach, 'onAllGenes_', nbClust.method, 'Clustering')
            } else input.data.name <- paste0(approach, 'onAllGenes_', clustering.methods, 'Clustering')
        } else if (!is.null(ncg)){
            printColoredMessage(
                message = paste0('Applying RLE on the transformed data using only ncg.'),
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
    printColoredMessage(
        message = '### Clustering the data',
        color = 'magenta',
        verbose = verbose)
    if(clustering.methods == 'kmeans'){
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
                        } else sub.input.data[index]
                        print(x)
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
    if(save.se.obj == TRUE){
        se.obj@metadata[['UV']][['Unknown']][[input.data.name]][['batches']] <- uv.sources
        se.obj@metadata[['UV']][['Unknown']][[input.data.name]][['input.data']] <- input.data
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
            message = '------------The indentifyUnknownUV function finished.',
            color = 'white',
            verbose = verbose)
        return(list(batches = uv.sources, input.data = input.data))
    }
}




