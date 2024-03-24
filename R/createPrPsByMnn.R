#' Create PRPS sets using k and mutual nearest neighbors in RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function uses the k and mutual nearest neighbors approaches to create PRPS in the RNA-seq data. This function
#' can be used in situation that the biological variation are not known. The function applies the 'findKnn' function
#' to find similar samples per batch and then average them to create pseudo-samples. The, function uses the 'findMnn' to
#' match up pseudo samples across batches to create pseudo-replicates.

#' @param se.obj A summarized experiment object.
#' @param assay.name Symbol. A symbol for the selection of the name of the assay in the SummarizedExperiment object to be
#' used to find k nearest neighbors data.
#' @param uv.variable Symbol. Indicates the name of a column in the sample annotation of the SummarizedExperiment object.
#' The 'uv.variable' can be either categorical and continuous. If 'uv.variable' is a continuous variable, this will be
#' divided into 'nb.clusters' groups using the 'clustering.method'.
#' @param data.input Symbol. Indicates which data should be used an input for finding the k nearest neighbors data. Options
#' include: 'expr' and 'pcs'. If 'pcs' is selected, the first PCs of the data will be used as input. If 'expr' is selected,
#' the data will be selected as input.
#' @param nb.pcs Numeric. Indicates the number PCs should be used as data input for finding the k nearest neighbors. The
#' 'nb.pcs' must be set when the "data.input = PCs". The default is 2.
#' @param center.pca Logical. Indicates whether to scale the data or not. If center is TRUE, then centering is done by
#' subtracting the column means of the assay from their corresponding columns. The default is TRUE.
#' @param scale.pca Logical. Indicates whether to scale the data or not before applying SVD.  If scale is TRUE, then scaling
#' is done by dividing the (centered) columns of the assays by their standard deviations if center is TRUE, and the root
#' mean square otherwise. The default is FALSE
#' @param svd.bsparam A BiocParallelParam object specifying how parallelization should be performed. The default is bsparam().
#' We refer to the 'runSVD' function from the BiocSingular R package.
#' @param clustering.method Symbol.A symbol indicating the choice of clustering method for grouping the 'uv.variable'
#' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
#' continuous variable. The default is 3.
#' @param k.nn Numeric.The maximum number of nearest neighbors to compute. The default is set 3.
#' @param hvg Vector. A vector of the names of the highly variable genes. These genes will be used to find the anchors
#' samples across the batches. The default is NULL.
#' @param normalization Symbol. Indicates which normalization methods should be applied before finding the knn. The default
#' is 'cpm'. If is set to NULL, no normalization will be applied.
#' @param regress.out.bio.variables Symbols. Indicates the columns names that contain biological variables in the
#' SummarizedExperiment object. These variables will be regressed out from the data before finding genes that are highly
#' affected by unwanted variation variable. The default is NULL, indicates the regression will not be applied.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is TRUE.
#' Please, note, any RNA-seq data (assays) must be in log scale before computing RLE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
#' log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.
#' @param mnn.bpparam Symbol. A BiocParallelParam object specifying how parallelization should be performed. The default
#' is SerialParam(). We refer to the 'findMutualNN' function from the BiocNeighbors R package for more details.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object
#'  or to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom Seurat VariableFeatures FindIntegrationAnchors
#' @importFrom stats setNames
#' @importFrom purrr map_df
#' @export

createPrPsByMnn <- function(
        se.obj,
        assay.name,
        uv.variable,
        data.input = 'expr',
        nb.pcs = 2,
        center.pca = TRUE,
        scale.pca = FALSE,
        svd.bsparam = bsparam(),
        clustering.method = 'kmeans',
        nb.clusters = 3,
        k.nn = 2,
        hvg = NULL,
        normalization = 'CPM',
        regress.out.bio.variables = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        mnn.bpparam = SerialParam(),
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        verbose = TRUE){
    printColoredMessage(message = '------------The unsupervisedPRPSmnn function starts:',
                        color = 'white',
                        verbose = verbose)
    # assess.se.obj #####
    if(assess.se.obj){
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = uv.variable,
            remove.na = remove.na,
            verbose = verbose)
    }
    # finding knn ####
    printColoredMessage(message = paste(
            'Applying the find_knn function.',
            'For individual samples per each group variable, ',
            k.nn,
            ' knn will be found.'),
            color = 'magenta',
            verbose = verbose)

    all.knn <- findKnn(
        se.obj = se.obj,
        assay.name = assay.name,
        uv.variable = uv.variable,
        data.input = data.input,
        nb.pcs = nb.pcs,
        center.pca = center.pca,
        scale.pca = scale.pca,
        svd.bsparam = svd.bsparam,
        clustering.method = clustering.method,
        nb.clusters = nb.clusters,
        k.nn = k.nn,
        hvg = hvg,
        normalization = normalization,
        regress.out.bio.variables = regress.out.bio.variables,
        apply.log = apply.log,
        pseudo.count = pseudo.count,
        assess.se.obj = assess.se.obj,
        remove.na = remove.na,
        save.se.obj = save.se.obj,
        verbose = verbose)
    if(isTRUE(save.se.obj))
        all.knn <- all.knn@metadata$PRPS$unsupervised$KnnMnn$knn[[1]]

    # finding mnn ####
    all.mnn <- findMnn(
        se.obj = se.obj,
        assay.name = assay.name,
        uv.variable = uv.variable,
        clustering.method = clustering.method,
        nb.clusters = nb.clusters,
        hvg = hvg,
        normalization = normalization,
        regress.out.bio.variables = regress.out.bio.variables,
        apply.log = apply.log,
        pseudo.count = pseudo.count,
        assess.se.obj = assess.se.obj,
        mnn.bpparam = SerialParam(),
        remove.na = remove.na,
        save.se.obj = save.se.obj,
        verbose = verbose)
    if(isTRUE(save.se.obj))
        all.mnn <- all.mnn@metadata$PRPS$unsupervised$KnnMnn$mnn[[1]]

    # find PRPS sets ####
    all.prps.sets <- lapply(
        1:nrow(all.mnn),
        function(x){
            # ps set 1
            ps.set.1 <- all.knn[ , c(1:c(k.nn+1)) ] == all.mnn$sample.no.1[x]
            ps.set.1 <- all.knn[rowSums(ps.set.1) > 0 , ]
            ps.set.1$mnn.sets <- paste0(
                sort(c(all.mnn[x , 3], all.mnn[x , 4])),
                collapse = '_')
            ps.set.1$mnn.sets.data <- paste0(
                sort(c(all.mnn[x , 1], all.mnn[x , 2])),
                collapse = '_')
            if(nrow(ps.set.1) > 1){
                ps.set.1 <- ps.set.1[ps.set.1$rank.aver.dist == min(ps.set.1$rank.aver.dist) , ]
                ps.set.1
            }
            # ps set 2
            ps.set.2 <- all.knn[ , c(1:c(k.nn+1)) ] == all.mnn$sample.no.2[x]
            ps.set.2 <- all.knn[rowSums(ps.set.2) > 0 , ]
            ps.set.2$mnn.sets <- paste0(
                sort(c(all.mnn[x , 3], all.mnn[x , 4])),
                collapse = '_')
            ps.set.2$mnn.sets.data <- paste0(
                sort(c(all.mnn[x , 1], all.mnn[x , 2])),
                collapse = '_')
            if(nrow(ps.set.2) > 1){
                ps.set.2 <- ps.set.2[ps.set.2$rank.aver.dist == min(ps.set.2$rank.aver.dist) , ]
                ps.set.2
            }
            prps.set <- rbind(ps.set.1, ps.set.2)
            prps.set
        })
    all.prps.sets <- do.call(rbind, all.prps.sets)
    all.prps.sets$aver.mnn.sets <- unlist(lapply(
        seq(1, nrow(all.prps.sets), 2),
        function(x) rep(mean(all.prps.sets$aver.dist[x:(x+1)]), 2)))
    all.prps.sets$rank.aver.mnn.sets <- rank(all.prps.sets$aver.mnn.sets)

    # filter PRPS sets ####
    all.prps.sets <- lapply(
        unique(all.prps.sets$mnn.sets.data),
        function(x){
            temp.prps.set <- all.prps.sets[all.prps.sets$mnn.sets.data == x , ]
            if(nrow(temp.prps.set) > 10*2){
                temp.prps.set <- temp.prps.set[order(temp.prps.set$rank.aver.mnn.sets) , ]
                temp.prps.set[1:c(10*2) , ]
            } else{
                temp.prps.set <- all.prps.sets[all.prps.sets$mnn.sets.data == x , ]
            }
        })
    all.prps.sets <- do.call(rbind, all.prps.sets)

    # check coverage ####
    groups <- sort(unique(unlist(strsplit(all.prps.sets$mnn.sets.data, '_'))))
    prps.coverage <- matrix('', nrow = nrow(all.prps.sets), ncol = length(groups))
    colnames(prps.coverage) <- groups
    prps.coverage <- t(sapply(1:nrow(all.prps.sets), function(i) {
        names(prps.coverage[i, ]) %in% unlist(strsplit(all.prps.sets$mnn.sets.data[i], '_'))
    }))
    colnames(prps.coverage) <- groups
    if(sum(colSums(prps.coverage) == 0)){
        printColoredMessage(
            message = paste(paste0( colnames(prps.coverage)[ colSums(prps.coverage) == 0], collapse = ' & '), 'are not covered by any PRPS set' ),
        color = 'blue',
        verbose = verbose)
    }
    # check connection
    printColoredMessage(message = '### Asseesing the connection between PRPS sets.',
                          color = 'blue',
                          verbose = verbose)
    prps.connection <- lapply(
        1:nrow(prps.coverage),
        function(y) {
            batch.names.a <- names(which(prps.coverage[y,]  == TRUE))
            con.prps <- lapply(
                c(1:nrow(prps.coverage))[-y],
                function(z) {
                    batch.names.b <- names(which(prps.coverage[z,] ==  TRUE))
                    inter.samples <- intersect(batch.names.a, batch.names.b)
                    if (length(inter.samples) > 0) {
                        sort(unique(c(batch.names.a, batch.names.b)), decreasing = FALSE)
                    } else {
                        sort(batch.names.a, decreasing = FALSE)
                    }
                })
            sort(unique(unlist(Filter(Negate(is.null), con.prps))), decreasing = FALSE)
        })
    all.covered.batches <- Filter(Negate(is.null), prps.connection)
    all.covered.batches <- unique(all.covered.batches)
    all.not.covered.batches <- unique(unlist(all.covered.batches))[unique(unlist(all.covered.batches)) %in% groups]
    if (length(unique(unlist(all.covered.batches))) == length(groups)) {
        printColoredMessage(message = 'The connections between PRPS sets can cover all the batches.',
                              color = 'white',
                              verbose = verbose)
    } else {
        printColoredMessage(message = 'All batches are not covered by PRPS.',
                              color = 'red',
                              verbose = verbose)
        printColoredMessage(
            message = paste0(
                'The PRPS sets donot cover ',
                paste0(all.not.covered.batches, collapse = ' & '),
                ' batches.'
            ),
            color = 'white',
            verbose = verbose
        )
    }
    # create PRPS data ####
    printColoredMessage(message = '-- Create PRPS data:',
        color = 'magenta',
        verbose = verbose)
    ## data transformation ####
    printColoredMessage(
        message = '- Apply log on the data before creating PRPS:',
        color = 'blue',
        verbose = verbose)
    ## apply log ####
    if (isTRUE(apply.log) & !is.null(pseudo.count)){
        printColoredMessage(
            message = paste0('Applying log2 + ', pseudo.count, ' (pseudo.count) on the ', assay.name, ' data.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else if (isTRUE(apply.log) & is.null(pseudo.count)){
        printColoredMessage(
            message = paste0('Applying log2 on the ', assay.name, ' data.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- log2(assay(x = se.obj, i = assay.name))
    } else if (isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0('The ', assay.name, ' data will be used without any log transformation.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- assay(x = se.obj, i = assay.name)
    }
    printColoredMessage(
        message = '- Aeverage samples to create psudo samples:',
        color = 'blue', verbose = verbose)
    prps.data <- lapply(
        unique(all.prps.sets$mnn.sets),
        function(x){
            temp.prps <- all.prps.sets[all.prps.sets$mnn.sets == x, ]
            index.a <- unlist(unname(temp.prps[1, grep('overal', colnames(temp.prps))]))
            index.b <- unlist(unname(temp.prps[2, grep('overal', colnames(temp.prps))]))
            prps.a <- rowMeans(expr.data[ , index.a])
            prps.b <- rowMeans(expr.data[ , index.b])
            prps <- cbind(prps.a, prps.b)
            colnames(prps) <- paste(uv.variable, temp.prps$mnn.sets, sep = '_')
            prps
        })
    prps.data <- do.call(cbind, prps.data)
    #sanity check ####
    if(!sum(table(colnames(prps.data)) == 2) == ncol(prps.data)/2){
        stop('There someting wrong with PRPS sets.')
    }

    # saving the outputs ####
    printColoredMessage(message = '-- Save the PRPS data',
                        color = 'magenta',
                        verbose = verbose)
    if (save.se.obj) {
        printColoredMessage(
            message = 'Save all the PRPS data into the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        output.name <- paste0(uv.variable, '||' , assay.name)
        se.obj@metadata[['PRPS']][['unsupervised']][['KnnMnn']][['prps.data']][[output.name]] <- prps.data
        printColoredMessage(message = '------------The unsupervisedPRPSmnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else{
        printColoredMessage(message = '------------The unsupervisedPRPSmnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(prps.data = prps.data)
    }

}

