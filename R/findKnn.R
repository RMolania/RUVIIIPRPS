#' find k-nearest neighbors in RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function finds k nearest neighbors samples in RNA-seq data. The k nearest neighbors will be used to create pseudo
#' sample within individual batches.

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
#' @param svd.bsparam Symbol. A BiocParallelParam object specifying how parallelization should be performed. The default is bsparam().
#' We refer to the 'runSVD' function from the BiocSingular R package.
#' @param clustering.method Symbol.A symbol indicating the choice of clustering method for grouping the 'uv.variable'
#' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
#' continuous variable. The default is 3.
#' @param k Numeric.The maximum number of nearest neighbors to compute. The default is set 3.
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
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object
#'  or to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.


#' @importFrom SummarizedExperiment assay colData
#' @importFrom RANN nn2
#' @importFrom stats dist
#' @importFrom utils txtProgressBar
#' @export

findKnn <- function(
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
        k = 2,
        hvg = NULL,
        normalization = 'CPM',
        regress.out.bio.variables = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The findKnn function starts:',
                        color = 'white',
                        verbose = verbose)
    # check input #####
    if (is.list(assay.name)) {
        stop('The "assay.name" cannot be a list.')
    } else if (length(assay.name) > 1) {
        stop('The "assay.name" must be the name of signle assay in the SummarizedExperiment object.')
    } else if (is.null(uv.variable)) {
        stop('The "uv.variable" variable cannot be empty.')
    } else if (length(uv.variable) > 1) {
        stop('The "uv.variable" must contain the name of signle variable in the SummarizedExperiment object.')
    } else if (!data.input %in% c('expr', 'pcs')) {
        stop('The "data.input" must be one of the "expr" or "pcs".')
    } else if (data.input == 'pcs' & is.null(nb.pcs)) {
        stop('The valuse of "nb.pcs" must be sepcified when the data.input = pcs.')
    } else if (k == 0) {
        stop('The k cannot be 0.')
    } else if (!uv.variable %in% colnames(colData(se.obj))) {
        stop('The "uv.variable" variable cannot be found in the SummarizedExperiment object.')
    }
    if(!is.null(hvg)){
        if (sum(hvg %in% row.names(se.obj)) != length(hvg))
            stop('All the hvg genes are not found in the SummarizedExperiment object.')
    }

    # assess summarized experiment ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = uv.variable,
            remove.na = remove.na,
            verbose = verbose)
    }
    all.samples.index <- c(1:ncol(se.obj))
    ini.variable <- se.obj[[uv.variable]]

    # check the uv variable ####
    if (class(se.obj[[uv.variable]]) %in% c('integer', 'numeric')) {
        printColoredMessage(
            message = paste0(
                'The "uv.variable" is a continouse variable, then this will be divided into ',
                nb.clusters,
                ' groups using ',
                clustering.method,
                ' clustering.'),
            color = 'blue',
            verbose = verbose
        )
        if (clustering.method == 'kmeans') {
            set.seed(3456)
            uv.cont.clusters <- kmeans(
                x = colData(se.obj)[[uv.variable]],
                centers = nb.clusters,
                iter.max = 1000
            )
            se.obj[[uv.variable]] <- factor(x = paste0(uv.variable, '_uv.variable', uv.cont.clusters$cluster))
        } else if (clustering.method == 'cut') {
            uv.cont.clusters <- as.numeric(cut(
                x = colData(se.obj)[[uv.variable]],
                breaks = nb.clusters,
                include.lowest = TRUE))
            se.obj[[uv.variable]] <- factor(x = paste0(uv.variable, '_group', uv.cont.clusters))
        } else if (clustering.method == 'quantile') {
            quantiles <- quantile(
                x = colData(se.obj)[[uv.variable]],
                probs = seq(0, 1, 1 / nb.clusters))
            uv.cont.clusters <- as.numeric(cut(
                x = colData(se.obj)[[uv.variable]],
                breaks = quantiles,
                include.lowest = TRUE))
            se.obj[[uv.variable]] <- factor(x = paste0(uv.variable, '_uv.variable', uv.cont.clusters))
        }
    } else if(is.factor(se.obj[[uv.variable]])){
        se.obj[[uv.variable]] <- factor(x = se.obj[[uv.variable]])
    }
    if (length(findRepeatingPatterns(vec = se.obj[[uv.variable]], n.repeat = k)) != length(unique(se.obj[[uv.variable]])))
        stop(paste0('Some sub-groups of the variable ', uv.variable, ' have less than ', k, ' samples. '))


    # data normalization and transformation and regression ####
    printColoredMessage(message = '-- Data normalization, transformation and regression:',
                        color = 'magenta',
                        verbose = verbose)
    groups <- levels(colData(se.obj)[[uv.variable]])
    all.norm.data <- lapply(
        groups,
        function(x) {
            selected.samples <- colData(se.obj)[[uv.variable]] == x
            if (!is.null(normalization) &is.null(regress.out.bio.variables)) {
                printColoredMessage(
                    message = paste0('Applying the ', normalization,' within samples from ', x, 'group.'),
                    color = 'blue',
                    verbose = verbose
                )
                norm.data <- applyOtherNormalizations(
                        se.obj = se.obj[, selected.samples],
                        assay.name = assay.name,
                        method = normalization,
                        apply.log = apply.log,
                        pseudo.count = pseudo.count,
                        assess.se.obj = FALSE,
                        save.se.obj = FALSE,
                        remove.na = 'none',
                        verbose = FALSE)
                norm.data
            } else if (!is.null(normalization) & !is.null(regress.out.bio.variables)) {
                printColoredMessage(
                    message = paste0('Applying the ', normalization,' within samples from ',
                        x, ' group, and then regressing out ', paste0(regress.out.bio.variables, collapse = '&'), '.'),
                    color = 'blue',
                    verbose = verbose)
                # normalization
                norm.data <- applyOtherNormalizations(
                        se.obj = se.obj[, selected.samples],
                        assay.name = assay.name,
                        method = normalization,
                        apply.log = apply.log,
                        pseudo.count = pseudo.count,
                        assess.se.obj = FALSE,
                        save.se.obj = FALSE,
                        remove.na = 'none',
                        verbose = FALSE)
                sample.info <- as.data.frame(colData(se.obj[ , selected.samples]))
                # regression
                norm.data <- t(norm.data)
                lm.formual <- paste('sample.info', regress.out.bio.variables, sep = '$')
                y <- lm(as.formula(paste(
                    'norm.data',
                    paste0(lm.formual, collapse = '+') ,
                    sep = '~'
                )))
                norm.data <- t(y$residuals)
                colnames(norm.data) <- colnames(norm.data)
                row.names(norm.data) <- row.names(norm.data)
                norm.data
            } else if (is.null(normalization) & is.null(regress.out.bio.variables & apply.log & !is.null(pseudo.count))) {
                printColoredMessage(
                    message = paste0('Applying log2 +  ', pseudo.count , '(pseudo.count) within samples from ', x, ' group.'),
                    color = 'blue',
                    verbose = verbose)
                norm.data <- log2(assay(x = se.obj[ , selected.samples], i = assay.name) + pseudo.count)
                norm.data
            } else if (is.null(normalization) & is.null(regress.out.bio.variables & apply.log & is.null(pseudo.count))) {
                printColoredMessage(
                    message = paste0('Applying log2 within samples from ', x, ' group.'),
                    color = 'blue',
                    verbose = verbose)
                norm.data <- log2(assay(x = se.obj[ , selected.samples], i = assay.name))
                norm.data
            } else if (is.null(normalization) & is.null(regress.out.bio.variables & !apply.log)) {
                printColoredMessage(
                    message = paste0('No any library size normalization and transformation ','within samples from ',
                        x, 'group.'),
                    color = 'blue',
                    verbose = verbose
                )
                norm.data <- assay(x = se.obj[ , selected.samples], i = assay.name)
                norm.data
            }
        })
    names(all.norm.data) <- groups
    printColoredMessage(message = '-- Finding k nearest neighbors.',
                        color = 'magenta',
                        verbose = verbose)
    printColoredMessage(
        message = paste0('For individual samples within each subgroup variable "', uv.variable,
            '", k = ', k, ' nearest neighbours will be found. This may take few minutes.'),
        color = 'blue',
        verbose = verbose
    )
    total <- length(unique(colData(se.obj)[[uv.variable]]))
    pb <- utils::txtProgressBar(
        min = 0,
        max = total,
        style = 3)
    all.knn <- lapply(
        1:total,
        function(x) {
            norm.data <- all.norm.data[[groups[x]]]
            # PCA #####
            if (data.input == 'expr' & !is.null(hvg)) {
                message(' ')
                printColoredMessage(
                    message = paste0('Finding knn in the data from ', groups[x],' using the highly variable genes.'),
                    color = 'blue',
                    verbose = verbose
                )
                norm.data <- t(norm.data[hvg, ])
            } else if (data.input == 'expr' & is.null(hvg)) {
                message(' ')
                printColoredMessage(
                    message = paste0('Finding knn in the data from ', groups[x],' using all genes.'),
                    color = 'blue',
                    verbose = verbose
                )
                norm.data <- t(norm.data)
            } else if (data.input == 'pcs' & !is.null(hvg)) {
                message(' ')
                printColoredMessage(
                    message = paste0('Finding knn using the first ', nb.pcs, ' PCs ', 'of the data from ',
                        groups[x], ' using the highly variable genes.'),
                    color = 'blue',
                    verbose = verbose
                )
                sv.dec <- BiocSingular::runSVD(
                    x = t(norm.data[hvg, ]),
                    k = nb.pcs,
                    BSPARAM = svd.bsparam,
                    center = center.pca,
                    scale = scale.pca
                )
                norm.data <- sv.dec$u
            } else if (data.input == 'pcs' & is.null(hvg)) {
                message(' ')
                printColoredMessage(
                    message = paste0('Finding knn using the first ', nb.pcs, ' PCs', ' of the data from ',
                        groups[x], ' using all genes.'),
                    color = 'blue',
                    verbose = verbose
                )
                sv.dec <- BiocSingular::runSVD(
                    x = t(norm.data),
                    k = nb.pcs,
                    BSPARAM = svd.bsparam,
                    center = center.pca,
                    scale = scale.pca
                )
                norm.data <- sv.dec$u
            }
            # find knn #####
            knn.samples <-RANN::nn2(data = norm.data, k = c(k + 1))
            # index
            knn.index <- as.data.frame(knn.samples$nn.idx)
            colnames(knn.index) <- paste0('dataset.index', 1:c(k + 1))
            selected.samples <- se.obj[[uv.variable]] == groups[x]
            ovral.cell.no <- as.data.frame(
                apply(
                    knn.index,
                    2,
                    function(z) all.samples.index[selected.samples][z]))
            colnames(ovral.cell.no) <- paste0('overal.index', 1:c(k + 1))
            knn.index <- as.data.frame(cbind(ovral.cell.no , knn.index))
            # distance
            knn.dis <- round(as.data.frame(knn.samples$nn.dists), digits = 3)
            knn.dis <- knn.dis[, -1, drop = FALSE]
            colnames(knn.dis) <- paste0('distance1_', 2:c(k + 1))
            knn.index.dist <- cbind(knn.index, knn.dis)
            if (k > 1) {
                all.comb <- combn(x = paste0('dataset.index', 2:c(k + 1)), m = 2)
                all.comb.names <- combn(x = 2:c(k + 1), m = 2)
                for (z in 1:ncol(all.comb)) {
                    pair.dist <- unlist(lapply(
                        1:nrow(knn.index.dist),
                        function(y) {
                            col1 <- knn.index.dist[, all.comb[, z][1]][y]
                            col2 <- knn.index.dist[, all.comb[, z][2]][y]
                            stats::dist(x = norm.data[c(col1, col2) ,],
                                method = "euclidean",
                                diag = FALSE,
                                upper = FALSE)
                        }))
                    name <- paste0('dist', all.comb.names[, z][1], '_', all.comb.names[, z][2])
                    knn.index.dist[[name]] <- pair.dist
                }
            }
            knn.index.dist$aver.dist <- rowMeans(knn.index.dist[, grep('distance', colnames(knn.index.dist)), drop = FALSE])
            knn.index.dist$rank.aver.dist <- rank(knn.index.dist$aver.dist)
            knn.index.dist$group <- groups[x]
            setTxtProgressBar(pb, x)
            return(knn.index.dist)
        })
    all.knn <- do.call(rbind, all.knn)

    # saving the output ####
    message(' ')
    printColoredMessage(message = '-- Save the results.',
                        color = 'magenta',
                        verbose = verbose)
    if (save.se.obj) {
        printColoredMessage(
            message = 'Save all the knn resulst into the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        output.name <- paste0(uv.variable, '||' , assay.name)
        se.obj@metadata[['PRPS']][['unsupervised']][['KnnMnn']][['knn']][[output.name]] <- all.knn
        return(se.obj)
    } else{
        printColoredMessage(
            message = 'All the knn results are outputed as list.',
            color = 'blue',
            verbose = verbose)
        return(all.knn)
    }
}
