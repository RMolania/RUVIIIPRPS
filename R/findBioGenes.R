#' find highly variable genes in RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function uses different approaches to find possible bioloical genes.


#' @param se.obj A SummarizedExperiment object.
#' @param assay.name  TTTTT
#' @param method  TTTTT
#' @param nb.bio.genes TTTTT
#' @param bio.variables Symbol. A symbol or a vector of symbols specifying the column names of unwanted variables in
#' the sample annotation of the SummarizedExperiment object. These 'uv.variables' can be either categorical or continuous
#' variables.
#' @param uv.variables Symbol. A symbol or a vector of symbols specifying the column names of unwanted variables in
#' the sample annotation of the SummarizedExperiment object. These 'uv.variables' can be either categorical or continuous
#' variables.
#' @param nb.pcs Numeric. A value determining the number of first principal components to be selected for clustering
#' when the approach is set to 'PCA'. The default is 2.
#' @param top.pca.loadings TTTTT
#' @param svd.bsparam TTTTT
#' @param center.pca Logical. Indicates whether to center the data or not before performing PCA. The default is 'TRUE'.
#' @param scale.pca Logical. Indicates whether to scale the data or not before performing PCA. The default is 'FASLE'.
#' @param normalization TTTTT
#' @param uv.clustering.method  TTTTT
#' @param nb.uv.clusters TTTTT
#' @param bio.clustering.method  TTTTT
#' @param nb.bio.clusters TTTTT
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to 'TRUE'.
#' The data must be in log transformation before comouting RLE or PCA.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param assess.se.obj Logical. Whether to assess the SummarizedExperiment object or not. If 'TRUE', the function
#' 'checkSeobj' will be applied. The default is 'TRUE'.
#' @param remove.na Symbol. Indicates whether to remove missing values from the 'uv.variables'. The options are
#' 'sample.annotation' or 'none'. The default is 'sample.annotation', indicating the missing values from the variables
#' will be removed.
#' @param save.se.obj Logical. Indicates whether to save the results to the metadata of the SummarizedExperiment object
#' or not. If 'TRUE', all the possible homogeneous groups will be saved into "metadata$HGgroups$UVgroups", otherwise
#' the results will outputted as a vector. The default is 'TRUE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom matrixStats rowMads

findBioGenes <- function(
        se.obj,
        assay.name,
        method = 'mad',
        nb.bio.genes = 10,
        bio.variables = NULL,
        uv.variables = NULL,
        nb.pcs = 3,
        top.pca.loadings = 1,
        svd.bsparam = bsparam(),
        center.pca = TRUE,
        scale.pca = FALSE,
        normalization = 'CPM',
        uv.clustering.method = 'kmeans',
        nb.uv.clusters = 3,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 3,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The findHVG function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if(!is.vector(assay.name)){
        stop('The "assay.name" must be a single assay name.')
    } else if (length(assay.name) > 1){
        stop('The "assay.name" must be a single assay name.')
    } else if(!method %in% c('pca', 'vst', 'mad', 'twoWayAnova')){
        stop('The "method" must be of one the "pca", "vst", "mad",  and "twoWayAnova".')
    }
    if(method == 'pca'){
        if(is.null(nb.pcs)){
            stop('The "nb.pcs" must be specified when the method = pca .')
        } else if (is.null(top.pca.loadings) | top.pca.loadings == 0){
            stop('The "top.pca.loadings" must be specified or more than 0 when the method = pca.')
        } else if(top.pca.loadings > 100 | top.pca.loadings < 0){
            stop('The "top.pca.loadings" must be between 0 and 100')
        }
    }
    if(method == 'twoWayAnova'){
        if(is.null(bio.variables)){
            stop('The "bio.variables" cannot be empty when the method = twoWayAnova.')
        } else if(is.null(uv.variables)){
            stop('The "uv.variables" cannot be empty when the method = twoWayAnova.')
        } else if(length(intersect(bio.variables, uv.variables)) > 0){
            stop('Variable cannot be in both "bio.variables" and "uv.variables".')
        } else if(is.null(uv.clustering.method)){
            stop('The "uv.clustering.method" cannot be empty when the method = twoWayAnova.')
        } else if(is.null(nb.uv.clusters)){
            stop('The "nb.uv.clusters" cannot be empty when the method = twoWayAnova.')
        } else if(is.null(bio.clustering.method)){
            stop('The "bio.clustering.method" cannot be empty when the method = twoWayAnova.')
        } else if(is.null(nb.bio.clusters)){
            stop('The "nb.bio.clusters" cannot be empty when the method = twoWayAnova.')
        }
    }
    if (nb.bio.genes >= 100 | nb.bio.genes <= 0){
        stop('The "nb.bio.genes" must be a positve value 0 < nb.bio.genes =< 100.')
    }

    # assess the SummarizedExperiment object #####
    if(assess.se.obj){
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(bio.variables, uv.variables),
            remove.na = remove.na,
            verbose = verbose)
    }
    # data transformation #####
    # data normalization and transformation and regression ####
    printColoredMessage(message = '-- Data normalization and transformation :',
                        color = 'magenta',
                        verbose = verbose)
    if (!is.null(normalization)) {
        printColoredMessage(
            message = paste0('Applying the ', normalization,' on the  ', assay.name, 'data.'),
            color = 'blue',
            verbose = verbose
        )
        temp.data <- applyOtherNormalizations(
            se.obj = se.obj,
            assay.name = assay.name,
            method = normalization,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            assess.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = FALSE)
    } else if (is.null(normalization) & !is.null(apply.log) & !is.null(pseudo.count) ) {
        printColoredMessage(
            message = paste0('Applying log2 +  ', pseudo.count , '(pseudo.count) on the  ', assay.name, ' data.'),
            color = 'blue',
            verbose = verbose)
        temp.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else if (is.null(normalization) & !is.null(apply.log) & is.null(pseudo.count)) {
        printColoredMessage(
            message = paste0('Applying log2 on the  ', assay.name, ' data.'),
            color = 'blue',
            verbose = verbose)
        temp.data <- log2(assay(x = se.obj, i = assay.name))
    } else if (is.null(normalization) & is.null(apply.log)) {
        printColoredMessage(
            message = paste0('No any library size normalization and transformation ','within samples from ', assay.name, ' data.'),
            color = 'blue',
            verbose = verbose
        )
        temp.data <- assay(x = se.obj, i = assay.name)
    }
    if(!is.null(uv.variables)){
        # create all possible homogeneous groups with respect to sources of unwanted variation ####
        printColoredMessage(
            message = '-- Create all possible major groups with respect to sources of unwanted variation:',
            color = 'magenta',
            verbose = verbose)
        all.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = uv.variables,
            nb.clusters = nb.uv.clusters,
            clustering.method = uv.clustering.method,
            assess.se.obj = FALSE,
            assess.variables = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose)
    }
    # find bio genes ####
    if (method == 'pca'){
        top.loadings <- round(x = top.loadings/100 * nrow(se.obj), digits = 0)
        if (is.null(uv.variables)){
            sv.dec <- runSVD(
                x = t(temp.data),
                k = nb.pcs,
                BSPARAM = svd.bsparam,
                center = center.pca,
                scale = scale.pca)
            row.names(sv.dec$v) <- row.names(se.obj)
            hvg <- unique(unlist(lapply(
                1:nb.pcs,
                function(x) names(sort(abs(sv.dec$v[,x]), decreasing = TRUE)[1:top.loadings]))))
        } else{
            uv.groups <- unique(all.uv.groups)
            hvg <- lapply(
                uv.groups,
                function(x){
                    index.samples <- all.uv.groups == x
                    sv.dec <- runSVD(
                        x = t(temp.data[ , index.samples]),
                        k = nb.pcs,
                        BSPARAM = svd.bsparam,
                        center = center.pca,
                        scale = scale.pca)
                    row.names(sv.dec$v) <- row.names(se.obj)
                    hvg <- unique(unlist(lapply(
                        1:nb.pcs,
                        function(x) names(sort(abs(sv.dec$v[,x]), decreasing = TRUE)[1:top.loadings]))))
                })
            hvg <- unique(unlist(hvg))
        }
    } else if (method == 'mad'){
        nb.bio.genes <- round(x = nb.bio.genes/100 * nrow(se.obj), digits = 0)
        if(is.null(uv.variables)){
            hvg <- matrixStats::rowMads(x = temp.data)
            names(hvg) <- row.names(se.obj)
            hvg <- names(sort(hvg, decreasing = TRUE)[1:nb.bio.genes])
        } else{
            uv.groups <- unique(all.uv.groups)
            hvg <- lapply(
                uv.groups,
                function(x){
                    index.samples <- all.uv.groups == x
                    hvg <- matrixStats::rowMads(x = temp.data[ , index.samples])
                    names(hvg) <- row.names(se.obj)
                    hvg <- names(sort(hvg, decreasing = TRUE)[1:nb.bio.genes])
                    hvg
                })
            hvg <- unique(unlist(hvg))
        }
    } else if (method == 'vst'){
        if(is.null(uv.variables)){
            hvg <- Seurat::FindVariableFeatures(
                object = temp.data,
                selection.method = 'vst',
                verbose = FALSE
            )
            hvg <- hvg$vst.variance.standardized
            names(hvg) <- row.names(se.obj)
            hvg <- names(sort(hvg, decreasing = TRUE)[1:nb.bio.genes])
        } else{
            uv.groups <- unique(all.uv.groups)
            hvg <- lapply(
                uv.groups,
                function(x){
                    index.samples <- all.uv.groups == x
                    hvg <- Seurat::FindVariableFeatures(
                        object = temp.data[ , index.samples],
                        selection.method = 'vst',
                        verbose = FALSE
                    )
                    hvg <- hvg$vst.variance.standardized
                    names(hvg) <- row.names(se.obj)
                    hvg <- names(sort(hvg, decreasing = TRUE)[1:nb.bio.genes])
                })
            hvg <- unique(unlist(hvg))
        }
    } else if (method == 'twoWayAnova'){
        printColoredMessage(
            message = '-- Perform two_way ANOVA:',
            color = 'magenta',
            verbose = verbose)
        # create all possible homogeneous biological groups ####
        printColoredMessage(
            message = '-- Create all possible major homogeneous biological groups:',
            color = 'magenta',
            verbose = verbose)
        all.bio.groups <- createHomogeneousBioGroups(
            se.obj = se.obj,
            bio.variables = bio.variables,
            nb.clusters = nb.bio.clusters,
            clustering.method = bio.clustering.method,
            assess.se.obj = FALSE,
            assess.variables = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose)
        printColoredMessage(
            message = paste0('This is between all individual gene expression and considering both biological and UV variables created above as factors.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
        all.aov <- aov(t(expr.data) ~ all.bio.groups + all.uv.groups)
        all.aov <- summary(all.aov)
        all.aov <- as.data.frame(
            t(sapply(c(1:nrow(se.obj)),
                     function(x) all.aov[[x]]$`F value`[1:2]))
        )
        colnames(all.aov) <- c('Biology', 'UV')
        row.names(all.aov) <- row.names(se.obj)
        all.aov <- all.aov[order(all.aov$Biology, decreasing = TRUE) , ]
        bio.genes <- row.names(all.aov)[1:500]

        printColoredMessage(
            message = paste0(sum(bio.genes), ' genes are selected as negative control genes.'),
            color = 'blue',
            verbose = verbose)
    }
    # return data ####
    # saving the output ####
    printColoredMessage(message = '------------The findKnn function finished.',
                        color = 'white',
                        verbose = verbose)
    if (save.se.obj) {
        se.obj@metadata[['HVG']][[paste0('Method:',method, '||','data:', assay.name)]] <- hvg
        return(se.obj)
    } else{
        return(hvg)
    }
}



