#' is used to compute the adjusted rand index (ARI).
#'
#'
#' @description
#' This functions computes the adjusted rand index using the first PCs of the assays in a SummarizedExperiment object
#' given a categorical variable.
#'
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to compute PCA. By default all the assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. Indicates the column name in the SummarizedExperiment object that contains a categorical
#' variable such as sample types or batches.
#' @param fast.pca Logical. Indicates whether to use the PCA calculated using a specific number of PCs instead of the
#' full range to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to be calculated for the fast pca process, by default is set to 3.
#' @param clustering.method Symbol. Indicates which clustering methods to be used. The function provides the mclust or
#' hclust methods. The default is hclust.
#' @param hclust.method Symbol. Indicate the agglomeration method to be used for the hclust methid. This should be
#' one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or
#' "centroid" (= UPGMC). See the hclust function for more details.
#' @param dist.measure Symbol. Indicates the distance measure to be used in the dist function. This must be one of
#' "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". See the dist function for more details.
#' @param plot.output Logical. Indicates whether to plot the ARI, by default it is set to FALSE.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object, by default it is
#' set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the
#' execution of the functions, by default it is set to TRUE.
#'
#'
#' @return A SummarizedExperiment object or a list that containing the computed ARI on the categorical variable.
#'
#'
#' @author Ramyar Molania
#'
#'
#' @importFrom SummarizedExperiment assays assay
#' @importFrom mclust mclustBIC Mclust adjustedRandIndex
#' @importFrom stats cutree hclust dist
#' @import ggplot2
#' @export

computeARI <- function(
        se.obj,
        assay.names = 'All',
        variable,
        fast.pca = TRUE,
        nb.pcs = 3,
        clustering.method = 'hclust',
        hclust.method = 'complete',
        dist.measure = 'euclidian',
        plot.output = FALSE,
        save.se.obj = TRUE,
        assess.se.obj = TRUE,
        verbose = TRUE
        ) {
    printColoredMessage(message = '------------The computeARI function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (length(assay.names) == 1 & assay.names != 'All') {
        if (!assay.names %in% names(assays(se.obj)))
            stop('The assay name cannot be found in the SummarizedExperiment object.')
    }
    if (length(assay.names) > 1) {
        if (sum(!assay.names %in% names(assays(se.obj))) > 0 )
            stop('The assay names cannot be found in the SummarizedExperiment object.')
    }
    if (is.null(assay.names)) {
        stop('The assay.names cannot be empty.')
    }
    if (is.null(nb.pcs)) {
        stop('The number of PCs (nb.pcs) must be specified')
    }
    if (is.null(variable)) {
        stop('The variable cannot be empty.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop('A categorical variable must be specified as variable.')
    }
    if(is.null(clustering.method)){
        stop('The clustering.method cannot be empty.')
    } else if(!clustering.method %in% c('mclust', 'hclust')){
        stop('The clustering.method method should be one of "mclust" or "hclust" .')
    }
    if(clustering.method == 'hclust'){
        if(is.null(dist.measure)){
            stop('The dist.measure cannot be empty when the clustering.method is hclust.')
        } else if (!dist.measure %in% c('euclidian',
                                       'maximum',
                                       'manhattan',
                                       'canberra',
                                       'binary',
                                       'minkowski')) {
            stop("The dist.measure should be one of the:'euclidean','maximum','manhattan','canberra','binary'or'minkowski'.")
        }
        if(is.null(hclust.method)){
            stop('The hclust.method cannot be when the clustering.method is hclust.')
        } else if (!hclust.method %in% c('complete',
                                       'ward.D',
                                       'ward.D2',
                                       'single',
                                       'average',
                                       'mcquitty',
                                       'median',
                                       'centroid')) {
            stop("The hclust.method should be one of the:'complete','ward.D','ward.D2','single','average', 'mcquitty', 'median' or 'centroid'.")
        }
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'All') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else  assay.names <- as.factor(unlist(assay.names))

    # assess the SummarizedExperiment ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = 'both',
            verbose = verbose)
    }
    # ari on all assays ####
    printColoredMessage(
        message = '-- Compute ARI: ',
        color = 'magenta',
        verbose = verbose)
    all.ari <- lapply(
        levels(assay.names),
        function(x) {
            if (fast.pca) {
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]])) {
                    stop('To compute the ARI the fast PCA must be computed first on the assay ', x, ' .')
                }
                ### Silhouette on PCs and categorical variable
                printColoredMessage(
                    message = paste0(
                        'Compute ARI using on the first ',
                        nb.pcs, '
                        PCs of ',
                        x,
                        ' data and the ',
                        variable,
                        ' variable.'),
                    color = 'blue',
                    verbose = verbose)
                pca.data <- se.obj@metadata[['metric']][[x]][['fastPCA']]$svd.dec$u[colnames(se.obj), ]
            } else {
                pca.data <- se.obj@metadata[['metric']][[x]][['PCA']]$svd.dec$u[colnames(se.obj), ]
            }
            if(clustering.method == 'mclust'){
                bic <- mclustBIC(data = pca.data)
                mod <- Mclust(data = pca.data, x = bic, G = length(unique(se.obj@colData[, variable])) )
                ari <- adjustedRandIndex(mod$classification, se.obj@colData[, variable])
            } else {
                clusters <- cutree(
                    tree = hclust(d = dist(x = pca.data, method = dist.measure), method = hclust.method),
                    k = length(unique(se.obj@colData[, variable])))
                ari <- adjustedRandIndex(clusters, se.obj@colData[, variable])
            }
            return(all.ari)
        })
    names(all.ari) <- levels(assay.names)
    # save the results ####
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(message = '-- Save the ARI results to the metadata of the SummarizedExperiment object:',
                            color = 'magenta',
                            verbose = verbose)

        for (x in levels(assay.names)) {
            ## check if metadata metric already exist
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!'metric' %in% names(se.obj@metadata)) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!x %in% names(se.obj@metadata[['metric']])) {
                se.obj@metadata[['metric']][[x]] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!'ari' %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['ari']] <- list()
            }
            if(clustering.method == 'mclust'){
                out.put.name <- 'ari.mclust'
            } else out.put.name <- paste0('ari.hclust.', hclust.method, '.', dist.measure)
            se.obj@metadata[['metric']][[x]][[out.put.name]][[variable]] <- all.ari[[x]]
        }
        printColoredMessage('The ARI results of induvial assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)
        ## Plot and save the plot into se.obj@metadata$plot
        if (plot.output == TRUE) {
            printColoredMessage(
                message = '-- Plot of the ARI results:',
                color = 'magenta',
                verbose = verbose
            )
            printColoredMessage(
                message = '-- A plot of the ARI results are saved to the metadata of the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose
            )
            se.obj <- plotMetric(
                se.obj,
                assay.names = assay.names,
                metric = 'ari',
                variable = variable,
                verbose = verbose)
        }
        printColoredMessage(message = '------------The computeARI function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)

        ## Return only the correlation result
    } else if (save.se.obj == FALSE) {
        printColoredMessage(message = '------------The computeARI function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.ari = all.ari)
    }
}
