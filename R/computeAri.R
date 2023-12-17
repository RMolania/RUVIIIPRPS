#' is used to compute the adjusted rand index (ARI) using the first PC of the assays in a SummarizedExperiment
#' object given a categorical variable.
#'
#' It can be used to assess how a group of biological samples are distributed across batches (i.e. example subtypes vs batch), and how batches ### update
#'
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the ARI. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param fast.pca Logical. Indicates whether to use the PCA calculated using a specific number of PCs instead of the full range
#' to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to be calculated for the fast pca process, by default is set to 3.
#' @param clustering.method Symbol. Indicates the clustering methods.
#' @param hclust.method Symbol.
#' @param dist.measure Indicates the clustering methods.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param plot.output Logical. Indicates whether to plot the ARI, by default it is set to FALSE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object, by default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the computed ARI
#' on the categorical variable.
#' @importFrom SummarizedExperiment assays assay
#' @importFrom mclust mclustBIC Mclust adjustedRandIndex
#' @importFrom stats hclust dist
#' @import ggplot2
#' @export
#'

computeARI <- function(
        se.obj,
        assay.names = 'All',
        variable,
        fast.pca = TRUE,
        nb.pcs = 3,
        clustering.method = 'hclust',
        hclust.method = 'complete',
        dist.measure = 'euclidian',
        save.se.obj = TRUE,
        plot.output = FALSE,
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
        if (!assay.names %in% names(assays(se.obj)))
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
            stop('The dist.measure cannot be when the clustering.method is hclust.')
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
                pca.data <- se.obj@metadata[['metric']][[x]][['fastPCA']]$sing.val$u[colnames(se.obj), ]
            } else {
                pca.data <- se.obj@metadata[['metric']][[x]][['PCA']]$sing.val$u[colnames(se.obj), ]
            }
            if(clustring.method == 'mclust'){
                bic <- mclustBIC(data = pca.data)
                mod <- Mclust(data = pca.data, x = bic, G = length(unique(se.obj@colData[, variable])) )
                ari <- adjustedRandIndex(mod$classification, se.obj@colData[, variable])
            } else {
                clusters <- cutree(
                    tree = hclust(d = dist(x = pca.dat, method = dist.measure), method = hclust.method),
                    k = length(unique(se.obj@colData[, variable])))
                ari <- adjustedRandIndex(clusters, se.obj@colData[, variable])
            }
            return(all.ari)
        })
    names(all.ari) <- levels(assay.names)
    # save the results ####
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(message = '-- Save the ARI to the metadata of the SummarizedExperiment object.',
                            color = 'blue',
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
            se.obj@metadata[['metric']][[x]][['ari']][[variable]] <- all.ari[[x]]

        }
        printColoredMessage( message = paste0(
            'The ARI is saved to metadata@metric$',
            x,
            '$ari$',
            variable,
            '.'),
            color = 'blue',
            verbose = verbose
        )
        ## Plot and save the plot into se.obj@metadata$plot
        if (plot.output == TRUE) {
            printColoredMessage(
                message = '### Plotting and Saving the ARI to the metadata of the SummarizedExperiment object.',
                color = 'magenta',
                verbose = verbose
            )
            se.obj = plotMetric(
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
