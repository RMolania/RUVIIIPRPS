#' is used to compute the average Silhouette coefficient width using principal components of assays in a SummarizedExperiment
#' object given a categorical variable.
#'
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the Silhouette coefficients. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param dist.measure A character string indicating which method
#' is to be used for the differential analysis: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
#' By default 'euclidean' will be selected.
#' @param fast.pca Logical. Indicates whether to use the PCA calculated using a specific number of PCs instead of the full range
#' to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to be calculated for the fast pca process, by default is set to 3.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param plot.output Logical. Indicates whether to plot the Silhouette coefficients, by default it is set to FALSE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param apply.round Logical. Indicates whether to round the ARI results,by default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the computed silhouette
#' on the categorical variable.
#'
#'
#' @importFrom SummarizedExperiment assays assay
#' @importFrom stats dist
#' @importFrom cluster silhouette
#' @import ggplot2
#' @export

computeSilhouette <- function(
        se.obj,
        assay.names = 'All',
        variable,
        dist.measure = 'euclidian',
        fast.pca = TRUE,
        nb.pcs = 3,
        save.se.obj = TRUE,
        plot.output = FALSE,
        assess.se.obj = TRUE,
        apply.round = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The computeSilhouette function starts:',
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
        stop('The number of PCs (nb.pcs) must be specified.')
    }
    if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0('The ', variable, 'must be a categorical variable'))
    }
    if (!dist.measure %in% c('euclidian',
                              'maximum',
                              'manhattan',
                              'canberra',
                              'binary',
                              'minkowski')) {
        stop("The dist.measure should be one of the: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.")
    }

    # find assays ####
    if (length(assay.names) == 1 && assay.names == 'All') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- as.factor(unlist(assay.names))

    # assess the se.obj ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = 'both',
            verbose = verbose)
    }

    # silhouette coefficients on all assays ####
    printColoredMessage(
        message = '-- Compute silhouette coefficient:',
        color = 'magenta',
        verbose = verbose)
    sil.coef <- lapply(
        levels(assay.names),
        function(x) {
            if (fast.pca) {
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To compute the Silhouette coefficient, the fast PCA must be computed first on the assay ', x, '.' )
                pca.data <- se.obj@metadata[['metric']][[x]][['fastPCA']]$svd$u
            } else {
                if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To compute the Silhouette coefficient, the PCA must be computed first on the assay ', x, ' .')
                pca.data <- se.obj@metadata[['metric']][[x]][['PCA']]$svd$u
            }
            d.matrix <- as.matrix(dist(pca.data[, seq_len(nb.pcs)], method = dist.measure))
            avg.width <- summary(silhouette(as.numeric(as.factor(se.obj@colData[, variable])), d.matrix))$avg.width
            return(avg.width)
        })
    names(sil.coef) <- levels(assay.names)
    # save the results ####
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = '-- Save the silhouette coefficients to the metadata of the SummarizedExperiment object.',
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
            if (!paste0('sil.', dist.measure) %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][[paste0('sil.', dist.measure)]] <- list()
            }
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][[paste0('sil.', dist.measure)]][[variable]] <- sil.coef[[x]]
        }
        printColoredMessage(
            message = 'The silhouette coefficients for individual assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)

        ## Plot and save the plot into se.obj@metadata$plot
        if (plot.output == TRUE) {
            printColoredMessage(
                message = '-- Plot the Silhouette coefficients:',
                color = 'magenta',
                verbose = verbose
            )
            se.obj = plotMetric(
                se.obj,
                assay.names = assay.names,
                metric = 'sil',
                variable = variable,
                verbose = verbose)
            printColoredMessage(
                message = '-- A Plot of the silhouette coefficients is save metadata@plot.',
                color = 'magenta',
                verbose = verbose
            )
        }
        printColoredMessage(message = '------------The computeSilhouette function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)

        ## Return only the correlation result
    } else if (save.se.obj == FALSE) {
        printColoredMessage(message = '------------The computeSilhouette function finished.',
                            color = 'white',
                            verbose = verbose)
        return(sil = sil.coef)
    }
}
