#' compute the average Silhouette coefficient.

#' @author Ramyar Molania

#' @description
#' This function calculates the mean Silhouette coefficient for categorical variables such as sample subtypes, batches,
#' etc., and a distance matrix based on the principal components of assay(s).

#' @details
#' We use silhouette coefficient analysis to assess the separation of biological populations and batch effects. The
#' silhouette function uses Euclidean distance to calculate both the similarity between one patient and the other patients
#' in each cluster and the separation between patients in different clusters. A better normalization method will lead to
#' higher and lower silhouette coefficients for biological and batch labels, respectively.

#' @param se.obj A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the Silhouette coefficients. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param dist.measure A character string indicating which method
#' is to be used for the differential analysis: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
#' By default 'euclidean' will be selected.
#' @param fast.pca Logical. Indicates whether to use the PCA calculated using a specific number of PCs instead of the
#' full range to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to be calculated for the fast pca process, by default is set to 3.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the
#' execution of the functions, by default it is set to TRUE.

#' @return SummarizedExperiment A SummarizedExperiment object containing the computed silhouette
#' on the categorical variable.

#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @importFrom SummarizedExperiment assays assay
#' @importFrom stats dist
#' @importFrom cluster silhouette
#' @import ggplot2
#' @export

computeSilhouette <- function(
        se.obj,
        assay.names = 'all',
        variable,
        dist.measure = 'euclidian',
        fast.pca = TRUE,
        nb.pcs = 3,
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The computeSilhouette function starts:',
                        color = 'white',
                        verbose = verbose)
    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (is.null(variable)){
        stop('The "variable" cannot be empty.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop('The "variable" must be categorical variable.')
    }
    if(fast.pca){
        if (is.null(nb.pcs)) {
            stop('The number of PCs (nb.pcs) must be specified.')
        }
    }
    if (!dist.measure %in% c('euclidian',
                              'maximum',
                              'manhattan',
                              'canberra',
                              'binary',
                              'minkowski')) {
        stop("The dist.measure should be one of the: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.")
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- as.factor(unlist(assay.names))
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # silhouette coefficients on all assays ####
    printColoredMessage(
        message = '-- Compute silhouette coefficient:',
        color = 'magenta',
        verbose = verbose)
    sil.coef <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('Compute silhouette coefficient for ', x, ' data.'),
                color = 'blue',
                verbose = verbose)

            printColoredMessage(
                message = paste0('-Obtain the first ', nb.pcs,' PCs.'),
                color = 'blue',
                verbose = verbose)
            if (fast.pca) {
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To compute the Silhouette coefficient, the fast PCA must be computed first on the assay ', x, '.' )
                pca.data <- se.obj@metadata[['metric']][[x]][['fastPCA']]$svd$u
            } else {
                if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To compute the Silhouette coefficient, the PCA must be computed first on the assay ', x, ' .')
                pca.data <- se.obj@metadata[['metric']][[x]][['PCA']]$svd$u
            }
            if(ncol(pca.data) < nb.pcs){
                printColoredMessage(
                    message = paste0('The number of PCs of the assay', x, 'are ', ncol(pca.data), '.'),
                    color = 'blue',
                    verbose = verbose)
                stop(paste0('The number of PCs of the assay ', x, ' are less than', nb.pcs, '.',
                            'Re-run the computePCA function with nb.pcs = ', nb.pcs, '.'))
            }
            pca.data <- pca.data[ , seq_len(nb.pcs)]
            if(!all.equal(row.names(pca.data), colnames(se.obj))){
                stop('The column names of the SummarizedExperiment object is not the same as row names of the PCA data.')
            }
            printColoredMessage(
                message = '-Calculate the distance matrix on the PCs.',
                color = 'blue',
                verbose = verbose)
            d.matrix <- as.matrix(dist(pca.data[, seq_len(nb.pcs)], method = dist.measure))
            printColoredMessage(
                message = '-Calculate the average Silhouette coefficient.',
                color = 'blue',
                verbose = verbose)
            avg.width <- summary(silhouette(as.numeric(as.factor(se.obj@colData[, variable])), d.matrix))$avg.width
            return(avg.width)
        })
    names(sil.coef) <- levels(assay.names)
    # save the results ####
    printColoredMessage(
        message = '-- Save all the Silhouette:',
        color = 'magenta',
        verbose = verbose)
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = '-- Save the silhouette coefficients to the metadata of the SummarizedExperiment object.',
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
            if (!'silhouette' %in% names(se.obj@metadata[['metric']][[x]] )) {
                se.obj@metadata[['metric']][[x]][['silhouette']] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!paste0('sil.', dist.measure) %in% names(se.obj@metadata[['metric']][[x]][['silhouette']])) {
                se.obj@metadata[['metric']][[x]][['silhouette']][[paste0('sil.', dist.measure)]] <- list()
            }
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['silhouette']][[paste0('sil.', dist.measure)]][[variable]]$silhouette <- sil.coef[[x]]
        }
        printColoredMessage(
            message = 'The silhouette coefficients for individual assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)

        printColoredMessage(message = '------------The computeSilhouette function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)

        ## return only the correlation result ####
    } else if (save.se.obj == FALSE) {
        printColoredMessage(message = '------------The computeSilhouette function finished.',
                            color = 'white',
                            verbose = verbose)
        return(sil.coef = sil.coef)
    }
}
