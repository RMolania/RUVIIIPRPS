#' compute the vector correlation.

#' @author Ramyar Molania

#' @description
#' This function calculates the the vector correlation between the first cumulative PCs of the gene expression (assay)
#' of a SummarizedExperiment object and a categorical variable (i.e. batch). Then, the functions generates a line-dot plot
#' between the first cumulative PCs and the correlation coefficient to see the relationship between different PCs with the
#' variable. An ideal normalization should results a low correlation with unwanted variation variables and high correlation
#' with known biology.

#' @details
#' We use the Rozeboom squared vector correlation to quantify the strength of (linear) relationships between two sets
#' of variables, such as the first k PCs (that is 1 ≤ k ≤ 10) and dummy variables representing time, batches, plates and
#' biological variables. Not only does this quantity summarize the full set of canonical correlations, but it also reduces
#' to the familiar R2 from multiple regression when one of the variable sets contains just one element.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to calculate RLE data, medians and interquartiles. The default is "all, which indicates all
#' the assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. Indicates a name of the column in the sample annotation of the SummarizedExperiment object.
#' The variable must be a categorical variable.
#' @param fast.pca Logical. Indicates whether to use the fast PCA or PCA results computed by the computePCA function. The
#' default is 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to use to calculate the vector correlation. The default is 10. This
#' number cannot be bigger that number of PCs calculated by the computePCA function.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param save.se.obj Logical. Indicates whether to save the vector correlation plots in the metadata of the SummarizedExperiment
#' object or to output the results as list. By default it is set to 'TRUE'.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object or a list that contains the vector correlation plots of individual assay(s) for
#' the categorical variable.

#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023

#' @importFrom SummarizedExperiment assays assay
#' @importFrom fastDummies dummy_cols
#' @importFrom stats cancor
#' @import ggplot2
#' @export

computePCVariableCorrelation <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        nb.pcs = 10,
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The computePCVariableCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (is.list(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s).')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    } else if(length(variable) > 1){
        stop('The "variable" must be the name of a variable.')
    } else if (!variable %in% colnames(se.obj@colData)){
        stop('The "variable" cannot be found in the SummarizedExperiment object.')
    } else if (class(se.obj@colData[[variable]]) %in% c('numeric', 'integer')) {
        stop(paste0('The "', variable, '" must be a categorical varible.'))
    } else if (length(unique(se.obj@colData[, variable])) < 2) {
        stop('The "variable" must have at least two levels.')
    } else if (sum(is.na(se.obj@colData[[variable]])) > 0){
        stop(paste0('The "', variable, '" contains NA. Re-run the computePCA with "remove.na = both"'))
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names, levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # create dummy variables ####
    printColoredMessage(
        message =  '-- Create dummy variables:',
        color = 'magenta',
        verbose = verbose)
    catvar.dummies <- fastDummies::dummy_cols(se.obj@colData[[variable]])
    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
    printColoredMessage(
        message =  paste0('-A design matrix with ', ncol(catvar.dummies), ' columns is created.'),
        color = 'blue',
        verbose = verbose)

    # compute vector correlation ####
    printColoredMessage(
        message =  '-- Compute vector correlation:',
        color = 'magenta',
        verbose = verbose)
    all.vec.corr <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('Compute vector correlation for ', x, ' data:'),
                color = 'blue',
                verbose = verbose)
            printColoredMessage(
                message = paste0('-Obtain the first ', nb.pcs, ' PCs.'),
                color = 'blue',
                verbose = verbose)
            if (fast.pca) {
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To compute the regression, the fast PCA must be computed first on the assay ', x, ' .' )
                pca.data <- se.obj@metadata[['metric']][[x]][['fastPCA']]$svd$u
            } else {
                if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To compute the regression, the PCA must be computed first on the assay ', x,' .')
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
            printColoredMessage(
                message = '-Calculate vector correlation.',
                color = 'blue',
                verbose = verbose)
            vector.corr <- sapply(
                1:nb.pcs,
                function(y) {
                    cca <- cancor(x = pca.data[, 1:y, drop = FALSE], y = catvar.dummies)
                    1 - prod(1 - cca$cor ^ 2)
                })
            return(vector.corr)
        })
    names(all.vec.corr) <- levels(assay.names)

    # add results to the SummarizedExperiment object ####
    printColoredMessage(
        message = '-- Save all the vector correlation results:',
        color = 'magenta',
        verbose = verbose)
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = '-Save the vector correlation results to the metadata of the SummarizedExperiment object.',
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
            if (!'pcs.vect.corr' %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['pcs.vect.corr']] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!variable %in% names(se.obj@metadata[['metric']][[x]][['pcs.vect.corr']])) {
                se.obj@metadata[['metric']][[x]][['pcs.vect.corr']][[variable]] <- list()
            }
            se.obj@metadata[['metric']][[x]][['pcs.vect.corr']][[variable]][['corrs']] <- list()
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['pcs.vect.corr']][[variable]][['corrs']] <- all.vec.corr[[x]]
        }
        printColoredMessage(
            message = 'The vector correlation for the individal assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The computePCVariableCorrelations function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)

        ## Return only the correlation result
    } else if (save.se.obj == FALSE) {
        printColoredMessage(
            message = 'The vector correlation for the individal assay(s) are outputed as list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The computePCVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.vec.corr = all.vec.corr)
    }
}
