#' compute the linear regression.

#' @author Ramyar Molania

#' @description
#' This functions calculates the linear regression between the the first cumulative PCs of the gene expression (assays)
#' of a SummarizedExperiment object and a continuous variable (i.e. library size).

#' @details
#' R2 values of fitted linear models are used to quantity the strength of the (linear) relationships between a single
#' quantitative source of unwanted variation, such as sample (log) library size or tumor purity, and global sample
#' summary statistics, such as the first k PCs (1 ≤ k ≤ 10).

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Symbol. A symbol or list of symbols for the selection of the name(s) of the assay(s) in the
#' SummarizedExperiment object to calculate regression analysis. The default is "all, which indicates all
#' the assays of the SummarizedExperiment object will be selected.
#' @param variable Symbol. Indicates a name of the column in the sample annotation of the SummarizedExperiment object.
#' The variable must be a continuous variable.
#' @param fast.pca Logical. Indicates whether to use the computed fast PCA or PCA results computed by the computePCA function. The
#' default is 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to use to plot the vector correlation. The default is 10.
#' @param save.se.obj Logical. Indicates whether to save the vector correlation plots to the meta data of the
#' SummarizedExperiment object or to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return A SummarizedExperiment object containing the computed regression for
#' the continuous variable and if requested the associated plot.

#' @importFrom stats lm var
#' @import ggplot2
#' @export

computePCVariableRegression <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        nb.pcs = 10,
        save.se.obj = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The computePCVariableRegression function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (!is.vector(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s).')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    } else if (length(variable) > 1){
        stop('The "variable" must contain only one variable.')
    } else if (!variable %in% colnames(se.obj@colData)){
        stop('The "variable" cannot be found in the SummarizedExperiment object.')
    } else if (var(se.obj[[variable]]) == 0) {
        stop('The "variable" must have some variation.')
    } else if (!class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop('The "variable" must be a continuous varible.')
    }
    if (sum(is.na(se.obj@colData[[variable]])) > 0){
        stop(paste0('The "', variable, '" contains NA.',
                    ' Run the checkSeObj function with "remove.na = both"',
                    ', then "computePCA"-->"computePCVariableRegression".'))
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names, levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Compute the regression on all assays ####
    printColoredMessage(
        message = '-- Compute regression:',
        color = 'magenta',
        verbose = verbose)
    all.r.squared <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('-- Compute regression for the', x , ' assay:'),
                color = 'magenta',
                verbose = verbose)
            printColoredMessage(
                message = paste0('-Obtain the first ', nb.pcs, ' PCs.'),
                color = 'blue',
                verbose = verbose)
            if (fast.pca) {
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To compute the regression,the fast PCA must be computed first on the assay ', x, ' .')
                pca.data <- se.obj@metadata[['metric']][[x]][['fastPCA']]$svd$u[colnames(se.obj),]
            } else {
                if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To compute the regression, the PCA must be computed first on the assay ', x, ' .')
                pca.data <- se.obj@metadata[['metric']][[x]][['PCA']]$svd$u[colnames(se.obj), ]
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
                message = '-Compute the R squared of the regression analysis.',
                color = 'blue',
                verbose = verbose)
            r.squared <- sapply(
                1:nb.pcs,
                function(y) lm.ls <- summary(lm(se.obj@colData[, variable] ~ pca.data[, 1:y]))$r.squared)
            return(r.squared)
        })
    names(all.r.squared) <- levels(assay.names)

    # save the results ####
    printColoredMessage(
        message = '-- Save all the regression r squared:',
        color = 'magenta',
        verbose = verbose)
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = '-Save the regression r squared to the metadata of the SummarizedExperiment object.',
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
            if (!'pcs.lm' %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['pcs.lm']] <- list()
            }
            if (!'rseq' %in% names(se.obj@metadata[['metric']][[x]][[variable]])) {
                se.obj@metadata[['metric']][[x]][['pcs.lm']][[variable]][['rseq']] <- list()
            }
            ## Check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['pcs.lm']][[variable]][['rseq']] <- all.r.squared[[x]]
        }
        printColoredMessage(
            message = 'The regression results for individal assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The computePCVariableRegression function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
        ## save the  regression r squared as a list ####
    } else if (save.se.obj == FALSE) {
        printColoredMessage(
            message = 'The regression r squared valuse for individual assay(s) are outputed as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The computePCVariableRegression function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.r.squared = all.r.squared)
    }
}
