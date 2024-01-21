#' is used to compute the linear regression.

#' @author Ramyar Molania

#' @description
#' This functions calculates the linear regression between the the first cumulative PCs of the gene expression (assays)
#' of a SummarizedExperiment object and a continuous variable (i.e. library size)

#' @details
#' R2 values of fitted linear models are used to quantity the strength of the (linear) relationships between a single
#' quantitative source of unwanted variation, such as sample (log) library size or tumor purity, and global sample
#' summary statistics, such as the first k PCs (1 ≤ k ≤ 10).

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the regression. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a continuous variable such as
#' library size from colData(se.obj).
#' @param fast.pca Logical. Indicates whether to use the PCA calculated using a specific number of PCs instead of the
#' full range to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of few first cumulative PCs, by default is set to 10.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na Symbol. To remove NA or missing values from the assay(s) or variable or both. The options are
#' "assays", "sample.annotation, "both" or "none. "See the checkSeObj function for more details.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the
#' execution of the functions, by default it is set to TRUE.

#' @return A SummarizedExperiment object containing the computed regression for
#' the continuous variable and if requested the associated plot.

#' @importFrom stats lm var
#' @import ggplot2
#' @export

## deal with PCA and remove NA from variable
computePCVariableRegression <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        nb.pcs = 10,
        save.se.obj = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'both',
        verbose = TRUE
) {
    printColoredMessage(message = '------------The computePCVariableRegression function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    } else if (!variable %in% colnames(se.obj@colData)){
        stop('The "variable" cannot be found in the SummarizedExperiment object.')
    } else if (var(se.obj[[variable]]) == 0) {
        stop('The "variable" must have some variation.')
    } else if (!class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop('The "variable" must be a continuous varible.')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- as.factor(unlist(assay.names))
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # assess the SummarizedExperiment ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = remove.na,
            verbose = verbose)
    }
    ### Compute the regression on all assays
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
            message = '-Save the regression r squared as a list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The computePCVariableRegression function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.r.squared = all.r.squared)
    }
}
