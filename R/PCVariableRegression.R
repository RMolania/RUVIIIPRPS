#' is used to compute the linear regression between the the first cumulative PCs
#' of the gene expression (assay) of a SummarizedExperiment class object and a continuous variable (i.e. library size)

#' @param se.obj A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the regression. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a continuous variable such as
#' library size from colData(se.obj).
#' @param fast.pca Logical. Indicates whether to use the PCA calculated using a specific number of PCs instead of the full range
#' to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of few first cumulative PCs, by default is set to 10.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param plot.output Logical. Indicates whether to plot the regression statistics, by default it is set to TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na TO BE DEFINED.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.

#' @return SummarizedExperiment A SummarizedExperiment object containing the computed regression for
#' the continuous variable and if requested the associated plot.

#' @importFrom stats lm var
#' @import ggplot2
#' @export

## deal with PCA and remove NA from variable
PCVariableRegression <- function(
        se.obj,
        assay.names = 'All',
        variable,
        fast.pca = TRUE,
        nb.pcs = 10,
        save.se.obj = TRUE,
        plot.output = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'both',
        verbose = TRUE
) {
    printColoredMessage(message = '------------The PCVariableRegression function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be null.')
    }
    if(length(assay.names) == 1 & assay.names!= 'All'){
        if(!assay.names %in% names(assays(se.obj)) )
            stop('The assay name cannot be found in the SummarizedExperiment object.')
    }
    if(length(assay.names) > 1){
        if(!assay.names %in% names(assays(se.obj)))
            stop('The assay names cannot be found in the SummarizedExperiment object.')
    }
    if (is.null(variable)) {
        stop('The variable cannot be empty.')
    } else if (!class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop( paste0( 'The ', variable, 'this should a continuous variable'))
    } else if (is.null(nb.pcs)) {
        stop('To compute the regression, the number of PCs (nb.pcs) must be specified.')
    }
    if (var(se.obj[[variable]]) == 0) {
        stop(paste0('The ', variable, ' appears to have no variation.'))
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'All') {
        assay.names = as.factor(names(assays(se.obj)))
    } else assay.names = as.factor(unlist(assay.names))

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
        message = '-- Apply regression:',
        color = 'magenta',
        verbose = verbose)
    all.r.squared <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0(
                    '-- Compute R^2 of the regression between PCs of the ',
                    x ,
                    ' and ',
                    variable,
                    ' variable.'),
                color = 'blue',
                verbose = verbose)
            if (fast.pca) {
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To compute the regression,the fast PCA must be computed first on the assay ', x, ' .')
                pca.data <-
                    se.obj@metadata[['metric']][[x]][['fastPCA']]$svd$u[colnames(se.obj),]
            } else {
                if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To compute the regression, the PCA must be computed first on the assay ', x, ' .')
                pca.data <- se.obj@metadata[['metric']][[x]][['PCA']]$svd$u[colnames(se.obj), ]
            }
            r.squared <- sapply(
                1:nb.pcs,
                function(y) lm.ls <- summary(lm(se.obj@colData[, variable] ~ pca.data[, 1:y]))$r.squared)
            return(r.squared)
        })
    names(all.r.squared) <- levels(assay.names)

    # save the results ####
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = '-- Save the regression results to the metadata of the SummarizedExperiment object.',
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
            if (!'pcs.lm' %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['pcs.lm']] <- list()
            }
            ## Check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['pcs.lm']][[variable]] <- all.r.squared[[x]]
        }
        printColoredMessage(
            message = 'The regression results for individal assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose
            )
        ## Plot and save the plot into se.obj@metadata$plot
        if (plot.output == TRUE) {
            printColoredMessage(
                message = '-- Plot the the regression results:',
                color = 'magenta',
                verbose = verbose)
            printColoredMessage(
                message = 'A plot of the R^2 of the regression are saved to metadata@plot.',
                color = 'blue',
                verbose = verbose)
            se.obj <- plotMetric(
                se.obj,
                assay.names = assay.names,
                metric = 'pcs.lm',
                variable = variable,
                verbose = verbose)
        }
        printColoredMessage(message = '------------The PCVariableRegression function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
        ## return only the regression results ####
    } else if (save.se.obj == FALSE) {
        printColoredMessage(message = '------------The PCVariableRegression function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.r.squared = all.r.squared)
    }
}
