#' is used to compute the vector correlation.

#' @author Ramyar Molania

#' @description
#' This function calculates the the vector correlation between the first cumulative PCs of the gene expression (assay)
#' of a SummarizedExperiment object and a categorical variable (i.e. batch).

#' @details
#' We used the Rozeboom squared vector correlation60 to quantify the strength of (linear) relationships between two sets
#' of variables, such as the first k PCs (that is 1 ≤ k ≤ 10) and dummy variables representing time, batches, plates and
#' biological variables. Not only does this quantity summarize the full set of canonical correlations, but it also reduces
#' to the familiar R2 from multiple regression (see below) when one of the variable sets contains just one element.
#'

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the correlation. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param fast.pca Logical. Indicates whether to use the PCA calculated using a specific number of PCs instead of the
#' full range to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of few first cumulative PCs, by default is set to 10.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result. By default it is set to TRUE.
#' @param plot.output Logical. Indicates whether to plot the correlation statistics, by default it is set to TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na Symbol. To remove NA or missing values from the assay(s) or variable or both. The options are
#' "assays", "sample.annotation, "both" or "none. "See the checkSeObj function for more details.
#' @param apply.round Logical. Indicates whether to round the ARI results, by default it is set to TRUE.
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.

#' @return SummarizedExperiment A SummarizedExperiment object containing the computed correlation for
#' the continuous variable and if requested the associated plot.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom fastDummies dummy_cols
#' @importFrom stats cancor
#' @import ggplot2
#' @export

PCVariableCorrelation <- function(
        se.obj,
        assay.names = 'All',
        variable,
        fast.pca = TRUE,
        nb.pcs = 10,
        save.se.obj = TRUE,
        plot.output = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'both',
        apply.round = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The PCVariableCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the inputs ####
    if (is.null(assay.names)) {
        stop('Please provide at least an assay name.')
    } else if (is.null(variable)) {
        stop('Please provide a variable.')
    } else if (length(unique(se.obj@colData[, variable])) < 2) {
        stop(paste0('The ', variable, ', contains only one variable.'))
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0(
            'The ',
            variable,
            ', is a numeric, but this should a categorical variable'
        ))
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'All') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else  assay.names <- as.factor(unlist(assay.names))

    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = remove.na,
            verbose = verbose)
    }
    ## create dummy variables ####
    printColoredMessage(
        message =  '-- Create dummy variables:',
        color = 'magenta',
        verbose = verbose)
    catvar.dummies <- dummy_cols(se.obj@colData[, variable])
    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]

    ### Compute the correlation on all assays
    ### Regression on PCs and continous variable
    printColoredMessage(
        message =  '-- Compute vector correlation:',
        color = 'magenta',
        verbose = verbose)
    all.vec.corr <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0('Obtain the first ', nb.pcs, ' PCs of ', x, ' data.'),
                color = 'blue',
                verbose = verbose)
            if (fast.pca) {
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To compute the regression, the fast PCA must be computed first on the assay ', x, ' .' )
                pca.data <- se.obj@metadata[['metric']][[x]][['fastPCA']]$sv.dec$u
            } else {
                if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]]))
                    stop('To compute the regression, the PCA must be computed first on the assay ', x,' .')
                pca.data <- se.obj@metadata[['metric']][[x]][['PCA']]$sv.dec$u
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
                message = 'Calculate vector correlation.',
                color = 'blue',
                verbose = verbose)
            cca.pcs <- sapply(
                1:nb.pcs,
                function(y) {
                    cca <- cancor(x = pca.data[, 1:y, drop = FALSE], y = catvar.dummies)
                    1 - prod(1 - cca$cor ^ 2)
                })
            return(cca.pcs)
        })
    names(all.vec.corr) <- levels(assay.names)

    # add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = '-- Save the vector correlation results to the metadata of the SummarizedExperiment object.',
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
            if (!'pcs.vect.corr' %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['pcs.vect.corr']] <- list()
            }
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['pcs.vect.corr']][[variable]] <- all.vec.corr[[x]]
        }
        printColoredMessage(
            message = 'The vector correlation for individal assays are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)
        ## Plot and save the plot into se.obj@metadata$plot ####
        if (plot.output == TRUE) {
            printColoredMessage(
                message = '-- Plot the the vector correlation results:',
                color = 'magenta',
                verbose = verbose)
            printColoredMessage(
                message = 'A plot of the vector correlations are saved to metadata@plot.',
                color = 'blue',
                verbose = verbose)
            se.obj <- plotMetric(
                se.obj,
                assay.names = assay.names,
                metric = 'pcs.vect.corr',
                variable = variable,
                verbose = verbose)
        }
        printColoredMessage(message = '------------The PCVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)

        ## Return only the correlation result
    } else if (save.se.obj == FALSE) {
        printColoredMessage(message = '------------The PCVariableCorrelation function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.vec.corr = all.vec.corr)
    }
}
