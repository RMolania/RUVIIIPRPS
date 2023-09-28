#' is used to compute the adjusted rand index (ARI) from the first PC of a SummarizedExperiment class
#' object given a categorical variable.
#'
#' It can be used to assess how a group of biological samples
#' are distributed across batches (i.e. example subtypes vs batch), and how batches ### update
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
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param plot.output Logical. Indicates whether to plot the ARI, by default it is set to FALSE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object, by default it is set to TRUE.
#' @param apply.round Logical. Indicates whether to round the ARI results, by default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the computed ARI
#' on the categorical variable.
#' @importFrom SummarizedExperiment assays assay
#' @importFrom mclust mclustBIC Mclust adjustedRandIndex
#' @import ggplot2
#' @export
#'

computeARI <-function(
        se.obj,
        assay.names = 'All',
        variable,
        fast.pca = TRUE,
        nb.pcs = 3,
        save.se.obj = TRUE,
        plot.output=FALSE,
        assess.se.obj = TRUE,
        apply.round = TRUE,
        verbose = TRUE
){

    printColoredMessage(message = '------------The computeAri function starts:',
                        color = 'white',
                        verbose = verbose)

    ### check the inputs
    if (is.null(nb.pcs)) {
        stop('To compute the ARI, the number of PCs (nb.pcs) must be specified.')
    } else if (is.null(assay.names)) {
        stop('Please provide at least an assay name.')
    } else if (is.null(variable)) {
        stop('Please provide a variable.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0('The ', variable,', is a numeric, but this should a categorical variable'))
    }

    ### Assess the se.obj
    if (assess.se.obj){
        se.obj <- checkSeObj(se.obj = se.obj,
                             assay.names = assay.names,
                             variables = variable,
                             remove.na = 'both',
                             verbose = verbose)
    }

    ## Assays
    if(length(assay.names) == 1 && assay.names=='All'){
        assay.names=as.factor(names(assays(se.obj)))
    }else {
        assay.names=as.factor(unlist(assay.names))
    }

    ## Categorical variable
    var=se.obj@colData[, variable]
    var.label=paste0(variable)

    if (fast.pca) {
        # ARI on all assays
        ari <- lapply(
            levels(assay.names),
            function(x){
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]])) {
                    stop('To compute the ARI,
                         the fast PCA must be computed first on the assay ', x, ' .')
                }
                ### Silhouette on PCs and categorical variable
                printColoredMessage(message = paste0(
                    '### Computing ARI based on PCs and the ',
                    variable,
                    ' variable.'
                ),
                color = 'magenta',
                verbose = verbose)
                pca_x <- se.obj@metadata[['metric']][[x]][['fastPCA']]$sing.val$u[colnames(se.obj),]
                BIC <- mclustBIC(data = pca_x)
                mod <- Mclust(data = pca_x, x = BIC)
                ari=adjustedRandIndex(
                    mod$classification,
                    var)
                return(ari)
        })
        names(ari) <- levels(assay.names)
    } else{
        # ARI on all assays
        ari <- lapply(
            levels(assay.names),
            function(x){
                if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]])) {
                    stop('To compute the ARI,
                         the PCA must be computed first on the assay ', x, ' .')
                }
                ### Silhouette on PCs and categorical variable
                printColoredMessage(message = paste0(
                    '### Computing ARI based on PCs and the ',
                    variable,
                    ' variable.'
                ),
                color = 'magenta',
                verbose = verbose)
                pca_x <- se.obj@metadata[['metric']][[x]][['PCA']]$sing.val$u[colnames(se.obj),]
                BIC <- mclustBIC(data = pca_x)
                mod <- Mclust(data = pca_x, x = BIC)
                ari=adjustedRandIndex(
                    mod$classification,
                    var)
                return(ari)
        })
        names(ari) <- levels(assay.names)
    }

    ## Round the regression statistic obtained to 4 digits
    if(apply.round){
        ari[] <- lapply(ari, round,4)
    }

    ### Add results to the SummarizedExperiment object
    if(save.se.obj == TRUE){
        printColoredMessage(message= '### Saving the ARI to the metadata of the SummarizedExperiment object.',
                            color = 'magenta',
                            verbose = verbose)

        for (x in levels(assay.names)){
            ## Check if metadata metric already exist
            if(length(se.obj@metadata)==0 ) {
                se.obj@metadata[['metric']] <- list()
            }
            ## Check if metadata metric already exist for this assay
            if(!'metric' %in% names(se.obj@metadata) ) {
                se.obj@metadata[['metric']] <- list()
            }
            ## Check if metadata metric already exist for this assay
            if(!x %in% names(se.obj@metadata[['metric']]) ) {
                se.obj@metadata[['metric']][[x]] <- list()
            }
            ## Check if metadata metric already exist for this assay and this metric
            if(!'sil' %in% names(se.obj@metadata[['metric']][[x]])  ) {
                se.obj@metadata[['metric']][[x]][['ari']] <- list()
            }
            ## Check if metadata metric already exist for this assay, this metric and this variable
            if(!variable %in% names(se.obj@metadata[['metric']][[x]][['ari']])  ) {
                se.obj@metadata[['metric']][[x]][['ari']][[variable]] <- ari[[x]]
            }
        }
        printColoredMessage(message= paste0(
            'The ARI is saved to metadata@metric$',
            x,
            '$ari$',
            variable, '.'),
            color = 'blue',
            verbose = verbose)


        ## Plot and save the plot into se.obj@metadata$plot
        if (plot.output==TRUE) {
            printColoredMessage(message= '### Plotting and Saving the ARI to the metadata of the SummarizedExperiment object.',
                                color = 'magenta',
                                verbose = verbose)

            se.obj=plotMetric(se.obj,
                              assay.names =assay.names,
                              metric='ari',
                              variable=variable,
                              verbose=verbose)
        }
        return(se.obj = se.obj)

        ## Return only the correlation result
    } else if(save.se.obj == FALSE){
        return(ari=ari)
    }

    printColoredMessage(message = '------------The computeAri function finished.',
                        color = 'white',
                        verbose = verbose)


}
