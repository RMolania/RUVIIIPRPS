#' is used to compute the Silhouette coefficient from the first PC of a SummarizedExperiment class
#' object given a categorical variable.
#'
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the Silhouette coefficients. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param fast.pca TO BE DEFINED.
#' @param nb.pcs TO BE DEFINED.
#' @param save.se.obj Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param plot.output Indicates whether to plot the Silhouette coefficients, by default it is set to FALSE.
#' @param assess.se.obj Indicates whether to assess the SummarizedExperiment class object.
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the computed silhouette
#' on the categorical variable.
#' @importFrom stats dist
#' @importFrom cluster silhouette
#' @import ggplot2
#' @export

## deal with PCA and remove NA from variable

computeSilhouette<-function(
        se.obj,
        assay.names = 'All',
        variable,
        fast.pca = TRUE,
        nb.pcs = 3,
        save.se.obj = TRUE,
        plot.output=FALSE,
        assess.se.obj = TRUE,
        verbose = TRUE
){

    printColoredMessage(message = '------------The computeSilhouette function starts:',
                        color = 'white',
                        verbose = verbose)

    ### check the inputs
    if (is.null(nb.pcs)) {
        stop('To compute the Silhouette coefficient, the number of PCs (nb.pcs) must be specified.')
    } else if (is.null(assay.names)) {
        stop('Please provide at least an assay name.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
            stop(paste0('The ', variable,', is a numeric, but this should a categorical variable'))
    }

    ### Assess the se.obj
    rm.pca=FALSE
    if (assess.se.obj){
        se.obj.orig=se.obj
        se.obj <- checkSeObj(se.obj = se.obj,
                             assay.names = assay.names,
                             variables = variable,
                             remove.na = 'both',
                             verbose = verbose)
        if (ncol(se.obj)!=ncol(se.obj.orig)){
            rm.pca=TRUE
            keep.samples <- complete.cases(colData(se.obj)[, variable, drop = FALSE])
        }
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
        # Silhouette coefficients on all assays
        silCoef <- lapply(
            levels(assay.names),
            function(x){
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]])) {
                    stop('To compute the Silhouette coefficient,
                         the fast PCA must be computed first on the assay ', x, ' .')
                }
                ### Silhouette on PCs and categorical variable
                printColoredMessage(message = paste0(
                    '### Computing Silhouette coefficient based on PCs and the ',
                    variable,
                    ' variable.'
                ),
                color = 'magenta',
                verbose = verbose)
                pca_x <- se.obj@metadata[['metric']][[x]][['fastPCA']]$sing.val$u
                    if (isTRUE(rm.pca)){
                        pca_x=pca_x[keep.samples]
                    }
                d.matrix <- as.matrix(dist(pca_x[, seq_len(nb.pcs)]))
                avg=summary(silhouette(
                    as.numeric(as.factor(var)),
                    d.matrix))$avg.width
                return(avg)
        })
        names(silCoef) <- levels(assay.names)
    } else{
        # Silhouette coefficients on all assays
        silCoef <- lapply(
            levels(assay.names),
            function(x){
                if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]])) {
                    stop('To compute the Silhouette coefficient,
                         the PCA must be computed first on the assay ', x, ' .')
                }
                ### Silhouette on PCs and categorical variable
                printColoredMessage(message = paste0(
                    '### Computing Silhouette coefficient based on PCs and the ',
                    variable,
                    ' variable.'
                ),
                color = 'magenta',
                verbose = verbose)
                pca_x <- se.obj@metadata[['metric']][[x]][['PCA']]
                d.matrix <- as.matrix(dist(pca_x[, seq_len(nb.pcs)]))
                avg=summary(silhouette(
                    as.numeric(as.factor(var)),
                    d.matrix))$avg.width
                return(avg)
            })
        names(silCoef) <- levels(assay.names)
    }

    ### Add results to the SummarizedExperiment object
    if(save.se.obj == TRUE){
        printColoredMessage(message= '### Saving the Silhouette coefficients to the metadata of the SummarizedExperiment object.',
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
                se.obj@metadata[['metric']][[x]][['sil']] <- list()
            }
            ## Check if metadata metric already exist for this assay, this metric and this variable
            if(!variable %in% names(se.obj@metadata[['metric']][[x]][['sil']])  ) {
                se.obj@metadata[['metric']][[x]][['sil']][[variable]] <- silCoef[[x]]
            }
        }
        printColoredMessage(message= paste0(
            'The Silhouette coefficients are saved to metadata@',
            x,
            '$sil$',
            variable, '.'),
            color = 'blue',
            verbose = verbose)


        ## Plot and save the plot into se.obj@metadata$plot
        if (plot.output==TRUE) {
            printColoredMessage(message= '### Plotting and Saving the Silhouette coefficients to the metadata of the SummarizedExperiment object.',
                                color = 'magenta',
                                verbose = verbose)

            se.obj=plotMetric(se.obj,
                              assay.names =assay.names,
                              metric='sil',
                              variable=variable,
                              verbose=verbose)
        }
        return(se.obj = se.obj)

        ## Return only the correlation result
    } else if(save.se.obj == FALSE){
        return(sil=silCoef)
    }

    printColoredMessage(message = '------------The computeSilhouette function finished.',
                        color = 'white',
                        verbose = verbose)

}

