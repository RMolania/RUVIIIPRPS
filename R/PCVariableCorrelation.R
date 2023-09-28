#' is used to compute the vector correlation between the first cumulative PCs of the gene expression (assay)
#' of a SummarizedExperiment class object and a categorical variable (i.e. batch).
#'
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the correlation. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param fast.pca Logical. Indicates whether to use the PCA calculated using a specific number of PCs instead of the full range
#' to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of few first cumulative PCs, by default is set to 10.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param plot.output Logical. Indicates whether to plot the correlation statistics, by default it is set to TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na TO BE DEFINED.
#' @param apply.round Logical. Indicates whether to round the ARI results, by default it is set to TRUE.
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the computed correlation for
#' the continuous variable and if requested the associated plot.
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom SummarizedExperiment assays assay
#' @importFrom matrixTests row_oneway_equalvar
#' @importFrom fastDummies dummy_cols
#' @importFrom wesanderson wes_palette
#' @importFrom stats cancor
#' @import ggplot2
#' @export

## deal with PCA and remove NA from variable
PCVariableCorrelation<-function(
        se.obj,
        assay.names = 'All',
        variable,
        fast.pca = TRUE,
        nb.pcs = 10,
        save.se.obj = TRUE,
        plot.output=TRUE,
        assess.se.obj = TRUE,
        remove.na = 'both',
        apply.round = TRUE,
        verbose = TRUE
){

    printColoredMessage(message = '------------The PCVariableCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)

    ### Check the inputs
    if (is.null(assay.names)) {
        stop('Please provide at least an assay name.')
    } else if (is.null(variable)) {
        stop('Please provide a variable.')
    } else if (length(unique(se.obj@colData[, variable])) < 2) {
        stop(paste0('The ', variable,', contains only one variable.'))
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0('The ', variable,', is a numeric, but this should a categorical variable'))
    }

    ### Assess the se.obj
    if (assess.se.obj){
        se.obj <- checkSeObj(se.obj = se.obj,
                             assay.names = assay.names,
                             variables = variable,
                             remove.na = remove.na,
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

    catvar.dummies <- dummy_cols(var)
    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]

    if (fast.pca) {
        ### Compute the correlation on all assays
        cca.all<- lapply(
            levels(assay.names),
            function(x){
                if (!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]])) {
                    stop('To compute the regression,
                         the fast PCA must be computed first on the assay ', x, ' .')
                }
                ### Regression on PCs and continous variable
                printColoredMessage(message = paste0(
                    '### Computing regression based on PCs and the ',
                    variable,
                    ' variable.'
                ),
                color = 'magenta',
                verbose = verbose)
                pca_x <- se.obj@metadata[['metric']][[x]][['fastPCA']]$sing.val$u[colnames(se.obj),]
                cca.pcs<- sapply(
                    1:nb.pcs,
                    function(y){
                        ## Canonical correlations
                        cca <- cancor(
                            x = pca_x[, 1:y, drop = FALSE],
                            y = catvar.dummies)
                        1 - prod(1 - cca$cor^2)
                    })
                return(cca.pcs)
            })
        names(cca.all)<- levels(assay.names)
    } else{
        ### Compute the correlation on all assays
        cca.all<- lapply(
            levels(assay.names),
            function(x){
                if (!'PCA' %in% names(se.obj@metadata[['metric']][[x]])) {
                    stop('To compute the regression,
                         the PCA must be computed first on the assay ', x, ' .')
                }
                ### Regression on PCs and continous variable
                printColoredMessage(message = paste0(
                    '### Computing regression based on PCs and the ',
                    variable,
                    ' variable.'
                ),
                color = 'magenta',
                verbose = verbose)
                pca_x <- se.obj@metadata[['metric']][[x]][['PCA']]$sing.val$u[colnames(se.obj),]
                cca.pcs<- sapply(
                    1:nb.pcs,
                    function(y){
                        ## Canonical correlations
                        cca <- cancor(
                            x = pca_x[, 1:y, drop = FALSE],
                            y = catvar.dummies)
                        1 - prod(1 - cca$cor^2)
                    })
                return(cca.pcs)
            })
        names(cca.all)<- levels(assay.names)
    }

    ## Round the regression statistic obtained to 2 digits
    if(apply.round){
        cca.all[] <- lapply(cca.all, round,3)
    }

    ### Add results to the SummarizedExperiment object
    if(save.se.obj == TRUE){
        printColoredMessage(message= '### Saving the correlation results to the metadata of the SummarizedExperiment object.',
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
            if(!'pcs.vect.corr' %in% names(se.obj@metadata[['metric']][[x]])  ) {
                se.obj@metadata[['metric']][[x]][['pcs.vect.corr']] <- list()
            }
            ## Check if metadata metric already exist for this assay, this metric and this variable
            if(!variable %in% names(se.obj@metadata[['metric']][[x]][['pcs.vect.corr']])  ) {
                se.obj@metadata[['metric']][[x]][['pcs.vect.corr']][[variable]] <- cca.all[[x]]
            }
        }
        printColoredMessage(message= paste0(
            'The correlation results are saved to metadata@metric$',
            x,
            '$pcs.vect.corr$',
            variable, '.'),
            color = 'blue',
            verbose = verbose)

        ## Plot and save the plot into se.obj@metadata$plot
        if (plot.output==TRUE) {
            printColoredMessage(message= '### Plotting and Saving the correlation results to the metadata of the SummarizedExperiment object.',
                                color = 'magenta',
                                verbose = verbose)

            se.obj=plotMetric(se.obj,
                              assay.names =assay.names,
                              metric='pcs.vect.corr',
                              variable=variable,
                              verbose=verbose)
        }
        return(se.obj = se.obj)

        ## Return only the correlation result
    } else if(save.se.obj == FALSE){
        return(cca.all=cca.all)
    }

    printColoredMessage(message = '------------The PCVariableCorrelation function finished.',
                        color = 'white',
                        verbose = verbose)

}
