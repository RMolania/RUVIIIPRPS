#' is used to compute the linear regression between the the first cumulative PCs
#' of the gene expression (assay) of a SummarizedExperiment class object and a continuous variable (i.e. library size)
#'
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the regression. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable String of the label of a continuous variable such as
#' library size from colData(se.obj).
#' @param fast.pca TO BE DEFINED.
#' @param nb.pcs TO BE DEFINED.
#' @param save.se.obj Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param plot.output Indicates whether to plot the F-test statistics, by default it is set to TRUE.
#' @param assess.se.obj Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na TO BE DEFINED.
#' @param apply.round TO BE DEFINED.
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the computed regression for
#' the continuous variable and if requested the associated plot.
#' @importFrom stats lm
#' @importFrom wesanderson wes_palette
#' @importFrom dplyr rename mutate
#' @importFrom tidyr pivot_longer %>%
#' @import ggplot2
#' @export

PCVariableRegression<-function(
        se.obj,
        assay.names = 'All',
        variable,
        fast.pca = FALSE,
        nb.pcs = 10,
        save.se.obj = TRUE,
        plot.output=TRUE,
        assess.se.obj = TRUE,
        remove.na = 'both',
        apply.round = TRUE,
        verbose = TRUE
){

    printColoredMessage(message = '------------The PCVariableRegression function starts:',
                        color = 'white',
                        verbose = verbose)

    ### check the inputs
    if (is.null(nb.pcs)) {
        stop('To compute the regression, the number of PCs (nb.pcs) must be specified.')
    } else if (is.null(assay.names)) {
        stop('Please provide at least an assay name.')
    }

    ### Assess the se.obj
    if (assess.se.obj){
        se.obj <- checkSeObj(se.obj = se.obj,
                             assay.names = assay.names,
                             variables = variable,
                             remove.na = remove.na,
                             verbose = verbose)
    }


    ### Check if the variable provided has a variance of 0:
    if(var(se.obj[[variable]]) == 0){
        stop(paste0('The ', variable, ' appears to have no variation.'))
    }
    var=se.obj@colData[, variable]

    ## Assays
    if(length(assay.names) == 1 && assay.names=='All'){
        assay.names=as.factor(names(assays(se.obj)))
    }else {
        assay.names=as.factor(unlist(assay.names))
    }

    if (fast.pca) {
        ### Compute the regression on all assays
        lreg.pcs<- lapply(
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
                pca_x <- se.obj@metadata[['metric']][[x]][['fastPCA']]$sing.val$u
                rSquared <- sapply(
                    1:nb.pcs,
                    function(y) {
                        lm.ls <- summary(lm(
                            var ~ pca_x[, 1:y])
                        )$r.squared
                    })
                return(rSquared)
            })
        names(lreg.pcs) <- levels(assay.names)
    } else{
        ### Compute the regression on all assays
        lreg.pcs<- lapply(
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
                pca_x <- se.obj@metadata[['metric']][[x]][['PCA']]$sing.val$u
                rSquared <- sapply(
                    1:nb.pcs,
                    function(y) {
                        lm.ls <- summary(lm(
                            cont_var ~ pca_x[, 1:y])
                        )$r.squared
                    })
                return(rSquared)
            })
        names(lreg.pcs) <- levels(assay.names)
    }

    ## Round the regression statistic obtained to 2 digits
    if(apply.round){
        lreg.pcs[] <- lapply(lreg.pcs, round,3)
    }


    ### Add results to the SummarizedExperiment object
    if(save.se.obj == TRUE){
        printColoredMessage(message= '### Saving the regression results to the metadata of the SummarizedExperiment object.',
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
            if(!'pcs.lm' %in% names(se.obj@metadata[['metric']][[x]])  ) {
                se.obj@metadata[['metric']][[x]][['pcs.lm']] <- list()
            }
            ## Check if metadata metric already exist for this assay, this metric and this variable
            if(!variable %in% names(se.obj@metadata[['metric']][[x]][['pcs.lm']])  ) {
                se.obj@metadata[['metric']][[x]][['pcs.lm']][[variable]] <- lreg.pcs[[x]]
            }
        }
        printColoredMessage(message= paste0(
            'The anova results are saved to metadata@',
            x,
            '$pcs.lm$',
            variable, '.'),
            color = 'blue',
            verbose = verbose)

        ## Plot and save the plot into se.obj@metadata$plot
        if (plot.output==TRUE) {
            printColoredMessage(message= '### Plotting and Saving the regression results to the metadata of the SummarizedExperiment object.',
                                color = 'magenta',
                                verbose = verbose)

            se.obj=plotMetric(se.obj,
                              assay.names =assay.names,
                              metric='pcs.lm',
                              variable=variable,
                              verbose=verbose)
        }
        return(se.obj = se.obj)

        ## Return only the correlation result
    } else if(save.se.obj == FALSE){
        return(gene.anova.var=lreg.pcs)
    }

    printColoredMessage(message = '------------The PCVariableRegression function finished.',
                        color = 'white',
                        verbose = verbose)

}
