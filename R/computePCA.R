#' is used to compute PCA of the gene expression (assay) of a SummarizedExperiment class object.
#'
#' @param se A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the correlation. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#' @param fast.pca TO BE DEFINED.
#' @param nb.pcs TO BE DEFINED.
#' @param return.pc.percentage TO BE DEFINED.
#' @param scale Either a logical value or a numeric-alike vector of length equal
#' to the number of columns of the gene expression (assay) of a SummarizedExperiment class object.
#' It is a generic function to scale the columns of a numeric matrix, by default it is set to 'FALSE'.
#' @param center Either a logical value or a numeric-alike vector of length equal
#' to the number of columns of the gene expression (assay) of a SummarizedExperiment class object.
#' It is a generic function to center the columns of a numeric matrix, by default is set to 'TRUE'.
#' @param BSPARAM TO BE DEFINED.
#' @param save.se.obj Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#' @param pseudo.count TO BE DEFINED.


#' @return SummarizedExperiment A SummarizedExperiment object containing the PCA.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocSingular bsparam
#' @importFrom Rfast transpose
#' @importFrom BiocSingular runSVD
#' @import ggplot2
#' @export


computePCA <- function(se.obj,
                       assay.names = 'All',
                       apply.log = TRUE,
                       fast.pca = TRUE,
                       nb.pcs = 5,
                       return.pc.percentage = TRUE,
                       scale = FALSE,
                       center = TRUE,
                       BSPARAM,
                       save.se.obj = TRUE,
                       verbose = TRUE,
                       pseudo.count = 1) {
    printColoredMessage(message = '------------The computePCA function starts:',
                          color = 'white',
                          verbose = verbose)
    ### check the inputs
    if (fast.pca & is.null(nb.pcs)) {
        stop('To perform fast PCA, the number of PCs (nb.pcs) must specified.')
    } else if (fast.pca & nb.pcs == 0) {
        stop('To perform fast PCA, the number of PCs (nb.pcs) must specified.')
    } else if (is.null(assay.names)) {
        stop('Please provide at least an assay name.')
    }

    ### PCA
    if (assay.names[1] == 'All') {
        printColoredMessage(message = 'The PC of all individual assays of the SummarizedExperiment will be computed.',
                            color = 'blue',
                            verbose = verbose)
        assay.names <- names(assays(se.obj))
    }
    if (scale) {
        printColoredMessage(message = 'We highly recommend not to scale the data before computing the PCA.',
                            color = 'blue',
                            verbose = verbose)
    }
    if (fast.pca & return.pc.percentage) {
        printColoredMessage(message = 'Please note: if you selected the fast PCA analysis, the percentage of variation of PCs will be
                            computed proportional to the highest selected number of PCs, not on all the PCs.',
                            color = 'red',
                            verbose = verbose)
    }

    ## Assays
    if(length(assay.names) == 1 && assay.names=='All'){
        assay.names=as.factor(names(assays(se.obj)))
    }else {
        assay.names=as.factor(unlist(assay.names))
    }


    if (fast.pca) {
        ### svd
        printColoredMessage(
            message = paste0(
                '### Performing singular value decomposition with scale=',
                scale,
                ' and center= ',
                center,
                '.'
            ),
            color = 'magenta',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0('Obtaning the first ', nb.pcs, ' PCs.'),
            color = 'blue',
            verbose = verbose
        )
        sv.dec.all <- lapply(
            levels(assay.names),
            function(x) {
                                 temp.data = as.matrix(assay(x = se.obj,i= x))
                                 if (apply.log == FALSE)
                                     temp.data <- temp.data
                                 else
                                     temp.data <- log2(temp.data + pseudo.count)
                                 sv.dec <- runSVD(
                                     x = transpose(temp.data),
                                     k = nb.pcs,
                                     BSPARAM = BSPARAM,
                                     center = center,
                                     scale = scale
                                 )
                                 rownames(sv.dec$u) <- colnames(se.obj)
                                 rownames(sv.dec$v) <- row.names(se.obj)
                                 if(return.pc.percentage){
                                     percentage  <- sv.dec$d ^ 2 / sum(sv.dec$d ^ 2) * 100
                                     percentage  <- sapply(seq_along(percentage), function(i)
                                             round(percentage [i], 1))
                                     pca = list(sing.val = sv.dec, variation = percentage)
                                     return(pca)
                                 } else{
                                     pca = list(sing.val = sv.dec)
                                 }

        })
        names(sv.dec.all) <- levels(assay.names)
    } else {
        printColoredMessage(
            message = paste0(
                '### Performing singular value decomposition with scale=',
                scale,
                ' and center= ',
                center,
                '.'
            ),
            color = 'magenta',
            verbose = verbose
        )
        sv.dec.all <- lapply(
            levels(assay.names),
                             function(x) {
                                 temp.data = as.matrix(assay(x = se.obj,i= x))
                                 if (apply.log == FALSE)
                                     temp.data <- temp.data
                                 else
                                     temp.data <-
                                         log2(temp.data + pseudo.count)
                                 sv.dec <-
                                     svd(scale(
                                         x = transpose(temp.data),
                                         center = center,
                                         scale = scale
                                     ))

                                 if(return.pc.percentage){
                                     percentage  <-
                                         sv.dec$d ^ 2 / sum(sv.dec$d ^ 2) * 100
                                     percentage  <-
                                         sapply(seq_along(percentage), function(i)
                                             round(percentage [i], 1))
                                     pca = list(sing.val = sv.dec, variation = percentage)
                                     return(pca)
                                 } else{
                                     pca = list(sing.val = sv.dec)
                                 }
                                 return(pca)
                             })
        names(sv.dec.all) <- levels(assay.names)
    }

    ### Add results to the SummarizedExperiment object
    if(save.se.obj == TRUE){
        printColoredMessage(message =
                                '### Saving the PCA results to the metadata of the SummarizedExperiment object.',
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


            if (fast.pca) {
                    ## Check if metadata metric already exist for this assay and this metric
                    if(!'fastPCA' %in% names(se.obj@metadata[['metric']][[x]])  ) {
                        se.obj@metadata[['metric']][[x]][['fastPCA']] <- sv.dec.all[[x]]
                    }

            } else {
                    ## Check if metadata metric already exist for this assay and this metric
                    if(!'PCA' %in% names(se.obj@metadata[['metric']][[x]])  ) {
                        se.obj@metadata[['metric']][[x]][['PCA']] <- sv.dec.all[[x]]
                    }
            }
                printColoredMessage(message= paste0(
                    'The PCA results are saved to metadata@',
                    x,
                    '$PCA$.'),
                    color = 'blue',
                    verbose = verbose)

        }
    ## Return only the correlation result
    }else if(save.se.obj == FALSE){
    return(PCA=sv.dec.all[[x]])
}
        printColoredMessage(message = '------------The computePCA function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
}
