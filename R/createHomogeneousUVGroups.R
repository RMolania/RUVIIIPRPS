#' Create all possible homogeneous groups with respect to unwanted variables.

#' @author Ramyar Molania

#' @description
#' This function generates all possible homogeneous sample groups based on the specified unwanted variables. If continuous
#' variables are provided, the function splits them into a number of groups determined by 'nb.clusters', using the method
#' specified in clustering.method'. In the end, the product of all groups is generated and treated as homogeneous sample
#' groups with respect to unwanted variables.

#' @param se.obj A SummarizedExperiment object.
#' @param uv.variables Symbol. A symbol or a vector of symbols specifying the column names of unwanted variables in
#' the sample annotation of the SummarizedExperiment object. These 'uv.variables' can be either categorical or continuous
#' variables.
#' @param clustering.method Symbol. A symbol specifying the clustering method to be applied for grouping each continuous
#' source of unwanted variables. Options include 'kmeans', 'cut', and 'quantile'. The default is 'kmeans' clustering.
#' @param nb.clusters Numeric. A value indicating the number of groups for continuous sources of unwanted variation.
#' The default is 3. This implies that each continuous source will be split into 3 groups using the specified
#' 'clustering.method' method.
#' @param assess.variables Logical. Indicates whether to assess correlations between the unwanted variables. If 'TRUE',
#' the function 'assessVariablesCorrelation' will be applied. For more details refer to the 'assessVariablesCorrelation'
#' function.The default is 'TRUE'.
#' @param cat.cor.coef Vector. A vector of two numerical values. Indicates the cut-off of the correlation coefficient
#' between each pair of categorical variables. The first one is between each pair of 'uv.variables' and the second one
#' is between each pair of 'bio.variables'. The correlation is computed by the function ContCoef from the DescTools
#' package. If the correlation of a pair of variable is higher than the cut-off, then only the variable that has the
#' highest number of factor will be kept and the other one will be excluded from the remaining analysis. By default they
#' are both set to 0.9.
#' @param cont.cor.coef Vector of two numerical values. Indicates the cut-off of the Spearman correlation coefficient
#' between each pair of continuous variables. The first one is between each pair of 'uv.variables' and the second one is
#' between each pair of 'bio.variables'. If the correlation of a pair of variable is higher than the cut-off, then only
#' the variable that has the highest variance will be kept and the other one will be excluded from the remaining analysis.
#' By default they are both set to 0.9.
#' @param assess.se.obj Logical. Whether to assess the SummarizedExperiment object or not. If 'TRUE', the function
#' 'checkSeobj' will be applied. The default is 'TRUE'.
#' @param remove.na Symbol. Indicates whether to remove missing values from the 'uv.variables'. The options are
#' 'sample.annotation' or 'none'. The default is 'sample.annotation', indicating the missing values from the variables
#' will be removed.
#' @param save.se.obj Logical. Indicates whether to save the results to the metadata of the SummarizedExperiment object
#' or not. If 'TRUE', all the possible homogeneous groups will be saved into "metadata$HGgroups$UVgroups", otherwise
#' the results will outputted as a vector. The default is 'TRUE'.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object containing the all possible homogeneous groups in the "metadata$HGgroups$UVgroups"
#' or a vector.

#' @importFrom SummarizedExperiment assay
#' @importFrom DescTools ContCoef
#' @importFrom stats kmeans quantile
#' @importFrom knitr kable

createHomogeneousUVGroups <- function(
        se.obj,
        uv.variables,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        assess.variables = TRUE,
        cat.cor.coef = c(0.9, 0.9),
        cont.cor.coef = c(0.9, 0.9),
        assess.se.obj = TRUE,
        remove.na = 'sample.annotation',
        save.se.obj = TRUE,
        verbose = TRUE) {
    printColoredMessage(message = '------------The createHomogeneousUVGroups function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if(is.null(uv.variables)){
        stop('The uv.variables cannot be empty.')
    } else if(!is.vector(uv.variables)){
        stop('The "uv.variables" must be a vector.')
    } else if(!clustering.method %in% c('kmeans', 'cut', 'quantile')){
        stop('The clustering.method should be one of: kmeans, cut or quantile.')
    } else if(max(cat.cor.coef) > 1){
        stop('The maximum value for the cat.cor.coef is 1.')
    } else if(max(cont.cor.coef) > 1){
        stop('The maximum value for the cont.cor.coef is 1.')
    } else if( remove.na %in% c('both', 'measuerments')){
        stop('The remove.na should be either sample.annotation or none.')
    }
    # assess variables correlation ####
    printColoredMessage(
        message = '-- Assessing the correlation between unwanted variation variables:',
        color = 'magenta',
        verbose = verbose)
    if (assess.variables) {
        se.obj <- assessVariablesCorrelation(
            se.obj = se.obj,
            bio.variables = NULL,
            uv.variables = uv.variables,
            cat.cor.coef = cat.cor.coef,
            cont.cor.coef = cont.cor.coef,
            assess.se.obj = TRUE,
            remove.na = remove.na,
            verbose = verbose)
        uv.variables <- se.obj$uv.variables
        se.obj <- se.obj$se.obj
    }
    # class of  uv variables ####
    class.uv.var <- unlist(lapply(
        uv.variables,
        function(x) class(se.obj[[x]])))
    categorical.uv <- uv.variables[class.uv.var %in% c('factor', 'character')]
    continuous.uv <- uv.variables[class.uv.var %in% c('numeric', 'integer')]
    ## check continuous variable ####
    if(length(continuous.uv) > 0){
        ## check clustering value for continuous variable
        if(nb.clusters == 1){
            stop('The value of "nb.clusters" must be bigger than 1.')
        } else if (is.null(nb.clusters)){
            stop('The "nb.clusters" cannot be empty.')
        } else if ( nb.clusters == 0) {
            stop('The value of "nb.clusters" cannot be 0.')
        }
        ### grammar
        if(length(continuous.uv) == 1){
            a = 'is'
            b = 'source'
        } else {
            a = 'are'
            b = 'sources'
        }
        ### clustering ####
        printColoredMessage(
            message = paste0('-- Clustering of the ', b, ' of unwanted variation:'),
            color = 'magenta',
            verbose = verbose)
            printColoredMessage(
                message = paste0(
                    'There ',
                    a,
                    ' ',
                    length(continuous.uv),
                    ' continuous ' ,
                    b ,
                    ' of unwanted variation:',
                    paste0(continuous.uv, collapse = ' & '),
                    '.'),
                color = 'blue',
                verbose = verbose)
            printColoredMessage(message = paste0(
                'Then, each source will be divided into ', nb.clusters, ' groups using ', clustering.method, '.'),
                color = 'blue',
                verbose = verbose)
            if (clustering.method == 'kmeans') {
                    set.seed(3456)
                    continuous.uv.groups <- lapply(
                        continuous.uv,
                        function(x){
                            uv.cont.clusters <- stats::kmeans(
                                x = colData(se.obj)[[x]],
                                centers = nb.clusters,
                                iter.max = 1000)
                            paste0(x, '_group', uv.cont.clusters$cluster)
                        })
                    names(continuous.uv.groups) <- continuous.uv
                    continuous.uv.groups <- as.data.frame(do.call(cbind, continuous.uv.groups))
                    continuous.uv.groups <- apply(
                        continuous.uv.groups,
                        1,
                        paste , collapse = "..")
            } else if (clustering.method == 'cut') {
                continuous.uv.groups <- lapply(
                    continuous.uv,
                    function(x) {
                        uv.cont.clusters <- as.numeric(
                            cut(x = colData(se.obj)[[x]],
                                breaks = nb.clusters,
                                include.lowest = TRUE))
                        paste0(x, '_group', uv.cont.clusters)
                    })
                names(continuous.uv.groups) <- continuous.uv
                continuous.uv.groups <- as.data.frame(do.call(cbind, continuous.uv.groups))
                continuous.uv.groups <- apply(
                    continuous.uv.groups,
                    1,
                    paste , collapse = "..")
            } else if (clustering.method == 'quantile') {
                continuous.uv.groups <- lapply(
                    continuous.uv,
                    function(x) {
                        quantiles <- quantile(x = colData(se.obj)[[x]], probs = seq(0, 1, 1 / nb.clusters))
                        uv.cont.clusters <- as.numeric(
                            cut(x = colData(se.obj)[[x]],
                                breaks = quantiles,
                                include.lowest = TRUE))
                        paste0(x, '_group', uv.cont.clusters)
                    })
                names(continuous.uv.groups) <- continuous.uv
                continuous.uv.groups <- as.data.frame(do.call(cbind, continuous.uv.groups))
                continuous.uv.groups <- apply(
                    continuous.uv.groups,
                    1,
                    paste , collapse = "..")
            }
        }
    ## check categorical variable ####
    if(length(categorical.uv) > 0){
        printColoredMessage(
            message = '-- Check the categorical sources:',
            color = 'magenta',
            verbose =  verbose
        )
        if (length(categorical.uv) == 1) {
            a = 'is'
            b = 'source'
        } else {
            a = 'are'
            b = 'sources'
        }
        printColoredMessage(
            message = paste0(
                'There ',
                a,
                ' ',
                length(categorical.uv),
                ' categorical ' ,
                b ,
                ' of unwanted variation.'
            ),
            color = 'blue',
            verbose = verbose
        )
        categorical.uv.groups <- colData(se.obj)[, categorical.uv, drop = FALSE]
        categorical.uv.groups <- as.data.frame(categorical.uv.groups)
    }
    # create all possible groups ####
    if(length(categorical.uv) > 0 & length(continuous.uv) > 0 ){
        printColoredMessage(
            message = '-- The products of the groups:',
            color = 'magenta',
            verbose =  verbose
        )
        all.groups <- categorical.uv.groups
        all.groups$cont <- continuous.uv.groups
        all.groups <- unname(apply(
            all.groups,
            1,
            paste , collapse = ".."
        ))
        printColoredMessage(
            message = paste0(
                'In totall, ',
                length(unique(all.groups)),
                ' homogeneous groups with respect to sources of unwanted variation are created:'),
            color = 'blue',
            verbose = verbose
        )
        if(verbose) print(kable(table(all.groups), caption = 'homogeneous groups'))
        all.groups <- gsub('_', '-', all.groups)
    } else if (length(categorical.uv) == 0 & length(continuous.uv) > 0){
        printColoredMessage(
            message = '-- The products of the groups:',
            color = 'magenta',
            verbose =  verbose
        )
        all.groups <- continuous.uv.groups
        printColoredMessage(
            message = paste0(
                'In totall, ',
                length(unique(all.groups)),
                ' homogeneous groups with respect to sources of unwanted variation are created:'),
            color = 'blue',
            verbose = verbose)
        if(verbose) print(kable(table(all.groups), caption = 'homogeneous groups'))
        all.groups <- gsub('_', '-', all.groups)
    } else if(length(categorical.uv) > 0 & length(continuous.uv) == 0){
        printColoredMessage(
            message = '-- The products of the groups:',
            color = 'magenta',
            verbose =  verbose
        )
        all.groups <- apply(
            categorical.uv.groups,
            1,
            paste , collapse = "..")
        printColoredMessage(
            message = paste0(
                'In totall, ',
                length(unique(all.groups)),
                ' homogeneous groups with respect to sources of unwanted variation are created:'),
            color = 'blue',
            verbose = verbose)
        if(verbose) print(kable(table(all.groups), caption = 'homogeneous groups'))
        all.groups <- gsub('_', '-', all.groups)
    }
    # add results to the SummarizedExperiment object ####
    out.put.name <- paste0(
        length(unique(all.groups)),
        ' groups||UV',
        paste0(uv.variables, collapse = '&'),
        '||Clus:',
        clustering.method,
        '_nb.clus:',
        nb.clusters,
        '.')
    printColoredMessage(message = '-- Save the results',
        color = 'magenta',
        verbose = verbose)
    if(save.se.obj == TRUE){
        printColoredMessage(
            message = '-- Save the homogeneous groups to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        ## Check if metadata NCG already exists
        if(length(se.obj@metadata$HGgroups) == 0 ) {
            se.obj@metadata[['HGgroups']] <- list()
        }
        se.obj@metadata[['HGgroups']][['UVgroups']][[out.put.name]] <- all.groups
        printColoredMessage(
            message = 'The homogeneous groups are saved to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The createHomogeneousUVGroups function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj)
    } else{
        printColoredMessage(
            message = '-- The homogeneous groups are outputed as a vector.',
            color = 'blue',
            verbose = verbose)
       printColoredMessage(
        message = '------------The createHomogeneousUVGroups function finished.',
        color = 'white',
        verbose = verbose)
    return(all.groups)
    }
}

