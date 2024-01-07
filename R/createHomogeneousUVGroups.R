#' This is used to create all possible groups with respect to unwanted variation variables in the SummarizedExperiment objct.
#'
#' @param se.obj A summarized experiment object.
#' @param uv.variables Symbol. Indicate of the columns names of unwanted variation variables to specify major groups
#' with respect to unwanted variation.
#' @param nb.clusters Numeric. A value to specify the number of groups for continuous sources of biological variation.
#' The default is 2. This means each continuous sources will be divided inot 2 groups using the 'clustering.method'.
#' @param clustering.method A clustering method to group each continuous sources of biological variation.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param assess.variables Logical. Indicates whether to assess association between the unwanted variables. For more
#' details refer to the 'variableCorrelation' function.
#' @param cat.cor.coef Numeric. Correlation coefficients (Cramer's V, Pearson's contingency coefficient) cut off for
#' assessing association between categorical sources of biological variation.
#' @param cont.cor.coef Spearman coefficients cut off for assessing association between continuous sources of biological
#' variation.
#' @param remove.na To remove NA or missing values from either the assays or sample annotation or both.
#' @param save.se.obj Logical. Indicates whether to save the results to SummarizedExperiment object or not.
#' @param verbose Logical. Whether to show the messages of the functions or not.

#' @return SummarizedExperiment A SummarizedExperiment object containing the a set of sutable negatve control genes
#' seleced by the two-way ANOVA approach.

#' @importFrom SummarizedExperiment assay
#' @importFrom DescTools ContCoef
#' @importFrom stats kmeans quantile
#' @importFrom knitr kable

createHomogeneousUVGroups <- function(
        se.obj,
        uv.variables,
        nb.clusters,
        clustering.method = 'kmeans',
        assess.se.obj = TRUE,
        assess.variables = TRUE,
        cat.cor.coef = c(0.7, 0.7),
        cont.cor.coef = c(0.7, 0.7),
        remove.na = 'sample.annotation',
        save.se.obj = TRUE,
        verbose = TRUE) {
    printColoredMessage(message = '------------The createHomogeneousUVGroups function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if(is.null(uv.variables)){
        stop('The uv.variables cannot be empty.')
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
        se.obj <- variablesCorrelation(
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
            stop('The value of "nb.clusters" should be bigger than 1.')
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
                            uv.cont.clusters <- kmeans(
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
    # create all possible biological groups ####
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
        'HUVG:',
        length(unique(all.groups)),
        ' groups||UVVariables_',
        paste0(uv.variables, collapse = '&'),
        '||Clusetring:',
        clustering.method,
        '_nb.clusters:',
        nb.clusters,
        '.')
    if(save.se.obj == TRUE){
        printColoredMessage(
            message = '-- Saving the homogeneous groups to the metadata of the SummarizedExperiment object.',
            color = 'magenta',
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
        message = '------------The createHomogeneousUVGroups function finished.',
        color = 'white',
        verbose = verbose)
    return(all.groups)
        }
}

