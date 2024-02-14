#' Create all possible homogeneous groups with respect to biological variables.

#' @author Ramyar Molania

#' @description
#' This function generates all possible homogeneous sample groups based on the specified biological variables. If continuous
#' variables are provided, the function splits them into a number of groups determined by 'nb.clusters', using the ckustering
#' method specified in 'clustering.method'. In the end, the product of all groups is generated and treated as homogeneous sample
#' groups with respect to biological variables.

#' @param se.obj A SummarizedExperiment object.
#' @param bio.variables Symbol. A symbol or a vector of symbols specifying the column names of biological variables in
#' the sample annotation of the SummarizedExperiment object. These 'bio.variables' can be either categorical or continuous
#' variables.
#' @param clustering.method Symbol. A symbol specifying the clustering method to be applied for grouping each continuous
#' source of biological variables. Options include 'kmeans', 'cut', and 'quantile'. The default is 'kmeans' clustering.
#' @param nb.clusters Numeric. A value indicating the number of groups for continuous sources of biological variation.
#' The default is 3. This implies that each continuous source will be split into 3 groups using the specified
#' 'clustering.method' method.
#' @param assess.variables Logical. Indicates whether to assess correlations between the biological variables. If 'TRUE',
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
#' or a vector of all possible homogeneous groups.

#' @importFrom SummarizedExperiment assay
#' @importFrom DescTools ContCoef
#' @importFrom stats kmeans quantile
#' @importFrom knitr kable

createHomogeneousBioGroups <- function(
        se.obj,
        bio.variables,
        clustering.method = 'kmeans',
        nb.clusters,
        assess.variables = TRUE,
        cat.cor.coef = c(0.9, 0.9),
        cont.cor.coef = c(0.9, 0.9),
        assess.se.obj = TRUE,
        remove.na = 'sample.annotation',
        save.se.obj = TRUE,
        verbose = TRUE) {
    printColoredMessage(message = '------------The createHomogeneousBioGroups function starts:',
                        color = 'white',
                        verbose = verbose)
    # check some inputs ###
    if (is.null(bio.variables)) {
        stop('The bio.variables cannot be empty.')
    } else if (!clustering.method %in% c('kmeans', 'cut', 'quantile')) {
        stop('The "clustering.method" should be one of: "kmeans", "cut" or "quantile".')
    } else if (max(cat.cor.coef) > 1 | min(cat.cor.coef) < 0) {
        stop('The maximum value for the "cat.cor.coef" cannot be more than 1 or negative.')
    } else if (max(cont.cor.coef) > 1 | min(cont.cor.coef) < 0) {
        stop('The maximum value for the "cont.cor.coef" cannot be more than 1 or negative.')
    } else if (remove.na %in% c('both', 'measuerments')) {
        stop('The remove.na should be either "sample.annotation" or "none".')
    }
    # assess correlation between the variables ####
    if (isTRUE(assess.variables)) {
        se.obj <- assessVariablesAssociation(
            se.obj = se.obj,
            bio.variables = bio.variables,
            uv.variables = NULL,
            cat.cor.coef = cat.cor.coef,
            cont.cor.coef = cont.cor.coef,
            assess.se.obj = TRUE,
            remove.na = remove.na,
            verbose = verbose
        )
        bio.variables <- se.obj$bio.variables
        se.obj <- se.obj$se.obj
    }
    class.bio.var <- unlist(lapply(
        bio.variables,
        function(x) class(se.obj[[x]])))
    categorical.bio.var <- bio.variables[class.bio.var %in% c('factor', 'character')]
    continuous.bio.var <- bio.variables[class.bio.var %in% c('numeric', 'integer')]
    # cluster continuous variables ####
    if (length(continuous.bio.var) > 0) {
        # grammar
        if (length(continuous.bio.var) == 1) {
            a = 'is'
            b = 'source'
        } else {
            a = 'are'
            b = 'sources'
        }
        printColoredMessage(
            message = paste0('-- Clustering of the ', b, ' of biological variation:'),
            color = 'magenta',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0(
                'There ',
                a,
                ' ',
                length(continuous.bio.var),
                ' continuous ' ,
                b ,
                ' of biological variation:',
                paste0(continuous.bio.var, collapse = ' & '),
                '.'
            ),
            color = 'blue',
            verbose = verbose
        )
        if (nb.clusters == 1) {
            stop('The value of "nb.clusters" should be more than 1.')
        } else if (is.null(nb.clusters)) {
            stop('The "nb.clusters" cannot be empty.')
        } else if (nb.clusters == 0) {
            stop('The "nb.clusters" cannot be 0.')
        }
        printColoredMessage(
            message = paste0(
                'Then, each source will be divided into ',
                nb.clusters,
                ' groups using ',
                clustering.method,
                '.'
            ),
            color = 'blue',
            verbose = verbose
        )
        if (clustering.method == 'kmeans') {
            set.seed(3456)
            bio.cont.groups <- lapply(
                continuous.bio.var,
                function(x) {
                    bio.cont.clusters <- kmeans(
                        x = colData(se.obj)[[x]],
                        centers = nb.clusters,
                        iter.max = 1000)
                    paste0(x, '_group', bio.cont.clusters$cluster)
                })
            names(bio.cont.groups) <- continuous.bio.var
        } else if (clustering.method == 'cut') {
            bio.cont.groups <- lapply(
                continuous.bio.var,
                function(x) {
                    bio.cont.clusters <- as.numeric(cut(
                        x = colData(se.obj)[[x]],
                        breaks = nb.clusters,
                        include.lowest = TRUE
                    ))
                    paste0(x, '_group', bio.cont.clusters)
                })
            names(bio.cont.groups) <- continuous.bio.var
        } else if (clustering.method == 'quantile') {
            bio.cont.groups <- lapply(
                continuous.bio.var,
                function(x) {
                    quantiles <- quantile(x = colData(se.obj)[[x]],
                                 probs = seq(0, 1, 1 / nb.clusters))
                    bio.cont.clusters <- as.numeric(cut(
                            x = colData(se.obj)[[x]],
                            breaks = quantiles,
                            include.lowest = TRUE
                        ))
                    paste0(x, '_group', bio.cont.clusters)
                })
            names(bio.cont.groups) <- continuous.bio.var
        }
        continuous.bio.groups <- as.data.frame(do.call(cbind, bio.cont.groups))
    }
    # categorical biological variables ####
    if (length(categorical.bio.var) > 0) {
        printColoredMessage(
            message = '-- Check the categorical sources:',
            color = 'magenta',
            verbose =  verbose
        )
        if (length(categorical.bio.var) == 1) {
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
                length(categorical.bio.var),
                ' categorical ' ,
                b ,
                ' of biological variation.'
            ),
            color = 'blue',
            verbose = verbose
        )
        categorical.bio.groups <- colData(se.obj)[, categorical.bio.var, drop = FALSE]
        categorical.bio.groups <- as.data.frame(categorical.bio.groups)
    }
    # create all possible biological groups ####
    if (length(categorical.bio.var) > 0 &
        length(continuous.bio.var) > 0) {
        printColoredMessage(
            message = '-- The products of the groups:',
            color = 'magenta',
            verbose =  verbose
        )
        all.groups <-
            cbind(continuous.bio.groups, categorical.bio.groups)
        all.groups <- unname(apply(all.groups,
                                   1,
                                   paste , collapse = ".."))
        printColoredMessage(
            message = paste0(
                'In totall, ',
                length(unique(all.groups)),
                ' homogeneous biological groups are created:'),
            color = 'blue',
            verbose = verbose
        )
        if (verbose)
            print(kable(table(all.groups), caption = 'Homogeneous biological groups'))
        all.groups <- gsub('_', '-', all.groups)
    } else if (length(categorical.bio.var) == 0 &
               length(continuous.bio.var) > 0) {
        printColoredMessage(
            message = '-- The products of the groups:',
            color = 'magenta',
            verbose =  verbose
        )
        all.groups <- unname(apply(
            continuous.bio.groups,
            1,
            paste , collapse = ".."))
        printColoredMessage(
            message = paste0(
                'In totall, ',
                length(unique(all.groups)),
                ' homogeneous biological groups are created:'),
            color = 'grey',
            verbose = verbose
        )
        if (verbose)
            print(kable(table(all.groups), caption = 'Homogeneous biological groups'))
        all.groups <- gsub('_', '-', all.groups)
    } else if (length(categorical.bio.var) > 0 &
               length(continuous.bio.var) == 0) {
        printColoredMessage(
            message = '-- The products of the groups:',
            color = 'magenta',
            verbose =  verbose
        )
        all.groups <- unname(apply(categorical.bio.groups,
                                   1,
                                   paste , collapse = ".."))
        printColoredMessage(
            message = paste0(
                'In totall, ',
                length(unique(all.groups)),
                ' homogeneous biological groups are created:'),
            color = 'grey',
            verbose = verbose
        )
        if (verbose)
            print(kable(table(all.groups), caption = 'Homogeneous biological groups'))
        all.groups <- gsub('_', '-', all.groups)
    }
    # add results to the SummarizedExperiment object ####
    out.put.name <- paste0(
        'HBIOG:',
        length(unique(all.groups)),
        ' Gro||Var_',
        paste0(bio.variables, collapse = '&'),
        '||Clust:',
        clustering.method,
        '_nb.clust:',
        nb.clusters,
        '.'
    )
    printColoredMessage(message = '-- Save the results',
                        color = 'magenta',
                        verbose = verbose)
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = '-- Saving the homogeneous groups to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        ## Check if metadata NCG already exists
        if (length(se.obj@metadata$HGgroups) == 0) {
            se.obj@metadata[['HGgroups']] <- list()
        }
        se.obj@metadata[['HGgroups']][['BioGroups']][[out.put.name]] <- all.groups
        printColoredMessage(message = 'The homogeneous groups are saved to the metadata of the SummarizedExperiment object.',
                            color = 'blue',
                            verbose = verbose)
        printColoredMessage(message = '------------The createHomogeneousBioGroups function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else{
        printColoredMessage(
            message = '-- The homogeneous groups are outputed as a vector.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The createHomogeneousBioGroups function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.groups)
    }
}
