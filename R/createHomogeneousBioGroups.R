#' This is used to create all possible homogeneous biological groups in the SummarizedExperiment objct.
#'
#'
#' @param se.obj A SummarizedExperiment object.
#' @param bio.variables Symbol. Indicates the column names biological variables the SummarizedExperiment object.
#' @param nb.clusters A value to specify the number of groups of continuous sources of biological variation..
#' @param clustering.method A clustering method to group each continuous sources of biological variation.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object.
#' @param assess.variables Logical. Indicates whether to assess association between the biological variables.
#' @param cat.cor.coef Numeric. Indicates correlation coefficients (Cramer's V, Pearson's contingency coefficient)
#' cut off for assessing association between categorical sources of biological variation.
#' @param cont.cor.coef Numeric. Indicates Spearman coefficients cut off for assessing association between continuous sources of biological variation.
#' @param remove.na Logical. To remove NA or missing values from sample annotation or not.
#' @param save.se.obj Logical. Indicates whether to save the results to SummarizedExperiment object or not.
#' @param verbose Logical. Whether to show the messages of the functions or not.

#' @importFrom SummarizedExperiment assay
#' @importFrom DescTools ContCoef
#' @importFrom stats kmeans quantile
#' @importFrom knitr kable

createHomogeneousBioGroups <- function(
        se.obj,
        bio.variables,
        nb.clusters,
        clustering.method = 'kmeans',
        assess.se.obj = TRUE,
        assess.variables = TRUE,
        cat.cor.coef = c(0.7, 0.7),
        cont.cor.coef = c(0.7, 0.7),
        remove.na = 'sample.annotation',
        save.se.obj = TRUE,
        verbose = TRUE) {
    printColoredMessage(message = '------------The createHomogeneousBioGroups function starts:',
                        color = 'white',
                        verbose = verbose)
    print('ffff')
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
    if (assess.variables) {
        se.obj <- variablesCorrelation(
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
    class.bio.var <- unlist(lapply(bio.variables,
                                   function(x)
                                       class(se.obj[[x]])))
    categorical.bio.var <-
        bio.variables[class.bio.var %in% c('factor', 'character')]
    continuous.bio.var <-
        bio.variables[class.bio.var %in% c('numeric', 'integer')]
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
            message = paste0('-- Clustering of the ', b, ' of unwanted variation:'),
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
            bio.cont.groups <- lapply(continuous.bio.var,
                                      function(x) {
                                          bio.cont.clusters <- kmeans(
                                              x = colData(se.obj)[[x]],
                                              centers = nb.clusters,
                                              iter.max = 1000
                                          )
                                          paste0(x, '_group', bio.cont.clusters$cluster)
                                      })
            names(bio.cont.groups) <- continuous.bio.var
        } else if (clustering.method == 'cut') {
            bio.cont.groups <- lapply(continuous.bio.var,
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
            bio.cont.groups <- lapply(continuous.bio.var,
                                      function(x) {
                                          quantiles <-
                                              quantile(x = colData(se.obj)[[x]],
                                                       probs = seq(0, 1, 1 / nb.clusters))
                                          bio.cont.clusters <-
                                              as.numeric(cut(
                                                  x = colData(se.obj)[[x]],
                                                  breaks = quantiles,
                                                  include.lowest = TRUE
                                              ))
                                          paste0(x, '_group', bio.cont.clusters)
                                      })
            names(bio.cont.groups) <- continuous.bio.var
        }
        continuous.bio.groups <-
            as.data.frame(do.call(cbind, bio.cont.groups))
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
        categorical.bio.groups <-
            colData(se.obj)[, categorical.bio.var, drop = FALSE]
        categorical.bio.groups <-
            as.data.frame(categorical.bio.groups)
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
                ' homogeneous biological groups are created:'
            ),
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
        all.groups <- unname(apply(continuous.bio.groups,
                                   1,
                                   paste , collapse = ".."))
        printColoredMessage(
            message = paste0(
                'In totall, ',
                length(unique(all.groups)),
                ' homogeneous biological groups are created:'
            ),
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
                ' homogeneous biological groups are created:'
            ),
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
    if (save.se.obj == TRUE) {
        printColoredMessage(
            message = '-- Saving the homogeneous groups to the metadata of the SummarizedExperiment object.',
            color = 'magenta',
            verbose = verbose
        )
        ## Check if metadata NCG already exists
        if (length(se.obj@metadata$HGgroups) == 0) {
            se.obj@metadata[['HGgroups']] <- list()
        }
        se.obj@metadata[['HGgroups']][['BioGroups']][[out.put.name]] <-
            all.groups
        printColoredMessage(message = 'The homogeneous groups are saved to the metadata of the SummarizedExperiment object.',
                            color = 'blue',
                            verbose = verbose)
        printColoredMessage(message = '------------The createHomogeneousBioGroups function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else{
        printColoredMessage(message = '------------The createHomogeneousBioGroups function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.groups)
    }
}
