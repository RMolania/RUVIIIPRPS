#' create supervised pseudo-replicates of pseudo samples (PRPS).

#' @author Ramyar Molania

#' @description
#' This function creates different PRPS sets when all sources of unwanted and biological variation are known. We will create
#' distinct group of pseudo-replicates for each source of unwanted variation defined in the 'uv.variables' argument. For
#' example to correct for batch effect if defined in the 'uv.variables' argument, several group of pseudo-samples will be
#' created by averaging the samples of the same biological subtype defined in 'bio.variables' in each batch. Then those
#' pseudo-samples will be defined as pseudo-replicates which constitutes a PRPS set. For example to correct for library
#' size if defined in the 'uv.variables' argument, several group of pseudo-samples will be created by averaging the top
#' and bottom-ranked samples by library size of the same biological subtype in each batch. Then those pseudo-samples will
#' be defined as pseudo-replicates which constitutes a PRPS set. Similarly to correct for purity if defined in the
#' 'uv.variables' argument, several group of pseudo-samples will be created by averaging the top and bottom-ranked samples
#' by purity of the same biological subtype in each batch. Then those pseudo-samples will be defined as pseudo-replicates
#' which constitutes a PRPS set.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name String for the selection of the name of the assay of the SummarizedExperiment class object.
#' @param uv.variables String of the label of a categorical variable
#' such as samples types or batch from colData(se) that will be used to define PRPS.
#' @param bio.variables String of the label of (a) categorical or continuous variable(s) used to define homogeneous
#' biological groups of samples such as samples types from colData(se).
#' @param apply.other.uv.variables String of the label of (a) categorical or continuous variable(s) used to define
#' homogeneous groups of samples such as library size, plates from colData(se). By default it is set to 'NULL' meaning the
#' samples will be assigned to homogeneous biological groups of samples using the 'bio.variables'.
#' @param min.sample.for.prps Numeric. Indicates the minimum number of samples to create one pseudo-sample,
#' by default it is set to 3.
#' to create a PRPS set for each continuous variable. The minimum should be '2*min.sample.for.prps'. By default it is set
#' to 6.
#' @param bio.clustering.method String of the clustering method to assign each sample to an homogeneous biological
#' clusters/group using the function 'kmeans', cut' from base or the function 'quantile'. By default it is to 'kmeans'.
#' @param other.uv.clustering.method String of the clustering method to assign each sample to an homogeneous clusters/group
#' of samples based on the 'other.uv.variables' using the function 'kmeans','cut' from base or the function 'quantile'.
#' By default it is to 'kmeans'.
#' @param nb.bio.clusters Numeric. A value to specify the number of homogeneous biological clusters/groups of samples.
#' By default it is set to 3.
#' @param nb.other.uv.clusters Numeric. A value to specify the number of groups for continuous sources of biological
#' variation. The default is 2. This means each continuous sources will be divided into 2 groups using the 'clustering.method'.
#' @param cat.cor.coef Vector of two numerical values. Indicates the cut-off of the correlation coefficient between each
#' pair of categorical variables. The first one is between each pair of 'uv.variables' and the second one is between each
#' pair of 'bio.variables'. The correlation is computed by the function ContCoef from the DescTools package. If the correlation
#' of a pair of variable is higher than the cut-off, then only the variable that has the highest number of factor will be
#' kept and the other one will be excluded from the remaining analysis. By default they are both set to 0.7.
#' @param cont.cor.coef Vector of two numerical values. Indicates the cut-off of the Spearman correlation coefficient between
#' each pair of continuous variables. The first one is between each pair of 'uv.variables' and the second one is between
#' each pair of 'bio.variables'. If the correlation of a pair of variable is higher than the cut-off, then only the variable
#' that has the highest variance will be kept and the other one will be excluded from the remaining analysis. By default
#' they are both set to 0.7.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' By default it is set to TRUE.
#' @param check.prps.connectedness TTTTT
#' @param remove.na String. Indicates whether to remove NA or missing values from either the 'assays', the 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be excluded. If 'sample.annotation' is selected, the
#' samples that contains NA or missing values for any 'bio.variables' and 'uv.variables' will be excluded. By default, it is set to
#' 'both'.
#' @param assess.variables TTTTTT
#' @param plot.output TTTTTTT
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result, by default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or
#' messages displayed during the execution of the functions, by default it is set to TRUE.

#' @return A SummarizedExperiment object or a list that contains all the PRPS sets.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr count
#' @importFrom tidyr %>%
#' @export

createSupervisedPRPS <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        apply.other.uv.variables = TRUE,
        min.sample.for.prps = 3,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 2,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 2,
        check.prps.connectedness = TRUE,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        assess.variables = FALSE,
        cat.cor.coef = c(0.7, 0.7),
        cont.cor.coef = c(0.7, 0.7),
        remove.na = 'both',
        save.se.obj = TRUE,
        plot.output = TRUE,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The supervisedPRPS function starts:',
                        color = 'white',
                        verbose = verbose)
    # define categorical and continuous variables ####
    uv.class <- sapply(uv.variables,
                       function(x) class(colData(se.obj)[[x]]))
    categorical.uv <- names(uv.class[which(uv.class %in% c('character', 'factor'))])
    continuous.uv <- uv.variables[!uv.variables %in% categorical.uv]

    # prps for categorical variables ####
    if (length(categorical.uv) > 0) {
        printColoredMessage(
            message = '-- Create PRPS for all categorical sources of unwanted variation:',
            color = 'magenta',
            verbose = verbose
        )
        if (!save.se.obj) {
            categorical.uv.prps <- lapply(
                categorical.uv,
                function(x) {
                    if (apply.other.uv.variables) {
                        other.uv.variables <- uv.variables[!uv.variables %in% x]
                    } else other.uv.variables <- NULL
                    createPRPSForCategoricalUV(
                        se.obj = se.obj,
                        assay.name = assay.name,
                        bio.variables = bio.variables,
                        main.uv.variable = x,
                        other.uv.variables = other.uv.variables,
                        min.sample.for.prps = min.sample.for.prps,
                        bio.clustering.method = bio.clustering.method,
                        nb.bio.clusters = nb.bio.clusters,
                        other.uv.clustering.method = other.uv.clustering.method,
                        nb.other.uv.clusters = nb.other.uv.clusters,
                        apply.log = apply.log,
                        pseudo.count = pseudo.count,
                        check.prps.connectedness = check.prps.connectedness,
                        assess.se.obj = FALSE,
                        save.se.obj = save.se.obj,
                        remove.na = remove.na,
                        verbose = verbose
                    )
                })
            names(categorical.uv.prps) <- categorical.uv
            categorical.uv.prps.all <- do.call(cbind, categorical.uv.prps)
        } else {
            for (x in categorical.uv) {
                if (apply.other.uv.variables) {
                    other.uv.variables <- uv.variables[!uv.variables %in% x]
                } else other.uv.variables <- NULL
                se.obj <- createPRPSForCategoricalUV(
                    se.obj = se.obj,
                    assay.name = assay.name,
                    bio.variables = bio.variables,
                    main.uv.variable = x,
                    other.uv.variables = other.uv.variables,
                    min.sample.for.prps = min.sample.for.prps,
                    bio.clustering.method = bio.clustering.method,
                    nb.bio.clusters = nb.bio.clusters,
                    other.uv.clustering.method = other.uv.clustering.method,
                    nb.other.uv.clusters = nb.other.uv.clusters,
                    apply.log = apply.log,
                    pseudo.count = pseudo.count,
                    cat.cor.coef = cat.cor.coef,
                    cont.cor.coef = cont.cor.coef,
                    check.prps.connectedness = check.prps.connectedness,
                    assess.se.obj = FALSE,
                    save.se.obj = save.se.obj,
                    remove.na = remove.na,
                    verbose = verbose
                )
            }
        }
    }
    if (length(continuous.uv) > 0) {
        printColoredMessage(message = '-- Create PRPS for all continuous sources of unwanted variation:',
                            color = 'magenta',
                            verbose = verbose)
        if (!save.se.obj) {
            continuous.uv.prps <- lapply(
                continuous.uv,
                function(x) {
                    if (apply.other.uv.variables) {
                        other.uv.variables <- uv.variables[!uv.variables %in% x]
                    } else other.uv.variables <- NULL
                    createPRPSForContinuousUV(
                        se.obj = se.obj,
                        assay.name = assay.name,
                        bio.variables = bio.variables,
                        main.uv.variable = x,
                        other.uv.variables = other.uv.variables,
                        min.sample.for.prps = min.sample.for.prps,
                        bio.clustering.method = bio.clustering.method,
                        nb.bio.clusters = nb.bio.clusters,
                        other.uv.clustering.method = other.uv.clustering.method,
                        nb.other.uv.clusters = nb.other.uv.clusters,
                        cat.cor.coef = cat.cor.coef,
                        cont.cor.coef = cont.cor.coef,
                        apply.log = apply.log,
                        assess.se.obj = FALSE,
                        save.se.obj = save.se.obj,
                        pseudo.count = pseudo.count,
                        remove.na = remove.na,
                        verbose = verbose
                    )
                })
            names(continuous.uv.prps) <- continuous.uv
            continuous.uv.prps.all <- do.call(cbind, continuous.uv.prps)
        } else{
            for (x in continuous.uv) {
                if (apply.other.uv.variables) {
                    other.uv.variables <- uv.variables[!uv.variables %in% x]
                } else other.uv.variables <- NULL
                se.obj <- createPRPSForContinuousUV(
                    se.obj = se.obj,
                    assay.name = assay.name,
                    bio.variables = bio.variables,
                    main.uv.variable = x,
                    other.uv.variables = other.uv.variables,
                    min.sample.for.prps = min.sample.for.prps,
                    bio.clustering.method = bio.clustering.method,
                    nb.bio.clusters = nb.bio.clusters,
                    other.uv.clustering.method = other.uv.clustering.method,
                    nb.other.uv.clusters = nb.other.uv.clusters,
                    cat.cor.coef = cat.cor.coef,
                    cont.cor.coef = cont.cor.coef,
                    apply.log = apply.log,
                    assess.se.obj = FALSE,
                    save.se.obj = save.se.obj,
                    pseudo.count = pseudo.count,
                    remove.na = remove.na,
                    verbose = verbose
                )
            }
        }
    }
    # save the output ####
    if (save.se.obj) {
        printColoredMessage(message = '------------The supervised.prps function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else{
        printColoredMessage(message = '------------The supervised.prps function finished.',
                            color = 'white',
                            verbose = verbose)
        if (length(continuous.uv) > 0 &
            length(categorical.uv) > 0) {
            return(
                list(
                    categorical.uv.prps = categorical.uv.prps,
                    continuous.uv.prps = continuous.uv.prps,
                    all.prps = cbind(categorical.uv.prps.all, continuous.uv.prps.all)
                )
            )
        } else if (length(continuous.uv) > 0 &
                   length(categorical.uv) == 0) {
            return(
                list(continuous.uv.prps = continuous.uv.prps,
                     all.prps = continuous.uv.prps.all)
            )
        } else if (length(continuous.uv) == 0 &
                   length(categorical.uv) > 0) {
            return(
                list(categorical.uv.prps = categorical.uv.prps,
                     all.prps = categorical.uv.prps.all)
            )
        }
    }

}
