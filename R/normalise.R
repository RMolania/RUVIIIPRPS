#' is used to normalise the gene expression (assay) of a SummarizedExperiment class object
#' using RUVIII-PRPS method.
#' The steps involves:
#' - Creation of the Pseudo-Replicates of Pseudo-Samples (PRPS)
#' - Define Negative Controls Genes (NCG)
#' - Run RUVIII-PRPS for multiple k values (the dimension of the unwanted variation).
#'
#'
#' @param se.obj A SummarizedExperiment object that will be used to computer fastRUV-III
#' @param assay.name String for the selection of the name of the assay data
#' of the SummarizedExperiment class object
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param bio.variable.prps String of the label of a categorical variable that specifies major biological groups
#' such as samples types from colData(se) that will be used to define PRPS.
#' @param uv.variables.prps String or vector of strings of the label of continuous or categorical variable(s)
#' such as samples types, batch or library size from colData(se) that will be used to define PRPS.
#' @param assess.cor.variables.prps Logical. Indicates whether to assess the assess the association between variables
#' using Spearman correlation to define PRPS.
#' @param min.sample.prps TO BE DEFINED that will be used to define PRPS.
#' @param bio.variables.ncg String or vector of strings of the label of a categorical variable that specifies major biological groups
#' such as samples types from colData(se) that will be used to find the negative controls.
#' @param uv.variables.ncg String or vector of strings of the label of continuous or categorical variable(s)
#' such as samples types, batch or library size from colData(se) that will be used to find the negative controls.
#' @param no.ncg Logical, TO BE BETTER DEFINED. if TRUE then a sample annotation the initially contains column names of the assays.???
#' @param regress.out.uv.variables.ncg TO BE DEFINED.
#' @param regress.out.bio.variables.ncg TO BE DEFINED.
#' @param apply.normalization.ncg Logical Indicates whether to apply a normalization method when providing the raw data assay to define NCG.
#' By default it is set to FALSE.
#' @param normalization.ncg String defining the normalization method to use from 'CPM', 'TMM', 'upper', 'full', 'median', 'VST',
#' and 'Quantile'to define NCG. By default it is set to 'CPM'.
#' @param k A single value or a vector of values containing a single k or a range of k - the number of unwanted factors - to be tested.
#' @param eta A matrix with n columns, gene-wise (as opposed to sample-wise) covariates. By default is set to NULL.
#' @param include.intercept Logical. Add an intercept term to eta if it does not include one already. By default is set to TRUE.
#' @param apply.average.rep Average replicates after adjustment. By default is set to FALSE.
#' @param fullalpha Can be included to speed up execution. By default is set to NULL.
#' @param inputcheck Perform a basic sanity check on the inputs, and issue a warning if there is a problem.
#' By default is set to TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' By default it is set to TRUE.
#' @param remove.na TO BE DEFINED.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the normalised gene expression
#' into a new assay (RUVIII_K), as well as the PRPS created and the negative control genes

#' @export
#'
#'
normalise <- function(
        se.obj,
        assay.name,
        apply.log = TRUE,
        pseudo.count = 1,
        bio.variable.prps,
        uv.variables.prps,
        min.sample.prps = 3,
        assess.cor.variables.prps = TRUE,
        bio.variables.ncg,
        uv.variables.ncg,
        no.ncg = 1000,
        regress.out.uv.variables.ncg = FALSE,
        regress.out.bio.variables.ncg = FALSE,
        apply.normalization.ncg=FALSE,
        normalization.ncg = 'CPM',
        k = NULL,
        eta = NULL,
        include.intercept = TRUE,
        apply.average.rep = FALSE,
        fullalpha = NULL,
        inputcheck = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'both',
        verbose = TRUE
) {
    printColoredMessage(message = '------------The normalise function starts.',
                        color = 'white',
                        verbose = verbose)

    ########### Creation of PRPS ###########
    se.obj= supervisedPRPS(se.obj=se.obj,
                       assay.name=assay.name,
                       bio.variable= bio.variable.prps,
                       uv.variables=uv.variables.prps,
                       min.sample.prps = min.sample.prps,
                       assess.cor.variables = assess.cor.variables.prps,
                       save.se.obj = TRUE,
                       verbose = verbose)


    ############## NCG ####################
    se.obj=supervisedFindNGC(se.obj=se.obj,
                         assay.name=assay.name,
                         bio.variables= bio.variables.ncg,
                         uv.variables= uv.variables.ncg,
                         no.ncg = no.ncg,
                         regress.out.uv.variables =regress.out.uv.variables.ncg,
                         regress.out.bio.variables = regress.out.bio.variables.ncg,
                         apply.normalization=apply.normalization.ncg,
                         normalization =  normalization.ncg,
                         save.se.obj = TRUE,
                         verbose = verbose)

    ############## RUVIII-PRPS ####################
    replicate.data=t(do.call(cbind,se.obj@metadata[['PRPS']][['supervised']]))
    ruvIIIMultipleK(se.obj=se.obj,
                    assay.name=assay.name,
                    replicate.data=replicate.data,
                    ctl=se.obj@metadata[['NCG']],
                    k = k,
                    eta = eta,
                    include.intercept = include.intercept,
                    apply.average.rep = apply.average.rep,
                    fullalpha = fullalpha,
                    return.info = FALSE,
                    inputcheck = inputcheck,
                    assess.se.obj = assess.se.obj,
                    remove.na = 'measurements',
                    save.se.obj = TRUE,
                    verbose = verbose
    )

    printColoredMessage(message = '------------The normalise function finished.',
                        color = 'white',
                        verbose = verbose)
}

