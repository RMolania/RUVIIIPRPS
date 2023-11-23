#' is used to normalise the gene expression (assay) of a SummarizedExperiment class object
#' using RUVIII-PRPS method.
#'
#' The steps involves:
#' - Creation of the Pseudo-Replicates of Pseudo-Samples (PRPS)
#' - Define Negative Controls Genes (NCG)
#' - Run RUVIII-PRPS for multiple k values (the dimension of the unwanted variation).
#'
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the RUV-III.
#' @param assay.name String for the selection of the name of the assay data of the SummarizedExperiment class object
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param bio.variable.prps String of the label of a categorical variable that specifies major biological groups
#' such as samples types from colData(se) that will be used to define PRPS.
#' @param uv.variables.prps String or vector of strings of the label of continuous or categorical variable(s)
#' such as samples types, batch or library size from colData(se) that will be used to define PRPS.
#' @param batch.variable.prps String of the label of a categorical variable that specifies major batch groups
#' such as plates from colData(se).
#' @param assess.cor.variables.prps Logical. Indicates whether to assess the association between pairs of categorical variables
#' and/or pairs of continuous variables in order to select only one of the variable of a pair of highly correlated
#' variables.
#' @param min.sample.for.prps Numeric. Indicates the minimum number of samples to create one pseudo-sample,
#' by default it is set to 3.
#' @param min.sample.per.batch.prps Numeric. Indicates the minimum number of homogeneous biological samples within each batch
#' to create a PRPS set for each continuous variable. The minimum should be '2*min.sample.for.prps'. By default it is set to 6.
#' @param norm.assay.name.ncg String for the selection of the name of the assay of the SummarizedExperiment class object to use
#' to define NCG. If you don't provide any assay, we recommend to set apply.normalization to TRUE.
#' @param apply.normalization.ncg Logical Indicates whether to apply a normalization method to define NCG if the 'norm.assay.name.ncg'
#' wasn't provided. By default it is set to FALSE.
#' @param normalization.ncg String defining the normalization method to use from 'CPM', 'TMM', 'upper', 'full', 'median', 'VST',
#' and 'Quantile'to define NCG. By default it is set to 'CPM'.
#' @param bio.variables.ncg String or vector of strings of the label of a categorical variable that specifies major biological groups
#' such as samples types from colData(se) that will be used to find the negative controls.
#' @param uv.variables.ncg String or vector of strings of the label of continuous or categorical variable(s)
#' such as samples types, batch or library size from colData(se) that will be used to find the negative controls.
#' @param no.ncg Logical, TO BE BETTER DEFINED. if TRUE then a sample annotation the initially contains column names of the assays.???
#' @param regress.out.uv.variables.ncg TO BE DEFINED.
#' @param regress.out.bio.variables.ncg TO BE DEFINED.
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
        batch.variable.prps,
        min.sample.for.prps = 3,
        min.sample.per.batch.prps=6,
        assess.cor.variables.prps = FALSE,
        norm.assay.name.ncg,
        apply.normalization.ncg=FALSE,
        normalization.ncg = 'CPM',
        bio.variables.ncg,
        uv.variables.ncg,
        no.ncg = 1000,
        regress.out.uv.variables.ncg = FALSE,
        regress.out.bio.variables.ncg = FALSE,
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

    ### Assess the input
    if(is.null(k)){
        stop('k cannot be 0. This means no adjustment will be made.')
    } else if(is.null(assay.name)){
        stop('No assay name has been provided.')
    }
    if (length(assay.name) > 1) {
        stop('The function can only take a single assay.name.')
    }

    if (missing(norm.assay.name.ncg)) {
        stop('norm.assay.name.ncg is empty. Please select an assay to use for the Negative Control Genes selection or
        set norm.assay.name.ncg to NULL and set apply.normalization.ncg to TRUE.')
    }


    ########### Creation of PRPS ###########
    se.obj= supervisedPRPS(se.obj=se.obj,
                       assay.name=assay.name,
                       bio.variable= bio.variable.prps,
                       uv.variables=uv.variables.prps,
                       batch.variable = batch.variable.prps,
                       min.sample.for.prps = min.sample.for.prps,
                       min.sample.per.batch=min.sample.per.batch.prps,
                       apply.log = apply.log,
                       pseudo.count = pseudo.count,
                       assess.se.obj = assess.se.obj,
                       assess.cor.variables = assess.cor.variables.prps,
                       save.se.obj = TRUE,
                       verbose = verbose)


    ############## NCG ####################
    se.obj=supervisedFindNCG(se.obj=se.obj,
                         assay.name=norm.assay.name.ncg,
                         bio.variables= bio.variables.ncg,
                         uv.variables= uv.variables.ncg,
                         no.ncg = no.ncg,
                         regress.out.uv.variables =regress.out.uv.variables.ncg,
                         regress.out.bio.variables = regress.out.bio.variables.ncg,
                         apply.normalization=apply.normalization.ncg,
                         normalization =  normalization.ncg,
                         assess.se.obj = assess.se.obj,
                         apply.log = apply.log,
                         pseudo.count = pseudo.count,
                         save.se.obj = TRUE,
                         verbose = verbose)

    ############## RUVIII-PRPS ####################
    replicate.data=t(do.call(cbind,se.obj@metadata[['PRPS']][['supervised']]))
    se.obj=ruvIIIMultipleK(se.obj=se.obj,
                    assay.name=assay.name,
                    apply.log=apply.log,
                    pseudo.count = pseudo.count,
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
    return(se.obj)
}

