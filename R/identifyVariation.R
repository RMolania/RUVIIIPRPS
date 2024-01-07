#' is used to identify the unwanted variation in RNA-seq data.
#'
#' Several assessment will be performed:
#' For each categorical variable:
#' - PCA plot of the categorical variable.
#' - Silhouette and ARI computed on the categorical variable.
#' - Differential analysis based ANOVA between the gene expression and the categorical variable.
#' - Vector correlation between the first cumulative PCs of the gene expression and the categorical variable.
#' For each continous variable:
#' - Linear regression between the first cumulative PC and continuous variable.
#' - Correlation between gene expression and continuous variable.
#'
#' It will output the following plots:
#' - PCA plot of each categorical variable.
#' - Boxplot of the F-test distribution from ANOVA between the gene expression and each categorical variable.
#' - Vector correlation between the first cumulative PCs of the gene expression and each categorical variable.
#' - Combined Silhouette plot of the combined pair of all categorical variables.
#' - Linear regression between the first cumulative PC and continuous variable.
#' - Boxplot of the correlation between gene expression and continuous variable.
#' - It will also output the RLE plot distribution.
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to identify the variation. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variables String or vector of strings of the label of continuous or categorical variable(s)
#' such as samples types, batch or library size from colData(se).
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#' @param output.file Path and name of the output file to save the assessments plots in a pdf format.
#' @param fast.pca logical. Indicates whether to calculate a specific number of PCs instead of the full range
#' to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to be calculated for the fast pca process, by default is set to 10.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.

#' @return  SummarizedExperiment A SummarizedExperiment object containing all the assessments plots and metrics.
#' If specified it will generate a pdf containing the assessments plots and metrics used for the assessment.

#' @importFrom gridExtra grid.arrange
#' @importFrom SummarizedExperiment assays colData
#' @export

identifyVariation <- function(
        se.obj,
        assay.names = 'All',
        variables = NULL,
        fast.pca = TRUE,
        nb.pcs = 10,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        output.file = NULL,
        verbose = TRUE
) {
    printColoredMessage(message = '------------The identifyVariation function starts.',
                        color = 'white',
                        verbose = verbose)

    se.obj <- RUVIIIPRPS::normAssessment(
        se.obj = se.obj,
        assay.names = assay.names,
        apply.log = apply.log,
        variables = variables,
        output.file = output.file,
        fast.pca = fast.pca,
        nb.pcs = nb.pcs,
        assess.se.obj = assess.se.obj,
        verbose = verbose,
        pseudo.count = pseudo.count)

    printColoredMessage(message = '------------The identifyVariation function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}
