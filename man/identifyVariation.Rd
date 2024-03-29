% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/identifyVariation.R
\name{identifyVariation}
\alias{identifyVariation}
\title{is used to identify the unwanted variation of a SummarizedExperiment class object.}
\usage{
identifyVariation(
  se.obj,
  assay.names = "All",
  apply.log = TRUE,
  variables = NULL,
  output.file = NULL,
  fast.pca = TRUE,
  nb.pcs = 10,
  assess.se.obj = TRUE,
  verbose = TRUE,
  pseudo.count = 1
)
}
\arguments{
\item{se.obj}{A SummarizedExperiment object that will be used to compute the PCA.}

\item{assay.names}{Optional string or list of strings for the selection of the name(s)
of the assay(s) of the SummarizedExperiment class object to identify the variation. By default}

\item{apply.log}{Logical. Indicates whether to apply a log-transformation to the data. By default
no transformation will be selected.}

\item{variables}{String or vector of strings of the label of continuous or categorical variable(s)
such as samples types, batch or library size from colData(se).}

\item{output.file}{Path and name of the output file to save the assessments plots in a pdf format.}

\item{fast.pca}{logical. Indicates whether to calculate a specific number of PCs instead of the full range
to speed up the process, by default is set to 'TRUE'.}

\item{nb.pcs}{Numeric. The number of first PCs to be calculated for the fast pca process, by default is set to 10.}

\item{assess.se.obj}{Logical. Indicates whether to assess the SummarizedExperiment class object.}

\item{verbose}{Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
of the functions, by default it is set to TRUE.}

\item{pseudo.count}{Numeric. A value as a pseudo count to be added to all measurements before log transformation,
by default it is set to 1.}
}
\value{
SummarizedExperiment A SummarizedExperiment object containing all the assessments plots and metrics.
If specified it will generate a pdf containing the assessments plots and metrics used for the assessment.
}
\description{
Several assessment will be performed:
For each categorical variable:
\itemize{
\item PCA plot of the categorical variable.
\item Silhouette and ARI computed on the categorical variable.
\item Differential analysis based ANOVA between the gene expression and the categorical variable.
\item Vector correlation between the first cumulative PCs of the gene expression and the categorical variable.
For each continous variable:
\item Linear regression between the first cumulative PC and continuous variable.
\item Correlation between gene expression and continuous variable.
}
}
\details{
It will output the following plots:
\itemize{
\item PCA plot of each categorical variable.
\item Boxplot of the F-test distribution from ANOVA between the gene expression and each categorical variable.
\item Vector correlation between the first cumulative PCs of the gene expression and each categorical variable.
\item Combined Silhouette plot of the combined pair of all categorical variables.
\item Linear regression between the first cumulative PC and continuous variable.
\item Boxplot of the correlation between gene expression and continuous variable.
\item It will also output the RLE plot distribution.
}
}
