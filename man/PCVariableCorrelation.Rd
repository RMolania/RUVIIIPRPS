% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCVariableCorrelation.R
\name{PCVariableCorrelation}
\alias{PCVariableCorrelation}
\title{is used to compute the vector correlation between the first cumulative PCs of the gene expression (assay)
of a SummarizedExperiment class object and a categorical variable (i.e. batch).}
\usage{
PCVariableCorrelation(
  se.obj,
  assay.names = "All",
  variable,
  fast.pca = TRUE,
  nb.pcs = 10,
  save.se.obj = TRUE,
  plot.output = TRUE,
  assess.se.obj = TRUE,
  remove.na = "both",
  apply.round = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{se.obj}{A SummarizedExperiment object that will be used to compute the PCA.}

\item{assay.names}{Optional string or list of strings for the selection of the name(s)
of the assay(s) of the SummarizedExperiment class object to compute the correlation. By default}

\item{variable}{String of the label of a categorical variable such as
sample types or batches from colData(se.obj).}

\item{fast.pca}{Logical. Indicates whether to use the PCA calculated using a specific number of PCs instead of the full range
to speed up the process, by default is set to 'TRUE'.}

\item{nb.pcs}{Numeric. The number of few first cumulative PCs, by default is set to 10.}

\item{save.se.obj}{Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
to output the result. By default it is set to TRUE.}

\item{plot.output}{Logical. Indicates whether to plot the correlation statistics, by default it is set to TRUE.}

\item{assess.se.obj}{Logical. Indicates whether to assess the SummarizedExperiment class object.}

\item{remove.na}{TO BE DEFINED.}

\item{apply.round}{Logical. Indicates whether to round the ARI results, by default it is set to TRUE.}

\item{verbose}{Indicates whether to show or reduce the level of output or messages displayed during the execution
of the functions, by default it is set to TRUE.}
}
\value{
SummarizedExperiment A SummarizedExperiment object containing the computed correlation for
the continuous variable and if requested the associated plot.
}
\description{
is used to compute the vector correlation between the first cumulative PCs of the gene expression (assay)
of a SummarizedExperiment class object and a categorical variable (i.e. batch).
}
