% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genesVariableCorrelation.R
\name{genesVariableCorrelation}
\alias{genesVariableCorrelation}
\title{is used to compute the correlation between the gene expression (assay)
of a SummarizedExperiment class object and a continuous variable (i.e. library size).}
\usage{
genesVariableCorrelation(
  se.obj,
  assay.names = "All",
  variable,
  apply.log = TRUE,
  method = "spearman",
  a = 0.05,
  rho = 0,
  save.se.obj = TRUE,
  boxplot.output = TRUE,
  plot.top.genes = FALSE,
  nb.top.genes = 3,
  assess.se.obj = TRUE,
  remove.na = "both",
  verbose = TRUE,
  pseudo.count = 1,
  apply.round = TRUE
)
}
\arguments{
\item{se.obj}{A SummarizedExperiment object that will be used to compute the correlation.}

\item{assay.names}{Optional string or list of strings for the selection of the name(s)
of the assay(s) of the SummarizedExperiment class object to compute the correlation. By default}

\item{variable}{String of the label of a continuous variable such as
library size from colData(se.obj).}

\item{apply.log}{Logical. Indicates whether to apply a log-transformation to the data. By default
the log transformation will be selected.}

\item{method}{A character string indicating which correlation coefficient
is to be used for the test: "pearson", "kendall", or "spearman". By default "spearman" will
be selected.}

\item{a}{The significance level used for the confidence intervals in the correlation,
by default it is set to 0.05.}

\item{rho}{The value of the hypothesised correlation to be used in the hypothesis testing,
by default it is set to 0.}

\item{save.se.obj}{Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
to output the result. By default it is set to TRUE.}

\item{boxplot.output}{Logical. Indicates whether to plot the boxplot of the correlation, by default it is set to TRUE.}

\item{plot.top.genes}{Logical. Indicates whether to plot the gene expression of the number of genes
from the high or low correlation, by default it is set to FALSE.}

\item{nb.top.genes}{Defines the number of genes from the high or low correlation to plot,
by default is set to 3.}

\item{assess.se.obj}{Logical. Indicates whether to assess the SummarizedExperiment class object.}

\item{remove.na}{TO BE DEFINED.}

\item{verbose}{Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
of the functions, by default it is set to TRUE.}

\item{pseudo.count}{Numeric. A value as a pseudo count to be added to all measurements before log transformation,
by default it is set to 1.}

\item{apply.round}{Logical. Indicates whether to round the ARI results, by default it is set to TRUE.}
}
\value{
SummarizedExperiment A SummarizedExperiment object containing the computed correlation on the continuous variable
and if requested the associated plot.
}
\description{
is used to compute the correlation between the gene expression (assay)
of a SummarizedExperiment class object and a continuous variable (i.e. library size).
}
