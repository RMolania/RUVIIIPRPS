#' is used to find a set of negative control genes  (NCG) of an assay in a SummarizedExperiment object.

#' @author Ramyar Molania

#' @description
#' This function uses the two-way ANOVA approach to find a set of genes as negative control genes (NCG) for RUV-III-PRPS
#' normalization. Sources of known biological and unwanted variation should be specified. First, the function creates all
#' possible sample groups with respect to biological and unwanted variation separately. Then, the function uses the groups
#' as factors in two-way ANOVA to find genes that are highly affected by biological and unwanted variation separately.
#' Finally, the function selects genes that have possible high F-statistics for the unwanted variation and low F-statistics
#' for the biological variation. The function uses different approaches to perform the final selection.


#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. Indicates a name of assay in the SummarizedExperiment object to be used to find a set of
#' negative control genes.
#' @param approach Symbol. A symbol indicates which method should be used for selection of negative control genes.
#' @param ncg.selection.method Symbol.Indicates how to select a set genes as NCG. For individual genes, the two-way ANOVA
#' calculates F-statistics for biological and unwanted variation factors separately. An ideal NCG set should have high
#' F-statistics for the unwanted variation variables and low F-statistics for the biological variables. The function ranks
#' the F-statistics obtained for the biological variable and negative of the F-statistics obtained for the unwanted variables.
#' Then this functions offers 5 ways to summarize the ranks of the two F-statistics. Prod' is the product of the ranks.
#' 'Sum', is the sum of the ranks. 'Average' is the average of the ranks. 'AbsNoneOverlap' is the none overlapped genes
#' of the 'top.rank.uv.genes' and 'top.rank.bio.genes'. 'noneOverlap' is the none overlapped genes of the 'top.rank.uv.genes'
#' and at least 'top.rank.bio.genes'. The F-statistics for biological and UV are first ranked.Then options are Prod
#' (product), Sum, Average.
#' @param bio.variables Symbols. Indicates the columns names that contain biological variables in the
#' SummarizedExperiment object.
#' @param uv.variables Symbols. Indicates the columns names that contains UV variables in the SummarizedExperiment object.
#' @param nb.ncg Numeric. Indicates the percentage of the total genes to be selected as a NCG set. The default is 10 percent.
#' @param grid.nb Numeric. Indicates the percentage for grid search when the ncg.selection.method is 'noneOverlap'. In the
#' 'noneOverlap' approach, the grid search starts with the initial top.rank.uv.genes' value an add the grid.nb in each loop
#' to find the 'nb.ncg'.
#' @param top.rank.bio.genes Numeric.Indicates the percentage of top ranked genes that are highly affected by the biological
#' variation. This is required to be specified when the 'ncg.selection.method' is either 'noneOverlap' or 'AbsNoneOverlap'.
#' @param top.rank.uv.genes Numeric.Indicates the percentage of top ranked genes that are highly affected by the unwanted
#' variation variables. This is required to be specified when the 'ncg.selection.method' is either 'noneOverlap' or
#' 'AbsNoneOverlap'.
#' @param bio.groups Symbols. a symbol or a vector of symbols indicating the columns names that contains biological variables
#' in the SummarizedExperiment object. If is not NULL, the 'bio.groups' will be used for grouping samples into different
#' homogeneous biological groups.
#' @param bio.clustering.method Symbols. Indicates which clustering methods should be used to group continuous sources of
#' biological variation if any is provided. The default is kmeans clustering.
#' @param nb.bio.clusters Numeric. Indicates the number of clusters for each continuous sources of biological variation.
#' The default is 3.
#' @param uv.groups Symbols. a symbol or a vector of symbols indicating the columns names that contains unwanted variation
#' variables in the SummarizedExperiment object. If is not NULL, the 'bio.groups' will be used for grouping samples into
#' different homogeneous biological groups.
#' @param uv.clustering.method Symbols.Indicates which clustering method should be used to group continuous sources of UV.
#' The default is kmeans clustering.
#' @param nb.uv.clusters Numeric. Indicates the number of clusters for each continuous sources of UV, by default it is
#' set to 2.
#' @param normalization Symbol. Indicates which normalization method should be applied to the data before finding genes
#' that are affected by biological variation. The default is CPM. We refer to the 'applyOtherNormalizations' function for
#' more details.
#' @param regress.out.bio.variables Symbols. Indicates the columns names that contain biological variables in the
#' SummarizedExperiment object. These variables will be regressed out from the data before finding genes that are highly
#' affected by unwanted variation variable. The default is NULL, indicates the regression will not be applied.
#' @param regress.out.uv.variables Symbols. Indicates the columns names that contain unwanted variation variables in the
#' SummarizedExperiment object. These variables will be regressed out from the data before finding genes that are highly
#' affected by biological variation. The default is NULL, indicates the regression will not be applied.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before perforing correlation and
#' ANOVA. The default is TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation.
#' @param min.sample.for.aov Numeric. Indicates the minimum number of samples to be present in each group before applying
#' the ANOVA. The default is 3.
#' @param min.sample.for.correlation Numeric. Indicates the minimum number of samples to be considered in each group before
#' applying the correlation analysis. The default is 10.
#' @param corr.method Symbol. Indicates which correlation method should be used to compute association between gene-level
#' expression and a continuous variable. The default is 'spearman'.
#' @param a The significance level used for the confidence intervals in the correlation, by default it is set to 0.05.
#' @param rho The value of the hypothesised correlation to be used in the hypothesis testing, by default it is set to 0.
#' @param anova.method Indicates which anova method should be used to compute association between gene-level
#' expression and a categorical variable. The default is 'aov'.
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected NCGs or not. This analysis involves
#' principal component analysis on the selected NCG and then explore the R^2 or vector correlation between the 'nb.pcs'
#' first principal components and with biological and unwanted variables.
#' @param variables.to.assess.ncg Symbols. Indicates the column names of the SummarizedExperiment object that contain
#' variables whose association with the selected genes as NCG. needs to be evaluated. The default is NULL. This means all
#' the variables in the 'bio.variables' and 'uv.variables' will be assessed.
#' @param nb.pcs Numeric. Indicates the number of the first principal components on selected NCG to be used to assess the
#' performance of NCGs.
#' @param center Logical. Indicates whether to center the data before applying principal component analysis or not.
#' The default is TRUE.
#' @param scale Logical. Indicates whether to scale the data before applying Principal component analysis. The default
#' is FALSE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object or not.
#' @param assess.variables Logical.Indicates whether to assess the SummarizedExperiment class object or not.
#' @param cat.cor.coef Vector of two numerical values. Indicates the cut-off of the correlation coefficient between each pair of
#' categorical variables. The first one is between each pair of 'uv.variables' and the second one is between each pair of 'bio.variables'.
#' The correlation is computed by the function ContCoef from the DescTools package. If the correlation of a pair of variable is higher than
#' the cut-off, then only the variable that has the highest number of factor will be kept and the other one will be excluded from the
#' remaining analysis. By default they are both set to 0.7.
#' @param cont.cor.coef Vector of two numerical values. Indicates the cut-off of the Spearman correlation coefficient between each pair of
#' continuous variables. The first one is between each pair of 'uv.variables' and the second one is between each pair of 'bio.variables'.
#' If the correlation of a pair of variable is higher than the cut-off, then only the variable that has the highest variance will
#' be kept and the other one will be excluded from the remaining analysis. By default they are both set to 0.7.
#' @param remove.na Symbol. Indicates whether to remove NA or missing values from either the 'assays', the 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be excluded. If
#' 'sample.annotation' is selected, the samples that contains NA or missing values for any 'bio.variables' and 'uv.variables'
#' will be excluded. By default, it is set to both'.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object
#' 'se.obj' or to output the result. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return a list of negative control genes.

#' @importFrom Matrix colSums
#' @importFrom dplyr left_join
#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @importFrom biomaRt getBM useMart useDataset
#' @importFrom S4Vectors DataFrame
#' @export

FindNCGsSupervised <- function(
    se.obj,
    assay.name,
    approach = 'TWAnova',
    ncg.selection.method = 'noneOverlap',
    bio.variables,
    uv.variables,
    nb.ncg = 10,
    grid.nb = 1,
    top.rank.bio.genes = 50,
    top.rank.uv.genes = 50,
    bio.groups = NULL,
    nb.bio.clusters = 3,
    bio.clustering.method = 'kmeans',
    uv.groups = NULL,
    nb.uv.clusters = 3,
    uv.clustering.method = 'kmeans',
    normalization = 'CPM',
    regress.out.uv.variables = NULL,
    regress.out.bio.variables = NULL,
    apply.log = TRUE,
    pseudo.count = 1,
    min.sample.for.aov = 3,
    min.sample.for.correlation = 10,
    corr.method = "spearman",
    a = 0.05,
    rho = 0,
    anova.method = 'aov',
    assess.ncg = TRUE,
    variables.to.assess.ncg = NULL,
    nb.pcs = 5,
    center = TRUE,
    scale = FALSE,
    assess.se.obj = TRUE,
    assess.variables = TRUE,
    cat.cor.coef = c(0.95, 0.95),
    cont.cor.coef = c(0.95, 0.95),
    remove.na = 'both',
    save.se.obj = TRUE,
    verbose = TRUE
    ){
    printColoredMessage(message = '------------The supervisedFindNCG function starts:',
                        color = 'white',
                        verbose = verbose)
    if(!approach %in% c('PbPbio', 'AnoCorrAs', 'TWAnova')){
        stop('The approach must be one of the "PbPbio", "AnoCorrAs" or "TWAnova".')
    }

    if(approach == 'PbPbio'){
        ncg.set <- supervisedFindNcgPbPbio(
            se.obj = se.obj,
            assay.name = assay.name,
            nb.ncg = nb.ncg,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            bio.variables = bio.variables,
            bio.groups = bio.groups,
            nb.bio.clusters = nb.bio.clusters,
            uv.variables = uv.variables,
            uv.groups = uv.groups,
            nb.uv.clusters = nb.uv.clusters,
            bio.clustering.method = bio.clustering.method,
            uv.clustering.method = uv.clustering.method,
            normalization = normalization,
            regress.out.uv.variables = regress.out.uv.variables,
            regress.out.bio.variables = regress.out.bio.variables,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            min.sample.for.aov = min.sample.for.aov,
            min.sample.for.correlation = min.sample.for.correlation,
            corr.method = corr.method,
            a = a,
            rho = rho,
            ncg.selection.method = ncg.selection.method,
            assess.ncg = assess.ncg,
            variables.to.assess.ncg = variables.to.assess.ncg,
            nb.pcs = nb.pcs,
            anova.method = anova.method,
            assess.se.obj = assess.se.obj,
            assess.variables = assess.variables,
            remove.na = remove.na,
            save.se.obj = save.se.obj,
            verbose = verbose)
    } else if(approach == 'AnoCorrAs'){
        ncg.set <- supervisedFindNcgAnoCorrAs(
            se.obj = se.obj,
            assay.name = assay.name,
            bio.variables = bio.variables,
            uv.variables = uv.variables,
            nb.ncg = nb.ncg,
            ncg.selection.method = ncg.selection.method,
            grid.nb = grid.nb,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            min.sample.for.aov = min.sample.for.aov,
            min.sample.for.correlation = min.sample.for.correlation,
            regress.out.uv.variables = regress.out.uv.variables,
            regress.out.bio.variables = regress.out.bio.variables,
            normalization = normalization,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            corr.method = corr.method,
            a = a,
            rho = rho,
            anova.method = anova.method,
            assess.ncg = assess.ncg,
            variables.to.assess.ncg = variables.to.assess.ncg,
            nb.pcs = nb.pcs,
            scale = scale,
            center = center,
            assess.se.obj = assess.se.obj,
            assess.variables = assess.variables,
            cat.cor.coef = cat.cor.coef,
            cont.cor.coef = cont.cor.coef,
            save.se.obj = save.se.obj,
            remove.na = remove.na,
            verbose = verbose)
    } else if(approach == 'TWAnova'){
        ncg.set <- supervisedFindNcgTWAnova(
            se.obj = se.obj,
            assay.name = assay.name,
            bio.variables = bio.variables,
            uv.variables = uv.variables,
            nb.ncg = nb.ncg,
            ncg.selection.method = ncg.selection.method,
            grid.nb = grid.nb,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            bio.clustering.method = bio.clustering.method,
            nb.bio.clusters = nb.bio.clusters,
            uv.clustering.method = uv.clustering.method,
            nb.uv.clusters = nb.uv.clusters,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            assess.ncg = assess.ncg,
            variables.to.assess.ncg = variables.to.assess.ncg,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            assess.se.obj = assess.se.obj,
            assess.variables = assess.variables,
            remove.na = remove.na,
            save.se.obj = save.se.obj,
            verbose = verbose)
    } else {
        stop('The approach should be one of PbPb or As')
    }
    printColoredMessage(message = '------------The supervisedFindNCG function finished.',
                        color = 'white',
                        verbose = verbose)
    return(ncg.set)
}
