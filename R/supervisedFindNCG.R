#' is used to find a set of negative control genes  (NCG) of an assay in a SummarizedExperiment object.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name A name of an assay data in the SummarizedExperiment object.
#' @param bio.variables numeric. In large sample situations, the minimum proportion of samples in a group that a gene needs to be expressed in. See Details below for the exact formula.
#' @param uv.variables logical, if TRUE then library size is calculated on the raw.count.assay.name.
#' @param no.ncg logical, if TRUE then a sample annotation the initially contains column names of the assays.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data.
#' @param pseudo.count Logical. Indicates whether to apply a log-transformation to the data.
#' @param regress.out.uv.variables Logical. Indicates whether to apply a log-transformation to the data.
#' @param regress.out.bio.variables Logical. Indicates whether to apply a log-transformation to the data.
#' @param normalization Logical. Indicates whether to apply a log-transformation to the data.
#' @param corr.method Logical. Indicates whether to apply a log-transformation to the data.
#' @param a Logical. Indicates whether to apply a log-transformation to the data.
#' @param rho Logical. Indicates whether to apply a log-transformation to the data.
#' @param anova.method Logical. Indicates whether to apply a log-transformation to the data.
#' @param assess.se.obj Logical. Indicates whether to apply a log-transformation to the data.
#' @param verbose Logical. Indicates whether to apply a log-transformation to the data.

#' @return a list of negative control genes.

#' @importFrom MAtrix colSums
#' @importFrom dplyr left_join
#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @importFrom biomaRt getBM useMart useDataset
#' @importFrom S4Vectors DataFrame transpose
#' @export

supervisedFindNCG <- function(
    se.obj,
    assay.name,
    approach = 'PbPb',
    nb.ncg = 0.1,
    top.rank.bio.genes = .3,
    top.rank.uv.genes = 1,
    bio.variables,
    bio.groups = NULL,
    nb.bio.clusters = 2,
    uv.variables,
    uv.groups = NULL,
    nb.uv.clusters = 2,
    bio.clustering.method = 'kmeans',
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
    grid.no = 10,
    ncg.selection.method = 'Prod',
    assess.ncg = TRUE,
    variables.to.assess.ncg = NULL,
    nb.pcs = 5,
    anova.method = 'aov',
    assess.se.obj = TRUE,
    assess.variables = TRUE,
    remove.na = 'both',
    save.se.obj = TRUE,
    verbose = TRUE
    ){
    printColoredMessage(message = '------------The supervisedFindNCG function starts:',
                        color = 'white',
                        verbose = verbose)
    if(!approach %in% c('PbPb', 'AnoCorrAs', 'TWAnova')){
        stop('The approach must be one of the "PbPb", "AnoCorrAs" or "TWAnova".')
    }

    if(approach == 'PbPb'){
        ncg.set <- supervisedFindNcgPbPb(
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
            ncg.selection.method = ncg.selection.method,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            nb.ncg = nb.ncg,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            min.sample.for.aov = min.sample.for.aov,
            min.sample.for.correlation = min.sample.for.correlation,
            regress.out.uv.variables = regress.out.uv.variables,
            regress.out.bio.variables = regress.out.bio.variables,
            normalization = normalization,
            corr.method = corr.method,
            a = a,
            rho = rho,
            variables.to.assess.ncg = variables.to.assess.ncg,
            anova.method = anova.method,
            assess.se.obj = assess.se.obj,
            assess.variables = assess.variables,
            remove.na = remove.na,
            assess.ncg = assess.ncg,
            nb.pcs = nb.pcs,
            save.se.obj = save.se.obj,
            verbose = verbose)
    } else if(approach == 'TWAnova'){
        ncg.set <- supervisedFindNcgTWAnova(
            se.obj = se.obj,
            assay.name = assay.name,
            nb.ncg = nb.ncg,
            ncg.selection.method = ncg.selection.method,
            grid.no = grid.no,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            bio.variables = bio.variables,
            bio.clustering.method = bio.clustering.method,
            nb.bio.clusters = nb.bio.clusters,
            uv.variables = uv.variables,
            nb.uv.clusters = nb.uv.clusters,
            uv.clustering.method = uv.clustering.method,
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
    printColoredMessage(message = '------------The supervisedFindNCG function finished:',
                        color = 'white',
                        verbose = verbose)
    return(ncg.set)
}
