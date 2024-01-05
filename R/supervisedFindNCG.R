#' is used to find a set of negative control genes  (NCG) of an assay in a SummarizedExperiment object.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name A name of an assay data in the SummarizedExperiment object.
#' @param approach TTTT
#' @param nb.ncg TTTTT
#' @param top.rank.bio.genes TTTTT
#' @param top.rank.uv.genes TTTTT
#' @param bio.variables TTTT
#' @param bio.groups TTTT
#' @param nb.bio.clusters TTTTT
#' @param uv.variables YYYY
#' @param uv.groups YYYY
#' @param nb.uv.clusters YYYYY
#' @param bio.clustering.method YYYY
#' @param uv.clustering.method UUUUU
#' @param normalization HHHH
#' @param regress.out.uv.variables UUUUU
#' @param regress.out.bio.variables UUUUUU
#' @param apply.log UUUUU
#' @param pseudo.count YYYYY
#' @param min.sample.for.aov UUUUU
#' @param min.sample.for.correlation UUUUU
#' @param corr.method Logical. Indicates whether to apply a log-transformation to the data.
#' @param a Logical. Indicates whether to apply a log-transformation to the data.
#' @param rho Logical. Indicates whether to apply a log-transformation to the data.
#' @param grid.nb NNNNNN
#' @param ncg.selection.method PPPP
#' @param assess.ncg PPPP
#' @param nb.pcs PPPPPP
#' @param center PP
#' @param scale PP
#' @param anova.method PPPPP
#' @param assess.variables PPPP
#' @param remove.na PPPPP
#' @param save.se.obj PPPP
#' @param cat.cor.coef PPPP
#' @param cont.cor.coef PPPP
#' @param variables.to.assess.ncg PPPP
#' @param assess.se.obj Logical. Indicates whether to apply a log-transformation to the data.
#' @param verbose Logical. Indicates whether to apply a log-transformation to the data.

#' @return a list of negative control genes.

#' @importFrom Matrix colSums
#' @importFrom dplyr left_join
#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @importFrom biomaRt getBM useMart useDataset
#' @importFrom S4Vectors DataFrame
#' @export

supervisedFindNCG <- function(
    se.obj,
    assay.name,
    approach = 'PbPbio',
    nb.ncg = 10,
    top.rank.bio.genes = 50,
    top.rank.uv.genes = 50,
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
    grid.nb = 10,
    anova.method = 'aov',
    ncg.selection.method = 'Prod',
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
    if(!approach %in% c('PbPb', 'AnoCorrAs', 'TWAnova')){
        stop('The approach must be one of the "PbPb", "AnoCorrAs" or "TWAnova".')
    }

    if(approach == 'PbPb'){
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
    printColoredMessage(message = '------------The supervisedFindNCG function finished:',
                        color = 'white',
                        verbose = verbose)
    return(ncg.set)
}
