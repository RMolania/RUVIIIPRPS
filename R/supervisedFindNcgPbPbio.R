#' is used to find a set of negative control genes using ANOVA and correlation analyses per homogeneous batches and per homogeneous biological groups.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol.Indicates a name of assay in the SummarizedExperiment object.
#' @param nb.ncg Numeric. Indicates the percentage of the total genes to be selected a NCG set,
#' by default it is set to 10 percent.
#' @param ncg.selection.method Symbol.Indicates how to select a set genes as NCG.
#' For individual genes, the two-way ANOVA calculates F-statistics for biological and
#' unwanted variation factors separately. An ideal NCG set should have high F-statistics
#' for the unwanted variation variables and low F-statistics for the biological variables.
#' The function ranks the F-statistics obtained for the biological variable and negative of
#' the F-statistics obtained for the unwanted variables. Then this functions offers 5 ways to
#' summarize the ranks of the two F-statistics. Prod' is the product of the ranks. 'Sum', is
#' the sum of the ranks. 'Average' is the average of the ranks. 'AbsNoneOverlap' is the none overlapped
#' genes of the 'top.rank.uv.genes' and 'top.rank.bio.genes'. 'noneOverlap' is the none overlapped genes
#' of the 'top.rank.uv.genes' and at least 'top.rank.bio.genes'. The F-statistics for biological and UV
#' are first ranked.Then options are Prod (product), Sum, Average,
#' @param grid.nb Numeric. Indicates the percentage for grid search when the ncg.selection.method is
#' 'noneOverlap'. In the 'noneOverlap' approach, the grid search starts with the initial
#' 'top.rank.uv.genes' value an add the grid.nb in each loop to find the 'nb.ncg'.
#' @param top.rank.bio.genes Numeric.Indicates the percentage of top ranked genes that are highly affected
#' by the biological variation. This is required to be specified when the 'ncg.selection.method' is
#' either 'noneOverlap' or 'AbsNoneOverlap'.
#' @param top.rank.uv.genes Numeric.Indicates the percentage of top ranked genes that are highly affected
#' by the unwanted variation variables. This is required to be specified when the 'ncg.selection.method' is
#' either 'noneOverlap' or 'AbsNoneOverlap'.
#' @param bio.variables Symbols.Indicates the columns names that contains biological variables in the SummarizedExperiment object.
#' @param bio.groups TBBB
#' @param bio.clustering.method Symbols.Indicates which clustering methods should be used to group continuous sources of biological variation.
#' @param nb.bio.clusters Numeric.Indicates the number of clusters for each continuous  sources of biological variation, by default it is set to 2.
#' @param uv.variables Symbols.Indicates the columns names that contains UV variables in the SummarizedExperiment object.
#' @param uv.groups TBBB
#' @param uv.clustering.method Symbols.Indicates which clustering method should be used to group continuous sources of UV.
#' @param nb.uv.clusters Numeric.Indicates the number of clusters for each continuous sources of UV, by default it is set to 2.
#' @param normalization TBBB
#' @param regress.out.uv.variables TBBB
#' @param regress.out.bio.variables TBBB
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation.
#' @param min.sample.for.aov TBBB
#' @param min.sample.for.correlation TBBB
#' @param corr.method TBBB
#' @param a The significance level used for the confidence intervals in the correlation, by default it is set to 0.05.
#' @param rho The value of the hypothesised correlation to be used in the hypothesis testing, by default it is set to 0.
#' @param anova.method TBBB
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected NCG or not.
#' This analysis involves principal component analysis on the selected NCG and
#' then explore the R^2 or vector correlation between the 'nb.pcs' first principal
#' components and with biological and unwanted variables.
#' @param variables.to.assess.ncg Symbols. Indicates the column names of the SummarizedExperiment object that
#' contain variables whose association with the selected genes as NCG.
#' needs to be evaluated. The default is NULL. This means all the variables in the 'bio.variables' and 'uv.variables' will be assessed.
#' @param nb.pcs Numeric. Indicates the number of the first principal components on selected NCG to be used to assess the performance of NCGs.
#' @param center Logical. Indicates whether to center the data before applying principal component analysis or not. The default is TRUE.
#' @param scale Logical. Indicates whether to scale the data before applying Principal component analysis. The default is FALSE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object or not.
#' @param assess.variables TBBB
#' @param remove.na Symbol. Indicates whether to remove NA or missing values from either the 'assays', the 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be excluded. If 'sample.annotation' is selected, the
#' samples that contains NA or missing values for any 'bio.variables' and 'uv.variables' will be excluded. By default, it is set to both'.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return Either the SummarizedExperiment object containing the a set of negative control genes
#' or a logical vector of the selected negative control genes.

#' @importFrom BiocSingular runSVD bsparam
#' @importFrom fastDummies dummy_cols
#' @importFrom dplyr mutate progress_estimated
#' @importFrom tidyr pivot_longer
#' @importFrom SummarizedExperiment assay
#' @importFrom knitr kable
#' @import ggplot2
#' @export

supervisedFindNcgPbPbio <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        nb.ncg = 10,
        ncg.selection.method = 'Prod',
        grid.nb = 0.05,
        top.rank.bio.genes = 30,
        top.rank.uv.genes = 40,
        bio.groups = NULL,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 2,
        uv.groups = NULL,
        uv.clustering.method = 'kmeans',
        nb.uv.clusters = 2,
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
        scale = FALSE,
        center = TRUE,
        assess.se.obj = TRUE,
        assess.variables = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        verbose = TRUE
){
    printColoredMessage(
        message = '------------The supervisedFindNcgPbPbio function starts:',
        color = 'white',
        verbose = verbose)
    # check some functions inputs ####
    if(length(assay.name) > 1){
        stop('Please provide a single assay name.')
    } else if(nb.ncg > 100 | nb.ncg <= 0){
        stop('The nb.ncg should be a positve value: 0 < nb.ncg < 100.')
    } else if(min.sample.for.aov <= 1){
        stop('The min.sample.for.aov should be at least 2.')
    } else if(min.sample.for.aov >= ncol(se.obj)){
        stop('The min.sample.for.aov should be less than the total number of samples in the data.')
    } else if(min.sample.for.correlation <= 2){
        stop('The min.sample.for.correlation should be at least 3.')
    } else if(min.sample.for.correlation >= ncol(se.obj)){
        stop('The min.sample.for.correlation should be less than the total number of samples in the data.')
    }
    # else if(!ncg.selection.method %in% c('Prod', 'Sum', 'Average', 'Intersect')){
    #     stop('The ncg.selection.method should be one of Prod, Sum or Average')
    # }
    # check the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = unique(c(bio.variables, uv.variables, bio.groups, uv.groups)),
            remove.na = remove.na,
            verbose = verbose)
    }
    # check the correlations between bio variables and uv variable separately ####
    if (assess.variables) {
        se.obj <- variablesCorrelation(
            se.obj = se.obj,
            bio.variables = unique(c(bio.variables,bio.groups)),
            uv.variables = unique(c(uv.variables,uv.groups)),
            cat.cor.coef = c(0.95, 0.95),
            cont.cor.coef = c(0.95, 0.95),
            assess.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose)
        bio.variables <- se.obj$bio.variables
        uv.variables <- se.obj$uv.variables
        se.obj <- se.obj$se.obj
    }
    # data transformation and normalization ####
    ## log transformation ####
    # data transformation and normalization ####
    printColoredMessage(
        message = '-- Data transformation and normalization',
        color = 'magenta',
        verbose = verbose)
    ## apply log ####
    if(apply.log){
        printColoredMessage(
            message = paste0(
                'Applying log2 + ',
                pseudo.count,
                '(pseudo.count) on the ',
                assay.name,
                ' data.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- log2(assay(se.obj, assay.name) + pseudo.count)
    } else {
        printColoredMessage(
            message = paste0('The ', assay.name, ' data will be used without any log transformation.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- assay(se.obj, assay.name)
    }
    ## regress out uv variables ####
    if(!is.null(normalization)){
        printColoredMessage(
            message = paste0('Applying the ',
                             normalization,
                             ' normalization on the ',
                             assay.name,
                             ' data.'),
            color = 'blue',
            verbose = verbose)
        expr.data.nor <- applyOtherNormalizations(
            se.obj = se.obj,
            assay.name = assay.name,
            method = normalization,
            pseudo.count = pseudo.count,
            apply.log = apply.log,
            assess.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose)
    }
    ### regress out unwanted variation variables
    if(!is.null(regress.out.uv.variables) & !is.null(normalization)){
        printColoredMessage(
            message = paste0(
                paste0(regress.out.uv.variables, collapse = ' & '),
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'We do not recommend regressing out ',
                paste0(regress.out.uv.variables, collapse = ' & '),
                ' if they are largely associated with the ',
                paste0(bio.variables, collapse = ' & '),
                '.'),
            color = 'red',
            verbose = verbose)
        expr.data.reg.uv <- t(expr.data.nor)
        uv.variables.all <- paste('se.obj', regress.out.uv.variables, sep = '$')
        adjusted.data <- lm(as.formula(paste(
            'expr.data.reg.uv',
            paste0(uv.variables.all, collapse = '+') ,
            sep = '~')))
        expr.data.reg.uv <- t(adjusted.data$residuals)
        colnames(expr.data.reg.uv) <- colnames(se.obj)
        row.names(expr.data.reg.uv) <- row.names(se.obj)
        rm(adjusted.data)
    }
    if(!is.null(regress.out.uv.variables) & is.null(normalization)){
        printColoredMessage(
            message = paste0(
                paste0(regress.out.uv.variables, collapse = ' & '),
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'We do not recommend regressing out ',
                paste0(regress.out.uv.variables, collapse = ' & '),
                ' if they are largely associated with the',
                paste0(bio.variables, collapse = ' & '),
                '.'),
            color = 'red',
            verbose = verbose)
        expr.data.reg.uv <- t(expr.data)
        uv.variables.all <- paste('se.obj', regress.out.uv.variables, sep = '$')
        adjusted.data <- lm(as.formula(paste(
            'expr.data.reg.uv',
            paste0(uv.variables.all, collapse = '+') ,
            sep = '~'
        )))
        expr.data.reg.uv <- t(adjusted.data$residuals)
        colnames(expr.data.reg.uv) <- colnames(se.obj)
        row.names(expr.data.reg.uv) <- row.names(se.obj)
        rm(adjusted.data)
    }
    ### regress out biological variables
    if(!is.null(regress.out.bio.variables)){
        printColoredMessage(
            message = paste0(
                paste0(regress.out.bio.variables, collapse = ' & '),
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'We do not recommend regressing out ',
                paste0(regress.out.bio.variables, collapse = ' & '),
                ' if they are largely associated with the',
                paste0(uv.variables, collapse = ' & '),
                '.'),
            color = 'red',
            verbose = verbose)
        expr.data.reg.bio <- t(expr.data)
        bio.variables.all <- paste('se.obj', regress.out.bio.variables, sep = '$')
        adjusted.data <- lm(as.formula(paste(
            'expr.data.reg.bio',
            paste0(bio.variables.all, collapse = '+') ,
            sep = '~')))
        expr.data.reg.bio <- t(adjusted.data$residuals)
        colnames(expr.data.reg.bio) <- colnames(se.obj)
        row.names(expr.data.reg.bio) <- row.names(se.obj)
        rm(adjusted.data)
    }
    # step 1: genes that are highly affected by different sources of unwanted variation ####
    ## step1.1 create all possible homogeneous biological groups ####
    printColoredMessage(
        message = '-- Finding genes that are highly affected by different sources of unwanted variation:',
        color = 'magenta',
        verbose = verbose)
    printColoredMessage(
        message = 'This step with be performed within each possible homogeneous biological groups.',
        color = 'red',
        verbose = verbose)
    printColoredMessage(
        message = 'Creating all possible major homogeneous biological groups:',
        color = 'blue',
        verbose = verbose)
    if(is.null(bio.groups) & length(bio.variables) == 1){
        printColoredMessage(
            message = paste0(
                'The ',
                bio.variables,
                ' variable will be used as a major homogeneous biological group.'),
            color = 'blue',
            verbose = verbose)
        if(class(se.obj[[bio.variables]]) %in% c('numeric', 'integer')){
            all.bio.groups <- createHomogeneousBioGroups(
                se.obj = se.obj,
                bio.variables = bio.variables,
                nb.clusters = nb.bio.clusters,
                clustering.method = bio.clustering.method,
                assess.se.obj = FALSE,
                assess.variables = FALSE,
                save.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose)
        } else {
            all.bio.groups <- se.obj[[bio.variables]]
            print(kable(table(se.obj[[bio.variables]])))
        }
    } else if(!is.null(bio.groups) & length(bio.groups) == 1){
        printColoredMessage(
            message = paste0(
                'The ',
                bio.groups,
                ' variable will be used as a major homogeneous biological group.'),
            color = 'blue',
            verbose = verbose)
        if(class(se.obj[[bio.groups]]) %in% c('numeric', 'integer')){
            all.bio.groups <- createHomogeneousBioGroups(
                se.obj = se.obj,
                bio.variables = bio.groups,
                nb.clusters = nb.bio.clusters,
                clustering.method = bio.clustering.method,
                assess.se.obj = FALSE,
                assess.variables = FALSE,
                save.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose)
        } else {
            all.bio.groups <- se.obj[[bio.groups]]
            print(kable(table(se.obj[[bio.groups]])))
        }
    } else if(!is.null(bio.groups) & length(bio.groups) > 1){
        printColoredMessage(
            message = paste0(
                'The',
                paste0(bio.groups, collapse = ' & '),
                ' variables will be used to create all possible major homogeneous biological groups.'),
            color = 'blue',
            verbose = verbose)
        all.bio.groups <- createHomogeneousBioGroups(
            se.obj = se.obj,
            bio.variables = bio.groups,
            nb.clusters = nb.bio.clusters,
            clustering.method = bio.clustering.method,
            assess.se.obj = FALSE,
            assess.variables = FALSE,
            save.se.obj = FALSE,
            verbose = verbose)
    } else if(is.null(bio.groups) & length(bio.variables) > 1){
        printColoredMessage(
            message = paste0(
                'The ',
                paste0(bio.variables, collapse = ' & '),
                ' variables will be used to create all possible major homogeneous biological groups.'),
            color = 'blue',
            verbose = verbose)
        all.bio.groups <- createHomogeneousBioGroups(
            se.obj = se.obj,
            bio.variables = bio.variables,
            nb.clusters = nb.bio.clusters,
            clustering.method = bio.clustering.method,
            save.se.obj = FALSE,
            assess.se.obj = FALSE,
            assess.variables = FALSE,
            verbose = verbose)
    }
    ## step 1.2 correlation between gene expression and all continuous source of unwanted variation with each biological groups ####
    uv.var.class <- unlist(lapply(
        uv.variables,
        function(x) class(colData(se.obj)[[x]])) )
    continuous.uv <- uv.variables[uv.var.class %in% c('numeric', 'integer')]
    if(length(continuous.uv) > 0){
        printColoredMessage(
            message = '-- Correlation analyses:',
            color = 'magenta',
            verbose = verbose)
        if(length(continuous.uv)==1){
            group = 'is'
            sour = 'source'
        } else {
            group = 'are'
            sour = 'sources'
        }
        printColoredMessage(
            message = paste0(
                length(continuous.uv),
                ' continuous ',
                sour,
                ' of unwanted variation:',
                paste0(continuous.uv, collapse = ' & '),
                ' ',
                group,
                ' found in the uv.variables.'),
            color = 'blue',
            verbose = verbose)
        selected.bio.groups <- findRepeatingPatterns(
            vec = all.bio.groups,
            n.repeat = min.sample.for.correlation
            )
        if(length(selected.bio.groups) > 0){
            if(length(selected.bio.groups) == 1){
                group = 'group has'
            } else group = 'groups have'
            printColoredMessage(
                message = paste0(
                    length(selected.bio.groups),
                    ' homogeneous biological ',
                    group,
                    ' at least ',
                    min.sample.for.correlation,
                    ' (min.sample.for.correlation) samples to pefrom correlation analysis between gene-level ',
                    'expression and all the continuous sources of unwanted variation.'),
                color = 'blue',
                verbose = verbose
            )
            if(is.null(regress.out.bio.variables)){
                data.to.use <- expr.data
            } else data.to.use <- expr.data.reg.bio
            corr.genes.uv <- lapply(
                continuous.uv,
                function(x) {
                    all.corr <- lapply(
                        selected.bio.groups,
                        function(y){
                            selected.samples <- all.bio.groups == y
                            corr.genes <- as.data.frame(correls(
                                y = se.obj@colData[, x][selected.samples],
                                x = t(data.to.use[ , selected.samples]),
                                type = corr.method,
                                a = a ,
                                rho = rho))
                            corr.genes <- cbind(
                                round(x = corr.genes[, 1:4], digits = 3),
                                corr.genes[, 5, drop = FALSE])
                            set.seed(2233)
                            corr.genes$ranked.genes <- rank(-abs(corr.genes[, 'correlation']), ties.method = 'random')
                            row.names(corr.genes) <- row.names(data.to.use)
                            corr.genes
                        })
                    names(all.corr) <- selected.bio.groups
                    all.corr
                })
            names(corr.genes.uv) <- continuous.uv
        } else {
            stop(paste0(
                'There are not homogeneous biological groups that have at least ',
                min.sample.for.correlation,
                ' (min.sample.for.correlation) samples for correlation analysis',
                ' between gene-level expression and all the continuous sources of unwanted variation.'))
            }
    } else corr.genes.uv <- NULL
    ## step1.3 anova between gene expression and all categorical source of variation with each biological groups ####
    categorical.uv <- uv.variables[uv.var.class %in% c('factor', 'character')]
    if(length(categorical.uv) > 0) {
        printColoredMessage(
            message = '-- ANOV analyses:',
            color = 'magenta',
            verbose = verbose)
        anova.genes.uv <- lapply(
            categorical.uv,
            function(x) {
                bio.batch <- table(all.bio.groups, colData(se.obj)[[x]])
                l.gorup <- length(unique(se.obj[[x]]))
                if (sum(rowSums(bio.batch >= min.sample.for.aov) == l.gorup ) > 0) {
                    if(sum(rowSums(bio.batch > min.sample.for.aov) == l.gorup) == 1){
                        group = 'group has'
                    } else group = 'groups have'
                    printColoredMessage(
                        message = paste0(
                            sum(rowSums(bio.batch > min.sample.for.aov) == l.gorup),
                            ' homogeneous biological ',
                            group,
                            ' at least ',
                            min.sample.for.aov,
                            ' (min.sample.for.aov) samples within individual batches of the ',
                            x,
                            ' variable.'),
                        color = 'blue',
                        verbose = verbose)
                } else {
                    printColoredMessage(
                        message = paste0(
                            'There are not homogeneous biological groups have at least ',
                            min.sample.for.aov ,
                            ' (min.sample.for.aov) samples within each batches of the ',
                            x,
                            ' variable. Checking the connection between batches:'),
                        color = 'red',
                        verbose = verbose)
                    connection.check <- lapply(
                        1:nrow(bio.batch),
                        function(y) {
                            batch.names.a <- names(which(bio.batch[y, ] >= min.sample.for.aov))
                            con.prps <- lapply(
                                    c(1:nrow(bio.batch))[-y],
                                    function(z) {
                                        batch.names.b <- names(which(bio.batch[z, ] >= min.sample.for.aov))
                                        inter.samples <- intersect(batch.names.a, batch.names.b)
                                        if (length(inter.samples) > 0) {
                                            sort(unique(c(batch.names.a, batch.names.b)), decreasing = FALSE)
                                        } else {
                                            sort(batch.names.a, decreasing = FALSE)
                                        }
                                    })
                            sort(unique(unlist(Filter(Negate(is.null), con.prps))), decreasing = FALSE)
                        })
                    covered.batches <- Filter(Negate(is.null), connection.check)
                    covered.batches <- unlist(lapply(covered.batches, length))
                    if (max(covered.batches) == length(unique(se.obj[[x]]))) {
                        printColoredMessage(
                            message = paste0(
                                'The connections between the homogeneous bioloical groups can cover all the batches of the ',
                                x,
                                ' variable.'),
                            color = 'blue',
                            verbose = verbose)
                    } else {
                        printColoredMessage(
                            message = paste0(
                                'Please note that the homogeneous bioloical groups do not cover all the batches of the ',
                                x,
                                ' variable.'),
                            color = 'blue',
                            verbose = verbose)
                    }
                }
                selected.bio.groups <- names(which(rowSums(bio.batch >= min.sample.for.aov) > 1))
                if(is.null(regress.out.bio.variables)){
                    data.to.use <- expr.data
                } else data.to.use <- expr.data.reg.bio
                all.anova <- lapply(
                    selected.bio.groups,
                    function(i){
                        selected.samples <- all.bio.groups == i
                        if(anova.method == 'aov'){
                            anova.genes.batch <- as.data.frame(
                                row_oneway_equalvar(
                                    x = data.to.use[ , selected.samples],
                                    g = se.obj@colData[, x][selected.samples]))
                        } else if(anova.method == 'welch.correction'){
                            anova.genes.batch <- as.data.frame(
                                row_oneway_welch(
                                    x = data.to.use[ , selected.samples],
                                    g = se.obj@colData[, x][selected.samples]))
                        }
                        set.seed(2233)
                        anova.genes.batch$ranked.genes <- rank(-anova.genes.batch[, 'statistic'], ties.method = 'random')
                        anova.genes.batch
                    })
                names(all.anova) <- selected.bio.groups
                all.anova
            })
        names(anova.genes.uv) <- categorical.uv
    } else anova.genes.uv <- NULL
    # step 2 genes that are not highly affected by biology ####
    printColoredMessage(
        message = '-- Step2: finding genes that are not highly affected by different sources of biological variation:',
        color = 'magenta',
        verbose = verbose)
    printColoredMessage(
        message = 'This step with be performed within each possible homogeneous unwanted groups.',
        color = 'red',
        verbose = verbose)
    ## step 2.1 create all possible homogeneous uv groups ####
    printColoredMessage(
        message = '### Creating all possible major groups with respect to sources of unwanted variation:',
        color = 'magenta',
        verbose = verbose)
    if(is.null(uv.groups) & length(uv.variables) == 1){
        printColoredMessage(
            message = paste0(
                'The ',
                uv.variables,
                ' will be used as a major source of unwanted variation',
                ' to create all possible groups.'),
            color = 'blue',
            verbose = verbose)
        if(class(se.obj[[uv.variables]]) %in% c('numeric', 'integer')){
            all.uv.groups <- createHomogeneousUVGroups(
                se.obj = se.obj,
                uv.variables = uv.variables,
                nb.clusters = nb.uv.clusters,
                clustering.method = uv.clustering.method,
                assess.se.obj = FALSE,
                assess.variables = FALSE,
                save.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose)
        } else {
            all.uv.groups <- se.obj[[uv.variables]]
            print(kable(table(se.obj[[uv.variables]])))
            }
    } else if(!is.null(uv.groups) & length(uv.groups) == 1){
        printColoredMessage(
            message = paste0(
                'The',
                uv.groups,
                ' variables will be used as a major source of unwanted variation',
                ' to create all possible groups.'),
            color = 'blue',
            verbose = verbose)
        if(class(se.obj[[uv.groups]]) %in% c('numeric', 'integer')){
            all.uv.groups <- createHomogeneousUVGroups(
                se.obj = se.obj,
                uv.variables = uv.groups,
                nb.clusters = nb.uv.clusters,
                clustering.method = uv.clustering.method,
                assess.se.obj = FALSE,
                assess.variables = FALSE,
                remove.na = 'none',
                verbose = verbose)
        } else{
            all.uv.groups <- se.obj[[uv.groups]]
            print(kable(table(se.obj[[uv.groups]])))
        }
    } else if(!is.null(uv.groups) & length(uv.groups) > 1){
        printColoredMessage(
            message = paste0(
                'The ',
                paste0(uv.variables, collapse = ' & '),
                ' variables will be used as a major sources of unwanted variation',
                ' to find all possible groups.'),
            color = 'blue',
            verbose = verbose)
        all.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = uv.groups,
            nb.clusters = nb.uv.clusters,
            clustering.method = uv.clustering.method,
            assess.se.obj = FALSE,
            assess.variables = FALSE,
            verbose = verbose)
    } else if(is.null(uv.groups) & length(uv.variables) > 1){
        printColoredMessage(
            message = paste0(
                'The ',
                paste0(uv.variables, collapse = ' & '),
                ' variables will be used as a major sources of unwanted variation',
                ' to find all possible groups.'),
            color = 'blue',
            verbose = verbose)
        all.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = uv.variables,
            nb.clusters = nb.uv.clusters,
            clustering.method = uv.clustering.method,
            save.se.obj = FALSE,
            assess.se.obj = FALSE,
            assess.variables = FALSE,
            verbose = verbose)
    }
    ## step 2.2 correlation between gene expression and all continuous source of biological variation with each uv groups ####
    bio.var.class <- unlist(lapply(
        bio.variables,
        function(x) class(colData(se.obj)[[x]])))
    continuous.bio <- bio.variables[bio.var.class %in% c('numeric', 'integer')]
    if(length(continuous.bio) > 0){
        printColoredMessage(
            message = '-- Correlation analyses:',
            color = 'magenta',
            verbose = verbose)
        if(length(continuous.bio)==1){
            group = 'is'
            sour = 'source'
        }else{
            group = 'are'
            sour = 'sources'
        }
        printColoredMessage(
            message = paste0(
                length(continuous.bio),
                ' continuous ',
                sour,
                ' of biological variation: ',
                paste0(continuous.bio, collapse = ' & '),
                ' ',
                group,
                ' found in the bio.variables.'),
            color = 'blue',
            verbose = verbose)
        selected.uv.groups <- findRepeatingPatterns(
            vec = all.uv.groups,
            n.repeat = min.sample.for.correlation)
        if(length(selected.uv.groups) > 0){
            if(length(selected.uv.groups) == 1){
                group = 'group has'
            }else group = 'groups have'
            printColoredMessage(
                message = paste0(
                    length(selected.uv.groups),
                    ' homogeneous groups with respect to the sources of unwanted variation ',
                    group,
                    ' at least ',
                    min.sample.for.correlation,
                    ' (min.sample.for.correlation) samples to pefrom correlation between gene-level',
                    'expression and all the continuous sources of bioloical variation.'),
                color = 'red',
                verbose = verbose)
            if(is.null(regress.out.uv.variables) & is.null(normalization)){
                data.to.use <- expr.data
            } else if(!is.null(regress.out.uv.variables) & !is.null(normalization)){
                data.to.use <- expr.data.reg.uv
            } else if(is.null(regress.out.uv.variables) & !is.null(normalization)){
                data.to.use <- expr.data.nor
            }
            corr.genes.bio <- lapply(
                continuous.bio,
                function(x) {
                    all.corr <- lapply(
                        selected.uv.groups,
                        function(y){
                            selected.samples <- all.uv.groups == y
                            corr.genes <- as.data.frame(correls(
                                y = se.obj@colData[, x][selected.samples],
                                x = t(data.to.use[ , selected.samples]),
                                type = corr.method,
                                a = a ,
                                rho = rho))
                            corr.genes <- cbind(
                                round(x = corr.genes[, 1:4], digits = 3),
                                corr.genes[, 5, drop = FALSE]
                            )
                            set.seed(2233)
                            corr.genes$ranked.genes <- rank(abs(corr.genes[, 'correlation']), ties.method = 'random')
                            row.names(corr.genes) <- row.names(data.to.use)
                            corr.genes
                        })
                    names(all.corr) <- selected.uv.groups
                    all.corr
                })
            names(corr.genes.bio) <- continuous.bio
        } else {
            stop(paste0(
                'There are not homogeneous groups with respect to sources of unwanted variation that have at least ',
                min.sample.for.correlation,
                ' (min.sample.for.correlation) samples for correlation analysis between gene-level expression and all',
                ' the continuous sources of bioloical variation.'))
        }
    } else corr.genes.bio <- NULL
    ## step 2.3 anova between gene expression and all categorical source of biological variation with each uv groups ####
    categorical.bio <- bio.variables[bio.var.class %in% c('factor', 'character')]
    if(length(categorical.bio) > 0) {
        printColoredMessage(
            message = '-- ANOV analyses:',
            color = 'magenta',
            verbose = verbose)
        anova.genes.bio <- lapply(
            categorical.bio,
            function(x) {
                bio.batch <- table(all.uv.groups, colData(se.obj)[[x]])
                l.gorup <- length(unique(se.obj[[x]]))
                if (sum(rowSums(bio.batch >= min.sample.for.aov) == l.gorup ) > 0) {
                    if(sum(rowSums(bio.batch > min.sample.for.aov) == l.gorup) == 1){
                        group = 'group has'
                    } else group = 'groups have'
                    printColoredMessage(
                        message = paste0(
                            sum(rowSums(bio.batch > min.sample.for.aov) == l.gorup),
                            ' homogeneous groups with respect to sources of unwanted variation ',
                            group,
                            'at least ',
                            min.sample.for.aov,
                            ' samples within individual groups of the ',
                            x,
                            ' variable.'),
                        color = 'blue',
                        verbose = verbose)
                } else {
                    printColoredMessage(
                        message = paste0(
                            'No homogeneous groups with respect to sources of unwanted variation have at least ',
                            min.sample.for.aov ,
                            ' samples within each group of the ',
                            x,
                            ' variable. Checking the connection between batches:'),
                        color = 'blue',
                        verbose = verbose
                    )
                    batch.connection <- lapply(
                        1:nrow(bio.batch),
                        function(y) {
                            batch.names.a <- names(which(bio.batch[y, ] >= min.sample.for.aov))
                            if(length(batch.names.a) > 0){
                                con.prps <- lapply(
                                    c(1:nrow(bio.batch))[-y],
                                    function(z) {
                                        batch.names.b <- names(which(bio.batch[z, ] >= min.sample.for.aov))
                                        inter.samples <- intersect(batch.names.a, batch.names.b)
                                        if (length(inter.samples) > 0) {
                                            sort(unique(c(batch.names.a, batch.names.b)), decreasing = FALSE)
                                        } else {
                                            sort(batch.names.a, decreasing = FALSE)
                                        }
                                    })
                                sort(unique(unlist(Filter( Negate(is.null), con.prps))), decreasing = FALSE)
                            }
                        })
                    covered.batches <- Filter(Negate(is.null), batch.connection)
                    covered.batches <- unlist(lapply(covered.batches, length))
                    if (max(covered.batches)== length(unique(se.obj[[x]]))) {
                        printColoredMessage(
                            message = paste0(
                                'The connections between the homogeneous groups with respect to unwanted',
                                'variation can cover all the groups of the ',
                                x,
                                ' variable.'),
                            color = 'white',
                            verbose = verbose)
                    } else {
                        printColoredMessage(
                            message = paste0(
                                'Please note that the homogeneous groups with respect to unwanted variation do not cover all the groups of the ',
                                x,
                                ' variable.'),
                            color = 'blue',
                            verbose = verbose)
                    }
                }
                selected.uv.groups <- names(which(rowSums(bio.batch >= min.sample.for.aov) > 1))
                if(length(selected.uv.groups) == 0){
                    stop(paste0(
                        'It seems there is complete association between ',
                        x,
                        ' homogeneous groups with respect to unwanted variation.'))
                }
                if(is.null(regress.out.uv.variables) & is.null(normalization)){
                    data.to.use <- expr.data
                } else if(!is.null(regress.out.uv.variables) & !is.null(normalization)){
                    data.to.use <- expr.data.reg.uv
                } else if(is.null(regress.out.uv.variables) & !is.null(normalization)){
                    data.to.use <- expr.data.nor
                }
                all.anova <- lapply(
                    selected.uv.groups,
                    function(i){
                        selected.samples <- all.uv.groups == i
                        anova.gene.bio <- as.data.frame(
                            row_oneway_equalvar(
                                x = data.to.use[ , selected.samples],
                                g = se.obj@colData[, x][selected.samples]))
                        set.seed(2233)
                        anova.gene.bio$ranked.genes <- rank(anova.gene.bio[, 'statistic'], ties.method = 'random')
                        anova.gene.bio
                    })
                names(all.anova) <- selected.uv.groups
                all.anova
            })
        names(anova.genes.bio) <- categorical.bio
    } else anova.genes.bio <- NULL
    # step3: final selection ####
    printColoredMessage(
        message = '-- Selection of as set of genes as NCG:',
        color = 'magenta',
        verbose = verbose)
    ## step 3.1 product of ranks ####
    if(ncg.selection.method %in% c('Prod', 'Sum', 'Average')){
        all.tests <- c('anova.genes.bio', 'corr.genes.bio', 'anova.genes.uv','corr.genes.uv')
        ncg.selected <- lapply(
            all.tests,
            function(x) {
                if (!is.null(x)) {
                    temp.data <- get(x)
                    ranks.data <- lapply(
                        names(temp.data),
                        function(y) {
                            all.ranks <- lapply(
                                names(temp.data[[y]]),
                                function(i) temp.data[[y]][[i]]$ranked.genes)
                            names(all.ranks) <- names(temp.data[[y]])
                            all.ranks <- do.call(cbind, all.ranks)
                            all.ranks
                        })
                    ranks.data <- do.call(cbind, ranks.data)
                    names(ranks.data) <- names(temp.data)
                    ranks.data
                }
            })
        ncg.selected <- do.call(cbind, ncg.selected)
        if(ncg.selection.method == 'Prod'){
            printColoredMessage(
                message = 'A set of NCG will be selected based on the product of ranks.',
                color = 'blue',
                verbose = verbose)
            stat.rank <- 10^(rowSums(log(ncg.selected, base = 10)))
            if(sum(is.infinite(stat.rank)) > 0){
                stop('The product of ranks results in infinity values.')
            }
        } else if(ncg.selection.method == 'Sum'){
            printColoredMessage(
                message = 'A set of NCG will be selected based on the sum of ranks.',
                color = 'blue',
                verbose = verbose)
            stat.rank <- rowSums(x = ncg.selected, na.rm = TRUE)
        } else if(ncg.selection.method == 'Average'){
            printColoredMessage(
                message = 'A set of NCG will be selected based on the average of ranks.',
                color = 'blue',
                verbose = verbose)
            stat.rank <- rowMeans(x = ncg.selected, na.rm = TRUE)
        }
        ncg.selected <- as.data.frame(ncg.selected)
        row.names(ncg.selected) <- row.names(se.obj)
        ncg.selected$stat.rank <- stat.rank
        ncg.selected$final.rank <- rank(ncg.selected$stat.rank)
        ncg.selected <- ncg.selected[order(ncg.selected$stat.rank, decreasing = FALSE) , ]
        ncg.selected <- row.names(ncg.selected[1:round(nb.ncg/100 * nrow(se.obj), digits = 0) , ])
        ncg.selected <- row.names(se.obj) %in% ncg.selected
    } else if(ncg.selection.method == 'noneOverlap'){
        printColoredMessage(
            message = 'A set of NCG is selected based on the noneOverlap approach.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'The non-overlap set of genes between top ',
                top.rank.bio.genes * 100,
                '% of highly affected genes by the bioloigcal variation and top ',
                top.rank.uv.genes * 100,
                '% of highly affected genes by the unwanted variation.'),
            color = 'blue',
            verbose = verbose)
        if(top.rank.bio.genes == 100){
            top.rank.bio.genes.nb <- top.rank.bio.genes
        } else{
            top.rank.bio.genes.nb <- round(top.rank.bio.genes/100 * nrow(se.obj), digits = 0)
            top.rank.bio.genes.nb <- c(nrow(se.obj) - top.rank.bio.genes.nb)
        }
        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
        top.bio.genes <- unique(unlist(lapply(
            all.bio.tests,
            function(x){
                if(!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            all.ranks <- sapply(
                                names(temp.data[[y]]),
                                function(z) temp.data[[y]][[z]]$ranked.genes)
                            set.seed(2233)
                            all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                            row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                        })))
                }
            })))
        top.rank.uv.genes <- round(top.rank.uv.genes/100 * nrow(se.obj), digits = 0)
        top.uv.genes <- c()
        ncg.selected <- c()
        lo <- nrow(se.obj) - top.rank.uv.genes
        grid.nb <- round(c(grid.nb/100) * nrow(se.obj), digits = 0)
        pro.bar <- progress_estimated(round(grid.nb, digits = 0) + 2)
        while(length(ncg.selected) < round(nb.ncg* nrow(se.obj), digits = 0)){
            pro.bar$pause(0.1)$tick()$print()
            all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
            top.uv.genes <- unique(unlist(lapply(
                all.uv.tests,
                function(x){
                    temp <- get(x)
                    if(length(names(temp))!=0){
                        ranks.data <- lapply(
                            names(temp),
                            function(y){
                                unlist(lapply(
                                    names(temp[[y]]),
                                    function(z){
                                        index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes
                                        row.names(temp[[y]][[z]])[index]
                                    }))
                            })
                    }
                })))
            top.rank.uv.genes <- top.rank.uv.genes + grid.nb
            if(top.rank.uv.genes > nrow(se.obj)){
                message(' ')
                printColoredMessage(message = paste0(length(ncg.selected), ' genes are found based on the current parameters.'),
                                    color = 'red', verbose = verbose)
                stop('The requested number of genes cannot be found. Please lower either the values of top.rank.bio.genes or nb.ncg.')
            }
            ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
            ncg.selected
        }
        ncg.selected <- row.names(se.obj) %in% ncg.selected
        if(top.rank.uv.genes >= 100){
            top.rank.uv.genes = 100
        }
        message(' ')
        printColoredMessage(
            message = paste0(
                'The non-overlap set of genes between top ',
                top.rank.bio.genes * 100,
                '% of highly affected genes by the bioloigcal variation and top ',
                top.rank.uv.genes,
                '% of highly affected genes by the unwanted variation.'),
            color = 'blue',
            verbose = verbose)
    } else if(ncg.selection.method == 'AbsNoneOverlap'){
        printColoredMessage(
            message = 'A set of NCG is selected based on the AbsNoneOverlap approach.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'The non-overlap set of genes between top ',
                top.rank.bio.genes * 100,
                '% of highly affected genes by the bioloigcal variation and top ',
                top.rank.uv.genes * 100,
                '% of highly affected genes by the unwanted variation.'),
            color = 'blue',
            verbose = verbose)
        if(top.rank.bio.genes == 100){
            top.rank.bio.genes.nb <- top.rank.bio.genes
        } else{
            top.rank.bio.genes <- round(top.rank.bio.genes/100 * nrow(se.obj), digits = 0)
            top.rank.bio.genes.nb <- c(nrow(se.obj) - top.rank.bio.genes)
        }
        all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
        top.bio.genes <- unique(unlist(lapply(
            all.bio.tests,
            function(x){
                if(!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            all.ranks <- sapply(
                                names(temp.data[[y]]),
                                function(z) temp.data[[y]][[z]]$ranked.genes)
                            set.seed(2233)
                            all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                            row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]

                        })))
                }
            })))
        top.rank.uv.genes <- round(top.rank.uv.genes/100 * nrow(se.obj), digits = 0)
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        top.uv.genes <- unique(unlist(lapply(
            all.uv.tests,
            function(x){
                if(!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            all.ranks <- sapply(
                                names(temp.data[[y]]),
                                function(z) temp.data[[y]][[z]]$ranked.genes)
                            set.seed(2233)
                            all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                            row.names(se.obj)[all.ranks < top.rank.uv.genes]

                        })))
                }
            })))
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        ncg.selected <- row.names(se.obj) %in% ncg.selected
    }
    printColoredMessage(
        message = paste0('A set of ', sum(ncg.selected), ' genes are selected for NCG.'),
        color = 'blue',
        verbose = verbose)
        # step 4: assessment of selected set of NCG  ####
        # pca
    if(assess.ncg){
        printColoredMessage(
            message = '-- Assess the performance of selected NCG set:',
            color = 'magenta',
            verbose = verbose)
        if(is.null(variables.to.assess.ncg)){
            all.variables <- c(bio.variables, uv.variables)
        } else all.variables <- variables.to.assess.ncg
        printColoredMessage(
            message = 'Peform PCA on only selected genes as NCG.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'Then, explore the association of the first ',
                nb.pcs,
                '  with the ',
                paste0(c(all.variables, uv.variables), collapse = ' & '),
                ' variables.'),
            color = 'blue',
            verbose = verbose)
        if(apply.log){
            temp.data <- log2(assay(se.obj, assay.name) + pseudo.count)
        } else{
            temp.data <- assay(se.obj, assay.name)
        }
        pca.data <- BiocSingular::runSVD(
            x = t(temp.data[ncg.selected, ]),
            k = nb.pcs,
            BSPARAM =  bsparam(),
            center = center,
            scale = scale)$u
        all.corr <- lapply(
            all.variables,
            function(x){
                if(class(se.obj[[x]]) %in% c('numeric', 'integer')){
                    rSquared <- sapply(
                        1:nb.pcs,
                        function(y) summary(lm(se.obj[[x]] ~ pca.data[, 1:y]))$r.squared)
                } else if(class(se.obj[[x]]) %in% c('factor', 'character')){
                    catvar.dummies <- dummy_cols(se.obj[[x]])
                    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                    cca.pcs <- sapply(
                        1:nb.pcs,
                        function(y){ cca <- cancor(
                            x = pca.data[, 1:y, drop = FALSE],
                            y = catvar.dummies)
                        1 - prod(1 - cca$cor^2)
                        })
                }
            })
        pcs <- groups <- NULL
        names(all.corr) <- all.variables
        pca.ncg <- as.data.frame(do.call(cbind, all.corr))
        pca.ncg['pcs'] <- c(1:nb.pcs)
        pca.ncg <- tidyr::pivot_longer(data = pca.ncg,-pcs,
                                       names_to = 'groups',
                                       values_to = 'ls')
        pca.ncg <- ggplot(pca.ncg, aes(x = pcs, y = ls, group = groups)) +
            geom_line(aes(color = groups), size = 1) +
            geom_point(aes(color = groups), size = 2) +
            xlab('PCs') +
            ylab (expression("Correlations")) +
            scale_x_continuous(
                breaks = (1:nb.pcs),
                labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
            scale_y_continuous(
                breaks = scales::pretty_breaks(n = nb.pcs),
                limits = c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 10, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 14),
                strip.text.x = element_text(size = 10),
                plot.title = element_text(size = 16)
            )
        if(verbose) print(pca.ncg)
    }
    # add results to the SummarizedExperiment object ####
    out.put.name <- paste0(
        sum(ncg.selected),
        '|',
        paste0(bio.variables, collapse = '&'),
        '|',
        paste0(uv.variables, collapse = '&'),
        '|PbPbio:',
        ncg.selection.method,
        '|',
        assay.name)
    if(save.se.obj == TRUE){
        printColoredMessage(
            message = '-- Saving a selected set of NCG to the metadata of the SummarizedExperiment object.',
            color = 'magenta',
            verbose = verbose)
        ## Check if metadata NCG already exists
        if(length(se.obj@metadata$NCG) == 0 ) {
            se.obj@metadata[['NCG']] <- list()
        }
        se.obj@metadata[['NCG']][['supervised']][[out.put.name]] <- ncg.selected
        printColoredMessage(
            message = 'The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The supervisedFindNcgPbPbio function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj)
    } else{
        printColoredMessage(
            message = '------------The supervisedFindNcgPbPbio function finished.',
            color = 'white',
            verbose = verbose)
        return(ncg.selected)
    }
}

