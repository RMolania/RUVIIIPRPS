#' is used to find a set of negative control genes (NCG) using ANOVA and correlation analyses across all samples.

#' @author Ramyar Molania

#' @description
#' This function uses the correlation and ANOVA analyses across all samples to find a set of genes as negative control
#' genes for RUV-III-PRPS normalization. The correlation and ANOVA are used to find genes that are highly affected by
#' continuous and categorical sources of variation respectively. The function selects genes as NCG that show possible
#' high correlation coefficients and F-statistics with the sources of unwanted variation and low correlation coefficients
#' and F-statistics with the sources of biological variation. The function uses different approaches to perform the final selection.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. Indicates a name of an assay in the SummarizedExperiment object. The selected assay should
#' be the one that will be used for RUV-III-PRPS normalization.
#' @param nb.ncg Numeric. Indicates how many genes should be selected as NCG. The value is the percentage of the total
#' genes in the SummarizedExperiment object. The default is 10 percent.
#' @param bio.variables Symbols. Indicates the columns names that contain biological variables in the
#' SummarizedExperiment object.
#' @param uv.variables Symbols. Indicates the columns names that contains UV variables in the SummarizedExperiment object.
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
#' @param min.sample.for.aov Numeric.Indicates the minimum number of samples to perform correlation analyses between continuous sources
#' of variation (biological and unwanted variation) with individual gene expression.
#' @param min.sample.for.correlation Numeric.Indicates the minimum number of samples to perform correlation analyses between continuous sources
#' of variation (biological and unwanted variation) with individual gene expression.
#' @param regress.out.uv.variables Symbol.Indicates the column names in the SummarizedExperiment object that will be regressed out from the data
#' before performing correlation and ANOVA. The default is NULL.
#' @param regress.out.bio.variables Symbol.Indicates the column names in the SummarizedExperiment object that will be regressed out from the data
#' before performing correlation and ANOVA. The default is NULL.
#' @param normalization Symbols.Indicates winch normalization method to use for library size adjustment before fining genes that are highly
#' affected by biological variation. The default is CPM. Refer to the 'applyOtherNormalization' function for further details.
#' @param corr.method Numeric.Indicates which correlation methods (pearson or spearman) should be used for the correlation analyses. The
#' default is 'spearman'.
#' @param a The significance level used for the confidence intervals in the correlation, by default it is set to 0.05.
#' @param rho The value of the hypothesised correlation to be used in the hypothesis testing, by default it is set to 0.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation.
#' @param anova.method Symbols. Indicates which ANOVA methods to use. The default in 'aov'. Refer to the function "genesVariableAnova" for more details.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected NCG or not.
#' This analysis involves principal component analysis on the selected NCG and
#' then explore the R^2 or vector correlation between the 'nb.pcs' first principal
#' components and with biological and unwanted variables.
#' @param variables.to.assess.ncg Symbols. Indicates the column names of the SummarizedExperiment object that
#' contain variables whose association with the selected genes as NCG.
#' needs to be evaluated. The default is NULL. This means all the variables in the 'bio.variables' and 'uv.variables' will be assessed.
#' @param nb.pcs Numeric. Indicates the number of the first principal components on selected NCG to be used to assess the performance of NCGs.
#' @param center Logical. Indicates whether to scale the data or not before applying SVD. If center is TRUE, then
#' centering is done by subtracting the column means of the assay from their corresponding columns. The default is TRUE.
#' @param scale Logical. Indicates whether to scale the data or not before applying SVD.  If scale is TRUE, then scaling
#' is done by dividing the (centered) columns of the assays by their standard deviations if center is TRUE, and the root
#' mean square otherwise. The default is FALSE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object or not.
#' @param assess.variables Logical. Indicates whether to assess the association between the biological and unwanted
#' variation variables separately. Refer to the 'assessVariablesCorrelation' function for more details. Th default is FALSE.
#' @param remove.na Symbol. Indicates whether to remove NA or missing values from either the 'assays', the 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be excluded. If 'sample.annotation' is selected, the
#' samples that contains NA or missing values for any 'bio.variables' and 'uv.variables' will be excluded. By default, it is set to both'.
#' @param cat.cor.coef Vector of two numerical values. Indicates the cut-off of the correlation coefficient between each pair of
#' categorical variables. The first one is between each pair of 'uv.variables' and the second one is between each pair of 'bio.variables'.
#' The correlation is computed by the function ContCoef from the DescTools package. If the correlation of a pair of variable is higher than
#' the cut-off, then only the variable that has the highest number of factor will be kept and the other one will be excluded from the
#' remaining analysis. By default they are both set to 0.7.
#' @param cont.cor.coef Vector of two numerical values. Indicates the cut-off of the Spearman correlation coefficient between each pair of
#' continuous variables. The first one is between each pair of 'uv.variables' and the second one is between each pair of 'bio.variables'.
#' If the correlation of a pair of variable is higher than the cut-off, then only the variable that has the highest variance will
#' be kept and the other one will be excluded from the remaining analysis. By default they are both set to 0.7.
#' @param save.se.obj Logical. Indicates whether to save the result of the function in the metadata of the SummarizedExperiment object or
#' to output the result. The default is TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return Either the SummarizedExperiment object containing the a set of negative control genes or a logical vector of
#' the selected negative control genes

#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @importFrom matrixStats rowProds
#' @export

findNCGsAcrossSamples <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        nb.ncg = 10,
        ncg.selection.method = 'AbsNoneOverlap',
        grid.nb = 1,
        top.rank.bio.genes = 50,
        top.rank.uv.genes = 50,
        min.sample.for.aov = 3,
        min.sample.for.correlation = 10,
        regress.out.uv.variables = NULL,
        regress.out.bio.variables = NULL,
        normalization = 'CPM',
        apply.log = TRUE,
        pseudo.count = 1,
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
        cat.cor.coef = c(0.95, 0.95),
        cont.cor.coef = c(0.95, 0.95),
        save.se.obj = TRUE,
        remove.na = 'both',
        verbose = TRUE
){
    printColoredMessage(message = '------------The supervisedFindNcgAnoCorrAs function starts:',
                        color = 'white',
                        verbose = verbose)
    # check several functions inputs ####
    if(length(assay.name) > 1){
        stop('Please provide a single assay name.')
    } else if(nb.ncg > 100 | nb.ncg <= 0){
        stop('The nb.ncg should be a positve value  0 < nb.ncg < 100.')
    } else if (is.null(bio.variables)){
        stop('The bio.variables cannot be empty.')
    } else if (length(intersect(bio.variables, uv.variables)) > 0){
        stop('The variable should be either bio.variables or uv.variables.')
    } else if (!ncg.selection.method %in% c('Prod', 'Sum', 'Average', 'noneOverlap', 'AbsNoneOverlap')){
        stop('The ncg.selection.method should be one of "Prod", "Sum", "Average", "noneOverlap" or "AbsNoneOverlap".')
    } else if (top.rank.bio.genes > 100 | top.rank.bio.genes <= 0){
        stop('The top.rank.bio.genes should be a positve value  0 < top.rank.bio.genes < 100.')
    } else if (top.rank.uv.genes > 100 | top.rank.uv.genes <= 0){
        stop('The top.rank.uv.genes should be a positve value  0 < top.rank.uv.genes < 100.')
    } else if (min.sample.for.aov <= 1){
        stop('The min.sample.for.aov should be at least 2.')
    } else if (is.null(min.sample.for.aov)){
        stop('The min.sample.for.aov cannot be empty.')
    } else if (min.sample.for.correlation >= ncol(se.obj) | min.sample.for.correlation < 3){
        stop('The min.sample.for.correlation should be more than 2 and less than the total number of samples in the data.')
    } else if (is.null(min.sample.for.correlation)){
        stop('The min.sample.for.correlation cannot be empty.')
    } else if (length(intersect(regress.out.bio.variables , regress.out.uv.variables)) > 0){
        stop('The variable to regress out should be either in regress.out.bio.variables or regress.out.uv.variables.')
    }
    if(!is.null(normalization)){
        if(!is.null(regress.out.uv.variables)){
            printColoredMessage(
                message = paste0('Both normalization and regress.out.uv.variables are selected.',
                'The function will perfom normalization first and the regression the UV variables.'),
                color = 'magenta',
                verbose = verbose)
        }
    }
    # check the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(bio.variables, uv.variables),
            remove.na = remove.na,
            verbose = verbose)
    }
    # check the variables correlations ####
    if (assess.variables) {
        se.obj <- assessVariablesCorrelation(
            se.obj = se.obj,
            bio.variables = bio.variables,
            uv.variables = uv.variables,
            cat.cor.coef = cat.cor.coef,
            cont.cor.coef = cont.cor.coef,
            assess.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose)
        bio.variables <- se.obj$bio.variables
        uv.variables <- se.obj$uv.variables
        se.obj <- se.obj$se.obj
    }
    # data transformation and normalization ####
    # data transformation and normalization ####
    printColoredMessage(
        message = '-- Data transformation and normalization:',
        color = 'magenta',
        verbose = verbose)
    ## apply log ####
    if (isTRUE(apply.log) & !is.null(pseudo.count)){
        printColoredMessage(
            message = paste0(
                'Applying log2 + ',
                pseudo.count,
                ' (pseudo.count) on the ',
                assay.name,
                ' data.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else if (isTRUE(apply.log) & is.null(pseudo.count)){
        printColoredMessage(
            message = paste0(
                'Applying log2 on the ',
                assay.name,
                ' data.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- log2(assay(x = se.obj, i = assay.name))
    } else if (isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0('The ', assay.name, ' data will be used without any log transformation.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- assay(x = se.obj, i = assay.name)
    }
    ## normalizations ####
    if(!is.null(normalization)){
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
    ## regress out unwanted variation ####
    if(!is.null(regress.out.uv.variables) & !is.null(normalization)){
        printColoredMessage(
            message = paste0(
                'The ',
                paste0(regress.out.uv.variables, collapse = ' & '),
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'Note, we do not recommend regressing out the ',
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
                'The',
                paste0(regress.out.uv.variables, collapse = ' & '),
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'Note: we do not recommend regressing out ',
                paste0(regress.out.uv.variables, collapse = ' & '),
                'if they are largely associated with the ',
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
    ## regress out biological variation ####
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
                'We do not recommend regressing out the ',
                paste0(regress.out.bio.variables, collapse = ' & '),
                'if they are largely associated with the ',
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
    # finding negative control genes ####
    ## step 1: highly affected by unwanted variation ####
    printColoredMessage(
        message = '-- Find genes that are highly affected by each sources of unwnated variation:',
        color = 'magenta',
        verbose = verbose)
    uv.var.class <- unlist(lapply(
        uv.variables,
        function(x) class(colData(se.obj)[[x]])))
    categorical.uv <- uv.variables[uv.var.class %in% c('factor', 'character')]
    continuous.uv <- uv.variables[uv.var.class %in% c('numeric', 'integer')]
    ### anova between genes and categorical sources of unwanted variation ####
    if(length(categorical.uv) > 0 ){
        if(!is.null(regress.out.bio.variables)){
            data.to.use <- expr.data.reg.bio
        } else data.to.use <- expr.data
        printColoredMessage(
            message = paste0(
                'Perform ANOVA between individual gene-level ',
                'expression and each categorical source of unwanted variation: ',
                paste0(categorical.uv, collapse = ' & '),
                '.'),
            color = 'blue',
            verbose = verbose
        )
        anova.genes.uv <- lapply(
            categorical.uv,
            function(x) {
                keep.samples <- findRepeatingPatterns(
                    vec = colData(se.obj)[[x]],
                    n.repeat = min.sample.for.aov)
                if(length(keep.samples) == 0){
                    stop(paste0(
                        'There are not enough samples to perfrom ANOVA between individual gene expression and the ',
                        x,
                        ' variable. Possible solutions is to lower min.sample.for.aov or remove',
                        x,
                        'from the uv.variables and re-run the function.'))
                } else if(length(keep.samples) == 1 ){
                    stop(paste0(
                        'There is only a single batch from ',
                        x,
                        ' that have enough samples ',
                        min.sample.for.aov,
                        '(min.sample.for.aov). Possible solutions is to lower min.sample.for.aov or remove'))
                } else if(length(keep.samples) != length(unique(colData(se.obj)[[x]])) ){
                    not.coverd <- unique(colData(se.obj)[[x]])[!unique(colData(se.obj)[[x]]) %in% keep.samples]
                    printColoredMessage(
                        message = paste0(
                            'Note, the ',
                            paste0(not.coverd, collapse = '&'),
                            ' batches do not have enough samples for the ANOVA analysis.'),
                        color = 'red',
                        verbose = verbose)
                }
                keep.samples <- colData(se.obj)[[x]] %in% keep.samples
                if(anova.method == 'aov'){
                    anova.gene.batch <- as.data.frame(row_oneway_equalvar(
                            x = data.to.use[ , keep.samples],
                            g = se.obj@colData[, x][keep.samples]))
                } else if(anova.method == 'welch.correction'){
                    anova.gene.batch <- as.data.frame(row_oneway_welch(
                            x = data.to.use[ , keep.samples],
                            g = se.obj@colData[, x][keep.samples]))
                }
                set.seed(2233)
                anova.gene.batch$ranked.genes <- rank(-anova.gene.batch[, 'statistic'], ties.method = 'random')
                anova.gene.batch
            })
        names(anova.genes.uv) <- categorical.uv
        rm(data.to.use)
    } else anova.genes.uv <- NULL
    ### correlation between genes and categorical sources of unwanted variation ####
    if(length(continuous.uv) > 0 ){
        if(!is.null(regress.out.bio.variables)){
            data.to.use <- expr.data.reg.bio
        } else data.to.use <- expr.data
        printColoredMessage(
            message = paste0(
                'Perform ',
                corr.method,
                ' correlation between individual gene-level ',
                'expression and each continuous source of unwanted variations: ',
                paste0(continuous.uv, collapse = '&'),
                '.'),
            color = 'blue',
            verbose = verbose
        )
        if(ncol(se.obj) <= min.sample.for.correlation){
            stop(paste0('There are not enough samples (min.sample.for.correlation:',
                       min.sample.for.correlation,
                       ') to perform correlation analysis.',
                       ' A possible soultion in to lower min.sample.for.correlation.'))
        }
        corr.genes.uv <- lapply(
            continuous.uv,
            function(x) {
                corr.genes.var <- as.data.frame(correls(
                    y = se.obj@colData[, x],
                    x = t(data.to.use),
                    type = corr.method,
                    a = a ,
                    rho = rho))
                corr.genes.var <- cbind(
                    round(x = corr.genes.var[, 1:4], digits = 3),
                    corr.genes.var[, 5, drop = FALSE])
                set.seed(2233)
                corr.genes.var$ranked.genes <- rank(-abs(corr.genes.var[, 'correlation']), ties.method = 'random')
                row.names(corr.genes.var) <- row.names(data.to.use)
                corr.genes.var
            })
        names(corr.genes.uv) <- continuous.uv
        rm(data.to.use)
    } else corr.genes.uv <- NULL
    ## step2: not highly affected by biology ####
    printColoredMessage(
        message = '-- Find genes that are not highly affected by sources of biological variation:',
        color = 'magenta',
        verbose = verbose)
    bio.var.class <- unlist(lapply(
        bio.variables,
        function(x) class(colData(se.obj)[[x]]) ))
    continuous.bio <- bio.variables[bio.var.class %in% c('numeric', 'integer')]
    categorical.bio <- bio.variables[bio.var.class %in% c('factor', 'character')]
    ### anova between genes and categorical sources of biological variation ####
    if(length(categorical.bio) > 0 ){
        if(!is.null(normalization) & is.null(regress.out.uv.variables)){
            data.to.use <- expr.data.nor
        } else if(!is.null(regress.out.uv.variables) & !is.null(normalization)){
            data.to.use <- expr.data.reg.uv
        } else if(!is.null(regress.out.uv.variables) & is.null(normalization)){
            data.to.use <- expr.data.reg.uv
        } else if (!is.null(regress.out.uv.variables) & !is.null(normalization)){
            data.to.use <- expr.data
        }
        printColoredMessage(
            message = paste0(
                'Perform ANOVA between individual gene-level ',
                'expression and each categorical source of biological variation: ',
                paste0(categorical.bio, collapse = ' & '),
                '.'),
            color = 'blue',
            verbose = verbose
        )
        anova.genes.bio <- lapply(
            categorical.bio,
            function(x) {
                keep.samples <- findRepeatingPatterns(
                    vec = colData(se.obj)[[x]],
                    n.repeat = min.sample.for.aov)
                if(length(keep.samples) == 0){
                    stop(paste0(
                        'There are not enough samples to perfrom ANOVA between individual genes expression and the ',
                        x,
                        ' variable. Possible solutions is to lower min.sample.for.aov or remove',
                        x,
                        'from the bio.variables and re-run the function.'))
                } else if(length(keep.samples) == 1 ){
                    stop(paste0(
                        'There is only a single batch from ',
                        x,
                        ' that have enough samples ',
                        min.sample.for.aov,
                        '(min.sample.for.aov). Possible solutions is to lower min.sample.for.aov or remove',
                        x,
                        'from the bio.variables and re-run the function'))
                } else if(length(keep.samples) != length(unique(colData(se.obj)[[x]])) ){
                    not.coverd <- unique(colData(se.obj)[[x]])[!unique(colData(se.obj)[[x]]) %in% keep.samples]
                    printColoredMessage(
                        message = paste0(
                            'Note, the',
                            paste0(not.coverd, collapse = '&'),
                            ' groups do not have enough samples for the ANOVA analysis.'),
                        color = 'red',
                        verbose = verbose)
                }
                keep.samples <- colData(se.obj)[[x]] %in% keep.samples
                if(anova.method == 'aov'){
                    anova.genes <- as.data.frame(row_oneway_equalvar(
                        x = data.to.use[ , keep.samples],
                        g = se.obj@colData[, x][keep.samples]))
                } else if(anova.method == 'welch.correction'){
                    anova.genes <- as.data.frame(row_oneway_equalvar(
                        x = data.to.use[ , keep.samples],
                        g = se.obj@colData[, x][keep.samples]))
                }
                set.seed(2233)
                anova.genes$ranked.genes <- rank(anova.genes[, 'statistic'], ties.method = 'random')
                anova.genes
            })
        names(anova.genes.bio) <- categorical.bio
        rm(data.to.use)
    } else anova.genes.bio <- NULL
    ### correlation between genes and continuous sources of biological variation ####
    if(length(continuous.bio) > 0 ){
        if(!is.null(normalization) & is.null(regress.out.uv.variables)){
            data.to.use <- expr.data.nor
        } else if(!is.null(regress.out.uv.variables) & !is.null(normalization)){
            data.to.use <- expr.data.reg.uv
        } else if(!is.null(regress.out.uv.variables) & is.null(normalization)){
            data.to.use <- expr.data.reg.uv
        } else if (!is.null(regress.out.uv.variables) & !is.null(normalization)){
            data.to.use <- expr.data
        }
        ### gene-batch anova
        printColoredMessage(
            message = paste0(
                'Perform ',
                corr.method,
                ' correlation between individual gene-level ',
                'expression and each continuous sources of biological variation: ',
                paste0(continuous.bio, collapse = '&'),
                '.'),
            color = 'blue',
            verbose = verbose)
        if(ncol(se.obj) <= min.sample.for.correlation){
            stop(paste0('There are not enough samples (min.sample.for.correlation:',
                        min.sample.for.correlation,
                        ') to perform correlation analysis.',
                        ' A possible soultion in to lower min.sample.for.correlation.'))
        }
        corr.genes.bio <- lapply(
            continuous.bio,
            function(x) {
                corr.genes.var <- as.data.frame(correls(
                    y = se.obj@colData[, x],
                    x = t(data.to.use),
                    type = corr.method,
                    a = a ,
                    rho = rho))
                corr.genes.var <- cbind(
                    round(x = corr.genes.var[, 1:4], digits = 3),
                    corr.genes.var[, 5, drop = FALSE])
                row.names(corr.genes.var) <- row.names(data.to.use)
                set.seed(2233)
                corr.genes.var$ranked.genes <- rank(abs(corr.genes.var[, 'correlation']), ties.method = 'random')
                corr.genes.var
            })
        names(corr.genes.bio) <- continuous.bio
    } else corr.genes.bio <- NULL

    # final selection ####
    printColoredMessage(
        message = '-- Selection of a set of genes as NCG:',
        color = 'magenta',
        verbose = verbose)
    ## prod, sum or average of ranks ####
    if(ncg.selection.method %in% c('Prod', 'Sum', 'Average')){
        all.tests <- c(
            'anova.genes.bio',
            'corr.genes.bio',
            'anova.genes.uv',
            'corr.genes.uv')
        ncg.selected <- lapply(
            all.tests,
            function(x){
                temp <- get(x)
                if(length(names(temp))!=0){
                    ranks.data <- lapply(
                        names(temp),
                        function(y){
                            temp[[y]]$ranked.genes
                        })
                    ranks.data <- do.call(cbind, ranks.data)
                    colnames(ranks.data) <- names(temp)
                    ranks.data
                }
            })
        ncg.selected <- do.call(cbind, ncg.selected)
        row.names(ncg.selected) <- row.names(se.obj)
        if(ncg.selection.method == 'Prod'){
            printColoredMessage(
                message = 'A set of NCG will be selected based on the product of ranks.',
                color = 'blue',
                verbose = verbose)
            stat.rank <- rowProds(ncg.selected)
            if(sum(is.infinite(stat.rank)) > 0){
                stop('The product of ranks results in infinity values.')
            }
        } else if(ncg.selection.method == 'Sum'){
            printColoredMessage(
                message = 'A set of NCG will be selected based on the sum of ranks.',
                color = 'blue',
                verbose = verbose)
            stat.rank <- rowSums(ncg.selected)
        } else if(ncg.selection.method == 'Average'){
            printColoredMessage(
                message = 'A set of NCG will be selected based on the average of ranks.',
                color = 'blue',
                verbose = verbose)
            stat.rank <- rowMeans(ncg.selected)
        }
        ncg.selected <- as.data.frame(ncg.selected)
        row.names(ncg.selected) <- row.names(se.obj)
        ncg.selected$stat.rank <- stat.rank
        ncg.selected$final.rank <- rank(ncg.selected$stat.rank)
        ncg.selected <- ncg.selected[order(ncg.selected$stat.rank, decreasing = FALSE) , ]
        ncg.selected <- row.names(ncg.selected[1:round(c(nb.ncg/100) * nrow(se.obj), digits = 0) , ])
        ncg.selected <- row.names(se.obj) %in% ncg.selected

        ## noneOverlap of ranks ####
    } else if (ncg.selection.method == 'noneOverlap'){
        printColoredMessage(
            message = 'A set of NCG is selected based on the noneOverlab approach.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'The non-overlap set of genes between top ',
                top.rank.bio.genes,
                '% of highly affected genes by the bioloigcal variation and top ',
                top.rank.uv.genes,
                '% of highly affected genes by the unwanted variation.'),
            color = 'blue',
            verbose = verbose)
        if(top.rank.bio.genes == 100){
            top.rank.bio.genes.nb <- nrow(se.obj)
        } else{
            top.rank.bio.genes.nb <- round(c(top.rank.bio.genes/100) * nrow(se.obj), digits = 0)
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
                            index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                            row.names(temp.data[[y]])[index] })))
                }
            })))
        top.rank.uv.genes <- round(c(top.rank.uv.genes/100) * nrow(se.obj), digits = 0)
        top.uv.genes <- unique(unlist(lapply(
            all.uv.tests,
            function(x){
                if(!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes
                            row.names(temp.data[[y]])[index] })))
                }
            })))
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        nb.ncg <- round(c(nb.ncg/100) * nrow(se.obj), digits = 0)
        if(length(ncg.selected) > nb.ncg){
            lo <- top.rank.uv.genes
            grid.nb <- round(c(grid.nb/100) * nrow(se.obj), digits = 0)
            pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
            while(length(ncg.selected) > nb.ncg | top.rank.uv.genes == 0){
                pro.bar$pause(0.1)$tick()$print()
                all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                top.uv.genes <- unique(unlist(lapply(
                    all.uv.tests,
                    function(x){
                        if(!is.null(x)){
                            temp.data <- get(x)
                            ranks.data <- unique(unlist(lapply(
                                names(temp.data),
                                function(y){
                                    index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes
                                    row.names(temp.data[[y]])[index] })))
                        }
                    })))
                top.rank.uv.genes <- top.rank.uv.genes - grid.nb
                if(top.rank.uv.genes > nrow(se.obj)){
                    message(' ')
                    printColoredMessage(
                        message = paste0(length(ncg.selected), ' genes are found based on the current parameters.'),
                        color = 'red', verbose = verbose)
                    stop('Any NCG cannot be found. Please lower either the values of top.rank.bio.genes or nb.ncg.')
                }
                ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                ncg.selected
            }
            ncg.selected <- row.names(se.obj) %in% ncg.selected
            top.rank.uv.genes <- round(top.rank.uv.genes/nrow(se.obj) * 100, digits = 0)
            if(top.rank.uv.genes >= 100){
                top.rank.uv.genes = 100
            }
            message(' ')
            printColoredMessage(
                message = paste0(
                    'The non-overlap set of genes between top ',
                    top.rank.bio.genes ,
                    '% of highly affected genes by the bioloigcal variation and top ',
                    top.rank.uv.genes,
                    '% of highly affected genes by the unwanted variation.'),
                color = 'blue',
                verbose = verbose)
        } else {
            message(' ')
            printColoredMessage(
                message = paste0(
                    length(ncg.selected),
                    ' genes are found based on the current parameters.'),
                color = 'blue',
                verbose = verbose)
            }
        ## AbsNoneOverlap of ranks ####
        } else if (ncg.selection.method == 'AbsNoneOverlap'){
        printColoredMessage(
            message = 'A set of NCG is selected based on the AbsNoneOverlab approach.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'The non-overlap set of genes between top ',
                top.rank.bio.genes,
                '% of highly affected genes by the bioloigcal variation and top ',
                top.rank.uv.genes,
                '% of highly affected genes by the unwanted variation.'),
            color = 'blue',
            verbose = verbose)
        if(top.rank.bio.genes == 100){
            top.rank.bio.genes.nb <- nrow(se.obj)
        } else{
            top.rank.bio.genes <- round(c(top.rank.bio.genes/100) * nrow(se.obj), digits = 0)
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
                            index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                            row.names(temp.data[[y]])[index] })))
                }
            })))
        top.rank.uv.genes <- round(c(top.rank.uv.genes/100) * nrow(se.obj), digits = 0)
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        top.uv.genes <- unique(unlist(lapply(
            all.uv.tests,
            function(x){
                if(!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes
                            row.names(temp.data[[y]])[index] })))
                }
            })))
        top.uv.genes <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        ncg.selected <- row.names(se.obj) %in% top.uv.genes
    }
    printColoredMessage(
        message = paste0('A set of ', sum(ncg.selected), ' genes are selected for NCG.'),
        color = 'blue',
        verbose = verbose)

    # assessment of the selected set of NCG ####
    ## pca on the NCG ####
    if(assess.ncg){
        printColoredMessage(
            message = '-- Assess the performance of the selected NCG set:',
            color = 'magenta',
            verbose = verbose)
        printColoredMessage(
            message = 'Perform PCA on only the selected genes as NCG.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'Explore the association of the first ',
                nb.pcs,
                '  with the ',
                paste0(c(bio.variables, uv.variables), collapse = ' & '),
                ' variables.'),
            color = 'blue',
            verbose = verbose)
        pca.data <- BiocSingular::runSVD(
            x = t(expr.data[ncg.selected, ]),
            k = nb.pcs,
            BSPARAM = bsparam(),
            center = TRUE,
            scale = FALSE)$u
        if(is.null(variables.to.assess.ncg)){
            variables.to.assess.ncg <- c(bio.variables, uv.variables)
        }
        all.corr <- lapply(
            variables.to.assess.ncg,
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
                        1 - prod(1 - cca$cor^2)})
                }
            })
        names(all.corr) <- variables.to.assess.ncg
        pcs <- Groups <- NULL
        pca.ncg <- as.data.frame(do.call(cbind, all.corr)) %>%
            dplyr::mutate(pcs = c(1:nb.pcs)) %>%
            tidyr::pivot_longer(
                -pcs,
                names_to = 'Groups',
                values_to = 'ls')
        pca.ncg <- ggplot(pca.ncg, aes(x = pcs, y = ls, group = Groups)) +
            geom_line(aes(color = Groups), size = 1) +
            geom_point(aes(color = Groups), size = 2) +
            xlab('PCs') +
            ylab (expression("Correlations")) +
            ggtitle('Assessment of the NCGs') +
            scale_x_continuous(breaks = (1:nb.pcs), labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', size = 1),
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
    printColoredMessage(
        message = '-- Save the NCGs:',
        color = 'magenta',
        verbose = verbose)
    out.put.name <- paste0(
        sum(ncg.selected),
        '|',
        paste0(bio.variables, collapse = '&'),
        '|',
        paste0(uv.variables, collapse = '&'),
        '|AnoCorrAs:',
        ncg.selection.method,
        '|',
        assay.name)
    if(save.se.obj == TRUE){
        printColoredMessage(
            message = '-- Save the selected set of NCG to the metadata of the SummarizedExperiment object.',
            color = 'blue',
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
            message = '------------The supervisedFindNcgAnoCorrAs function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj)
    } else{
        printColoredMessage(
            message = '-- The NCGs are outpputed as a logical vector.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The supervisedFindNcgAnoCorrAs function finished.',
            color = 'white',
            verbose = verbose)
        return(ncg.selected)
    }
}
