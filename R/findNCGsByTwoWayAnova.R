#' is used to find a set of negative control genes (NCG) using two-way ANOVA.

#' @author Ramyar Molania

#' @description
#' This function uses the two-way ANOVA approach to find a set of genes as negative control genes (NCG) for RUV-III-PRPS
#' normalization. Sources of known biological and unwanted variation should be specified. First, the function creates all
#' possible sample groups with respect to biological and unwanted variation separately. Then, the function uses the groups
#' as factors in two-way ANOVA to find genes that are highly affected by biological and unwanted variation separately.
#' Finally, the function selects genes that have possible high F-statistics for the unwanted variation and low F-statistics
#' for the biological variation. The function uses different approaches to perform the final selection.

#' @details
#' The function uses 5 ways to summarize two gene-level F-statistics obtained for the biological and unwanted variation
#' separately. The options are 'Prod', 'Sum', 'Average', 'AbsNoneOverlap' or 'noneOverlap'. The 'Prod', 'Sum', 'Average'
#' are based on the rank of the gene-level F-statistics. The 'prod' stands for the product of of the ranks, 'Sum' stands
#' for the sum of the ranks and the 'Average' stands for the average of the rank.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. Indicates a name of an assay in the SummarizedExperiment object. The selected assay should
#' be the one that will be used for RUV-III-PRPS normalization.
#' @param bio.variables Symbols. Indicates the columns names that contain known biological variables in the
#' SummarizedExperiment object. The biological variables cab be categorical and continuous. The continuous variable will
#' be grouped into 'nb.bio.clusters' based on a clustering method selected in 'bio.clustering.method'.
#' @param uv.variables Symbols. Indicates the columns names that contains UV variables in the SummarizedExperiment object.
#' The unwanted variables can be categorical and continuous. The continuous variables will
#' be grouped into 'nb.uv.clusters' based on a clustering method selected in 'uv.clustering.method'.
#' @param nb.ncg Numeric. Indicates how many genes should be selected as NCG. The value is the percentage of the total
#' genes in the SummarizedExperiment object. The default is 10 percent.
#' @param ncg.selection.method Symbol. Indicates how to select a set genes as NCG after performing two-way ANOVA. For individual
#' genes, the two-way ANOVA calculates F-statistics for biological and unwanted variation factors separately. An ideal NCG
#' set should have high F-statistics for the unwanted variation variables and low F-statistics for the biological variables.
#' The function ranks the F-statistics obtained for the biological variable and negative of the F-statistics obtained for
#' the unwanted variables. Then this functions offers 5 ways to summarize the ranks of the two F-statistics. Prod' is the
#' product of the ranks. 'Sum', is the sum of the ranks. 'Average' is the average of the ranks. 'AbsNoneOverlap' is the
#' none overlapped genes of the 'top.rank.uv.genes' and 'top.rank.bio.genes'. 'noneOverlap' is the none overlapped genes
#' of the 'top.rank.uv.genes' and at least 'top.rank.bio.genes'.
#' @param grid.nb Numeric. Indicates the percentage for grid search when the ncg.selection.method is noneOverlap'. In the
#' 'noneOverlap' approach, the grid search starts with the initial top.rank.uv.genes' value an add the grid.nb in each
#' loop to find the 'nb.ncg'.
#' @param top.rank.bio.genes Numeric. Indicates the percentage of top ranked genes that are highly affected by the biological
#' variation. This is required to be specified when the 'ncg.selection.method' is either 'noneOverlap' or 'AbsNoneOverlap'.
#' @param top.rank.uv.genes Numeric. Indicates the percentage of top ranked genes that are highly affected by the unwanted
#' variation variables. This is required to be specified when the 'ncg.selection.method' is either 'noneOverlap' or
#' 'AbsNoneOverlap'.
#' @param bio.clustering.method Symbols. Indicates which clustering methods should be used to group continuous sources
#' of biological variation.
#' @param nb.bio.clusters Numeric.Indicates the number of clusters for each continuous  sources of biological variation,
#' by default it is set to 2.
#' @param uv.clustering.method Symbols.Indicates which clustering method should be used to group continuous sources of UV.
#' @param nb.uv.clusters Numeric.Indicates the number of clusters for each continuous sources of UV, by default it is set to 2.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation. The
#' default is 1.
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected NCG or not. This analysis involves
#' principal component analysis on only the selected NCG and then explore the R^2 or vector correlation between the 'nb.pcs'
#' first principal components and with the biological and unwanted variables.
#' @param variables.to.assess.ncg Symbols. Indicates the column names of the SummarizedExperiment object that contain
#' variables whose association with the selected genes as NCG needs to be evaluated. The default is NULL. This means all
#' the variables in the 'bio.variables' and 'uv.variables' will be assessed.
#' @param nb.pcs Numeric. Indicates the number of the first principal components on the selected NCG to be used to assess
#' the performance of NCGs. The default is 5.
#' @param center Logical. Indicates whether to center the data before applying principal component analysis or not. The
#' default is TRUE.
#' @param scale Logical. Indicates whether to scale the data before applying Principal component analysis. The default is FALSE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object or not.
#' @param assess.variables Logical. Indicates whether to assess the correlation between biological and unwanted variation
#' variables separately.
#' @param remove.na Symbol. Indicates whether to remove NA or missing values from either the 'assays', the 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be excluded. If
#' 'sample.annotation' is selected, the samples that contains NA or missing values for any 'bio.variables' and
#' 'uv.variables' will be excluded. By default, it is set to both'.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment object or
#' to output the result. By default it is set to TRUE.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @return Either the SummarizedExperiment object containing the a set of negative control genes in the metadata  or a
#' logical vector of the selected negative control genes.

#' @importFrom SummarizedExperiment assay
#' @importFrom BiocSingular runSVD bsparam
#' @importFrom fastDummies dummy_cols
#' @importFrom dplyr mutate progress_estimated
#' @importFrom tidyr pivot_longer
#' @importFrom stats aov
#' @import ggplot2
#' @export

supervisedFindNcgTWAnova <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        nb.ncg = 10,
        ncg.selection.method = 'noneOverlap',
        grid.nb = 1,
        top.rank.bio.genes = 70,
        top.rank.uv.genes = 70,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 3,
        uv.clustering.method = 'kmeans',
        nb.uv.clusters = 3,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.ncg = TRUE,
        variables.to.assess.ncg = NULL,
        nb.pcs = 5,
        center = TRUE,
        scale = FALSE,
        assess.se.obj = TRUE,
        assess.variables = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(
        message = '------------The supervisedFindNcgTWAnova function starts:',
        color = 'white',
        verbose = verbose)

    # check inputs ####
    if(!is.vector(assay.name)){
        stop('The "assay.name" must be a single assay name.')
    } else if (length(assay.name) > 1){
        stop('The "assay.name" must be a single assay name.')
    } else if (is.null(bio.variables)){
        stop('The "bio.variables" cannot be empty.')
    } else if (is.null(uv.variables)){
        stop('The "uv.variables" cannot be empty.')
    } else if(!is.vector(bio.variables) | !is.vector(uv.variables) ){
        stop('The "uv.variables" and "bio.variables" must be a vector of variables names.')
    } else if (length(intersect(bio.variables, uv.variables)) > 0){
        stop('Any variables must be either in "bio.variables" or "uv.variables".')
    } else if (nb.ncg >= 100 | nb.ncg <= 0){
        stop('The "nb.ncg" must be a positve value 0 < nb.ncg =< 100.')
    } else if (!ncg.selection.method %in% c('Prod', 'Sum', 'Average', 'AbsNoneOverlap', 'noneOverlap')){
        stop('The ncg.selection.method should be one of "Prod", "Sum", "Average", "AbsNoneOverlap" or "noneOverlap".')
    } else if (top.rank.bio.genes > 100 | top.rank.bio.genes <= 0){
        stop('The "top.rank.bio.genes" should be a positve value  0 < top.rank.bio.genes < 100.')
    } else if (top.rank.uv.genes > 100 | top.rank.uv.genes <= 0){
        stop('The "top.rank.uv.genes" should be a positve value  0 < top.rank.uv.genes < 100.')
    }

    # check the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = unique(c(bio.variables, uv.variables, variables.to.assess.ncg)),
            remove.na = remove.na,
            verbose = verbose)
    }

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

    # create all possible homogeneous biological groups ####
    printColoredMessage(
        message = '-- Create all possible major homogeneous biological groups:',
        color = 'magenta',
        verbose = verbose)
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

    # create all possible homogeneous groups with respect to sources of unwanted variation ####
    printColoredMessage(
        message = '-- Create all possible major groups with respect to sources of unwanted variation:',
        color = 'magenta',
        verbose = verbose)
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

    # perform Two_way ANOVA between ####
    printColoredMessage(
        message = '-- Perform two_way ANOVA:',
        color = 'magenta',
        verbose = verbose)
    printColoredMessage(
        message = paste0('This is between all individual gene expression and considering both biological and UV variables created above as factors.'),
        color = 'blue',
        verbose = verbose)
    all.aov <- aov(t(expr.data) ~ all.bio.groups + all.uv.groups)
    all.aov <- summary(all.aov)
    all.aov <- as.data.frame(
        t(sapply(c(1:nrow(se.obj)),
                 function(x) all.aov[[x]]$`F value`[1:2]))
        )
    colnames(all.aov) <- c('Biology', 'UV')
    row.names(all.aov) <- row.names(se.obj)
    all.aov$bio.rank <- rank(all.aov$Biology)
    all.aov$uv.rank <- rank(-all.aov$UV)
    # selection of NCG ####
    printColoredMessage(
        message = '-- Selection of a set of genes as NCG:',
        color = 'magenta',
        verbose = verbose)
    if (ncg.selection.method %in% c('Prod', 'Average', 'Sum')){
        if(ncg.selection.method == 'Prod'){
            printColoredMessage(
                message = 'A set of NCG will be selected based on the product of ranks.',
                color = 'blue',
                verbose = verbose)
            all.aov$all.rank <- all.aov$bio.rank * all.aov$uv.rank
            if(sum(is.infinite(all.aov$all.rank)) > 0){
                stop('The product of ranks results in infinity values.')
            }
        } else if (ncg.selection.method == 'Average'){
            printColoredMessage(
                message = 'A set of NCG will be selected based on the average of ranks.',
                color = 'blue',
                verbose = verbose)
            all.aov$all.rank <- rowMeans(all.aov[ , c('bio.rank', 'uv.rank')])
        } else if(ncg.selection.method == 'Sum'){
            printColoredMessage(
                message = 'A set of NCG will be selected based on the sum of ranks.',
                color = 'blue',
                verbose = verbose)
            all.aov$all.rank <- rowSums(all.aov[ , c('bio.rank', 'uv.rank')])
        }
        all.aov <- all.aov[order(all.aov$all.rank, decreasing = FALSE), ]
        ncg.selected <- row.names(all.aov)[1:round(c(nb.ncg/100) * nrow(se.obj), digits = 0)]
        ncg.selected <- row.names(se.obj) %in% ncg.selected
    } else if (ncg.selection.method == 'AbsNoneOverlap'){
        printColoredMessage(
            message = 'A set of NCG is selected based on the AbsNoneOverlap approach.',
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
        top.rank.bio.genes <- round(c(top.rank.bio.genes/100) * nrow(se.obj), digits = 0)
        top.bio.genes <- all.aov$bio.rank > c(nrow(se.obj) - top.rank.bio.genes)
        top.bio.genes <- row.names(all.aov)[top.bio.genes]
        top.rank.uv.genes <- round(c(top.rank.uv.genes/100) * nrow(se.obj), digits = 0)
        top.uv.genes <- all.aov$uv.rank <  top.rank.uv.genes
        top.uv.genes <- row.names(all.aov)[top.uv.genes]
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        ncg.selected <- row.names(se.obj) %in% ncg.selected
        if(sum(ncg.selected) == 0){
            stop(paste0('The non-oevrlab genes between genes that are highly affected by batch and not biology is 0.',
                        'Please increase either "top.rank.bio.genes" or "top.rank.uv.genes" values'))
        }
    } else if (ncg.selection.method == 'noneOverlap'){
        printColoredMessage(
            message = 'A set of NCG is selected based on the noneOverlap approach.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0(
                'The non-overlap set of genes between top ',
                top.rank.bio.genes ,
                '% of highly affected genes by the bioloigcal variation and at least top ',
                top.rank.uv.genes,
                '% of highly affected genes by the unwanted variation.'),
            color = 'blue',
            verbose = verbose)
        nb.ncg <- round(c(nb.ncg/100) * nrow(se.obj), digits = 0)
        top.rank.bio.genes.no <- round(c(top.rank.bio.genes/100) * nrow(se.obj), digits = 0)
        top.bio.genes <- all.aov$bio.rank > c(nrow(se.obj) - top.rank.bio.genes.no)
        top.bio.genes <- row.names(all.aov)[top.bio.genes]

        top.rank.uv.genes <- round(c(top.rank.uv.genes/100) * nrow(se.obj), digits = 0)
        top.uv.genes <- all.aov$uv.rank <  top.rank.uv.genes
        top.uv.genes <- row.names(all.aov)[top.uv.genes]
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        if(length(ncg.selected) >= nb.ncg){
            grid.nb <- round(c(grid.nb/100) * nrow(se.obj), digits = 0)
            pro.bar <- progress_estimated(round(top.rank.uv.genes/grid.nb, digits = 0) + 1)
            while(length(ncg.selected) > nb.ncg | top.rank.uv.genes == 0){
                pro.bar$pause(0.1)$tick()$print()
                top.uv.genes <- all.aov$uv.rank <  top.rank.uv.genes
                top.uv.genes <- row.names(all.aov)[top.uv.genes]
                top.rank.uv.genes <- top.rank.uv.genes - grid.nb
                ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
            }
            if(top.rank.uv.genes == 0){
                top.rank.uv.genes = 0
            } else {
                top.rank.uv.genes <- round(top.rank.uv.genes/nrow(se.obj) * 100, digits = 0)
                if(top.rank.uv.genes >= 100){
                    top.rank.uv.genes = 100
                }
            }
            message(' ')
            printColoredMessage(
                message = paste0(
                    'The non-overlap set of genes between top ',
                    top.rank.bio.genes,
                    '% of highly affected genes by the bioloigcal variation and top ',
                    top.rank.uv.genes,
                    '% of highly affected genes by the unwanted variation.'),
                color = 'blue',
                verbose = verbose)
            ncg.selected <- row.names(se.obj) %in% ncg.selected
        } else {
            message(' ')
            printColoredMessage(
                    message = paste0(
                        length(ncg.selected),
                        ' genes are found based on the current parameters.'),
                    color = 'red',
                    verbose = verbose)
        }
    }
    printColoredMessage(
        message = paste0(sum(ncg.selected), ' genes are selected as negative control genes.'),
        color = 'blue',
        verbose = verbose)

    # assessment of selected set of NCG  ####
    if(is.null(variables.to.assess.ncg)){
        variables.to.assess.ncg <- c(bio.variables, uv.variables)
    }
    printColoredMessage(
        message = '-- Assess the performance of selected NCG set:',
        color = 'magenta',
        verbose = verbose)
    printColoredMessage(
        message = paste0('Perform PCA on only selected genes as NCG.'),
        color = 'blue',
        verbose = verbose)
    printColoredMessage(
        message = paste0(
            'Explore the association of the first ',
            nb.pcs,
            '  PCs with the ',
            paste0(variables.to.assess.ncg, collapse = ' & '),
            ' variables.'),
        color = 'blue',
        verbose = verbose)
    pca.data <- BiocSingular::runSVD(
        x = t(expr.data[ncg.selected, ]),
        k = nb.pcs,
        BSPARAM = bsparam(),
        center = TRUE,
        scale = FALSE)$u
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
                    1 - prod(1 - cca$cor^2)
                    })
            }
        })
    names(all.corr) <- variables.to.assess.ncg
    pca.ncg <- as.data.frame(do.call(cbind, all.corr))
    pcs <- Groups <- NULL
    pca.ncg['pcs'] <- c(1:nb.pcs)
    pca.ncg <- tidyr::pivot_longer(
        data = pca.ncg,
        -pcs,
        names_to = 'Groups',
        values_to = 'ls')
    pca.ncg <- ggplot(pca.ncg, aes(x = pcs, y = ls, group = Groups)) +
        geom_line(aes(color = Groups), linewidth = 1) +
        geom_point(aes(color = Groups), size = 2) +
        xlab('PCs') +
        ylab (expression("Correlations")) +
        ggtitle('Assessment of the NCGs') +
        scale_x_continuous(breaks = (1:nb.pcs),labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = nb.pcs), limits = c(0,1)) +
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', linewidth = 1),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 10, angle = 25, hjust = 1),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 14),
            plot.title = element_text(size = 16)
        )
    if(verbose)  print(pca.ncg)
    # add results to the SummarizedExperiment object ####
    out.put.name <- paste0(
        sum(ncg.selected),
        '|',
        paste0(bio.variables, collapse = '&'),
        '|',
        paste0(uv.variables, collapse = '&'),
        '|TWAnova:',
        ncg.selection.method,
        '|',
        assay.name)
    if(save.se.obj == TRUE){
        printColoredMessage(
            message = '-- Saving the selected NCG to the metadata of the SummarizedExperiment object.',
            color = 'magenta',
            verbose = verbose)
        ## Check if metadata NCG already exists
        if(length(se.obj@metadata$NCG) == 0 ) {
            se.obj@metadata[['NCG']] <- list()
        }
        se.obj@metadata[['NCG']][['Supervised']][[out.put.name]] <- ncg.selected
        printColoredMessage(
            message = 'The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The supervisedFindNcgTWAnova function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj)
    } else{
        printColoredMessage(
            message = '-- The NCGs are outpputed a logical vectors.',
            color = 'magenta',
            verbose = verbose)
        printColoredMessage(
            message = '------------The supervisedFindNcgTWAnova function finished.',
            color = 'white',
            verbose = verbose)
        return(ncg.selected)
    }
}




