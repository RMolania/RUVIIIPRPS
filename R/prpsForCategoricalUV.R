#' is used to create pseudo-replicates of pseudo-samples (PRPS) for a categorical variable as source of unwanted variation.

#' @description
#' Distinct group of PRPS are created for each source of unwanted variation defined in the 'main.uv.variable' argument within
#' an homogeneous group of samples. The grouping of samples are created based on each biological subtype defined using the 'bio.variables'
#' using the 'bio.clustering.method' selected and it might be also be combined with 'other.uv.variables' if requested.
#' By default 'other.uv.variables' is set up to NULL, meaning that the homogeneous grouping of samples created is based solely on the biological subtype(s)
#' defined in 'bio.variables'. If 'other.uv.variables' is not set to NULL, the creation of homogeneous groups of samples is based on the biological
#' subtype defined in 'bio.variables' combined with 'other.uv.variables' using the 'other.uv.clustering.method' selected.
#' A pseudo-sample will be created for each homogeneous group of samples by averaging all the samples within that group of samples.
#' For example to correct for batch effect arising related to a batch variable defined in the 'main.uv.variable' argument, several group
#' of pseudo-samples will be created, one for each homogeneous group of samples and for each batch.
#' All those pseudo-samples belonging to the same group across batches will be defined as pseudo-replicates which constitutes a PRPS set.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. Indicates a name of an assay in the SummarizedExperiment object. The selected assay should
#' be the one that will be used for RUV-III-PRPS normalization.
#' @param bio.variables Symbol. A symbol or a list of symbols indicating the column names that contain biological variables.
#' All categorical and continuous variables can be provided. If more than one symbols are specified, all possible
#' homogeneous biological groups will be created by the createHomogeneousBioGroups function.
#' @param main.uv.variable Symbol. Indicates a name of the column that contains a source of categorical variable of
#' unwanted variation that PRPS should be created for.
#' @param other.uv.variables String of the label of (a) categorical or continuous variable(s) used to define homogeneous
#' groups of samples such as library size, plates from colData(se). By default it is set to 'NULL' meaning the samples
#' will be assigned to homogeneous biological groups of samples using the 'bio.variables'.
#' @param min.sample.for.prps Numeric. Indicates the minimum number of samples to be averaged to create one pseudo-sample,
#'  by default is 3.
#' to create a PRPS set for each continuous variable. The minimum should be '2*min.sample.for.prps'. By default it is set to 6.
#' @param bio.clustering.method String of the clustering method to assign each sample to an homogeneous biological clusters/group using the function 'kmeans',
#' 'cut' from base or the function 'quantile'. By default it is to 'kmeans'.
#' @param other.uv.clustering.method String of the clustering method to assign each sample to an homogeneous clusters/group of samples based on the 'other.uv.variables'
#' using the function 'kmeans','cut' from base or the function 'quantile'.By default it is to 'kmeans'.
#' @param nb.bio.clusters Numeric. A value to specify the number of homogeneous biological clusters/groups of samples. By default it is set to 3.
#' @param nb.other.uv.clusters Numeric. A value to specify the number of groups for continuous sources of biological variation. The default is 2.
#' This means each continuous sources will be divided into 2 groups using the 'clustering.method'.
#' @param cat.cor.coef Vector of two numerical values. Indicates the cut-off of the correlation coefficient between each pair of
#' categorical variables. The first one is between each pair of 'uv.variables' and the second one is between each pair of 'bio.variables'.
#' The correlation is computed by the function ContCoef from the DescTools package. If the correlation of a pair of variable is higher than
#' the cut-off, then only the variable that has the highest number of factor will be kept and the other one will be excluded from the
#' remaining analysis. By default they are both set to 0.7.
#' @param cont.cor.coef Vector of two numerical values. Indicates the cut-off of the Spearman correlation coefficient between each pair of
#' continuous variables. The first one is between each pair of 'uv.variables' and the second one is between each pair of 'bio.variables'.
#' If the correlation of a pair of variable is higher than the cut-off, then only the variable that has the highest variance will
#' be kept and the other one will be excluded from the remaining analysis. By default they are both set to 0.7.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' By default it is set to TRUE.
#' @param check.prps.connectedness TTTTT
#' @param plot.output TTTTTTT
#' @param remove.na String. Indicates whether to remove NA or missing values from either the 'assays', the 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be excluded. If 'sample.annotation' is selected, the
#' samples that contains NA or missing values for any 'bio.variables' and 'uv.variables' will be excluded. By default, it is set to
#' 'both'.
#' @param assess.variables TTTTT
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result, by default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or
#' messages displayed during the execution of the functions, by default it is set to TRUE.

#' @return A SummarizedExperiment object or a list that containing the PRPS data.

#' @author Ramyar Molania

#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr count
#' @importFrom tidyr %>%
#' @import ggplot2
#' @export


prpsForCategoricalUV <- function(
        se.obj,
        assay.name,
        bio.variables,
        main.uv.variable,
        other.uv.variables = NULL,
        min.sample.for.prps = 3,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 2,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 2,
        check.prps.connectedness = TRUE,
        apply.log = TRUE,
        pseudo.count = 1,
        remove.na = 'both',
        assess.se.obj = TRUE,
        assess.variables = FALSE,
        cat.cor.coef = c(0.95, 0.95),
        cont.cor.coef = c(0.95, 0.95),
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE) {
    printColoredMessage(message = '------------The prpsForCategoricalUV function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs of the function ####
    if (length(assay.name) > 1) {
        stop('A single assay name should be provided.')
    } else if (is.null(assay.name)) {
        stop('The assay.name cannot be empty.')
    } else if(is.null(bio.variables)){
        stop('The "bio.variables" cannot be empty')
    } else if (length(main.uv.variable) > 1) {
        stop('The function can only take a single categorical uv variable for the "main.uv.variable" argument.')
    } else if (is.null(main.uv.variable)) {
        stop('The main.uv.variable cannot be empty.')
    } else if (!class(se.obj[[main.uv.variable]]) %in% c('character', 'factor')) {
        stop('The main.uv.variable should be a categorical source of unwanted variation, e.g., platform effects')
    } else if (length(unique(se.obj[[main.uv.variable]])) == 1) {
        stop('The main.uv.variable should have at least two levels.')
    } else if (min.sample.for.prps <= 1) {
        stop('The minimum value for the min.sample.for.prps is 2.')
    } else if(main.uv.variable %in% other.uv.variables){
        stop('The main.uv.variable should not be in the other.uv.variables.')
    } else if (pseudo.count < 0){
        stop('The value for pseudo.count can not be negative.')
    } else if (max(cat.cor.coef) > 1 | max(cont.cor.coef) > 1){
        stop('The maximum value for cat.cor.coef or cont.cor.coef cannot be more than 1.')
    }
    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(main.uv.variable, bio.variables, other.uv.variables),
            remove.na = remove.na)
    }

    # data transformation ####
    printColoredMessage(
        message = '-- Data transformation:',
        color = 'magenta',
        verbose = verbose)
    if (apply.log) {
        printColoredMessage(
            message = paste0(
                'Applying log2 + ',
                pseudo.count,
                ' (pseudo.count',
                ') on the ',
                assay.name,
                ' data before creating PRPS.'),
            color = 'blue',
            verbose = verbose)
        expre.data <- log2(assay(se.obj, assay.name) + pseudo.count)
    } else{
        printColoredMessage(
            message = paste0(
                'The ',
                assay.name,
                ' data will be used without log transformation for creating PRPS.'),
            color = 'blue',
            verbose = verbose)
        expre.data <- assay(se.obj, assay.name)
    }
    # assign homogeneous biological groups of samples to each sample ####
    if(length(bio.variables) > 1){
        printColoredMessage(
            message = '-- Create homogeneous biological group of samples:',
            color = 'magenta',
            verbose = verbose )
        homo.bio.groups <- createHomogeneousBioGroups(
            se.obj = se.obj,
            bio.variables = bio.variables,
            nb.clusters = nb.bio.clusters,
            clustering.method = bio.clustering.method,
            assess.se.obj = assess.se.obj,
            assess.variables = assess.variables,
            cont.cor.coef = cont.cor.coef,
            cat.cor.coef = cat.cor.coef,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose)
        if(sum(table(homo.bio.groups) == 1) == length(unique(homo.bio.groups))){
            stop('All the homogeneous biological group of samples have only 1 sample. PRPS cannot be created.')
        }
    } else{
        homo.bio.groups <- colData(se.obj)[[bio.variables]]
        if(length(unique(homo.bio.groups)) == 1 )
            printColoredMessage(
                message = 'The level of the "bio.variables" is only 1.',
                color = 'red',
                verbose = verbose)
    }

    # PRPS within homogeneous biological * uv populations ####
    if(!is.null(other.uv.variables)){
        printColoredMessage(
            message = '-- Create homogeneous sample groups with respect to the other.uv.variables:',
            color = 'magenta',
            verbose = verbose )
        ## create all possible sample groups with respect to other.uv.variables ####
        homo.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = other.uv.variables,
            nb.clusters = nb.other.uv.clusters,
            clustering.method = other.uv.clustering.method,
            assess.se.obj = FALSE,
            assess.variables = FALSE,
            cont.cor.coef = cont.cor.coef,
            cat.cor.coef = cat.cor.coef,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose)
        all.groups <- data.frame(
            homo.uv.groups = homo.uv.groups,
            homo.bio.groups = homo.bio.groups)
        printColoredMessage(
            message = '-- Combine all homogeneous samples groups with respect to both biological and unwanted variables:',
            color = 'magenta',
            verbose = verbose)
        bio.batch <- uv.group <- NULL
        all.groups$bio.batch <- paste(
            all.groups$homo.uv.groups,
            all.groups$homo.bio.groups,
            sep = '||')
        all.groups$uv.group <- se.obj[[main.uv.variable]]
        printColoredMessage(
            message = paste0(
                length(unique(all.groups$bio.batch)),
                ' of ',
                c(length(unique(all.groups$homo.bio.groups)) *length(unique(all.groups$homo.uv.groups)) ),
                ' maximum possible groups (',
                length(unique(all.groups$homo.bio.groups)),
                ' biological groups * ',
                length(unique(all.groups$homo.uv.groups)),
                ' unwanted variation groups) ',
                'groups are created.'),
            color = 'blue',
            verbose = verbose)
        samples.dis <- table(all.groups$bio.batch, all.groups$uv.group)
        selected.groups <- sum(rowSums(samples.dis >= min.sample.for.prps) > 1)
        printColoredMessage(
            message = paste0(
                selected.groups,
                ' sample groups of ',
                length(unique(all.groups$bio.batch)),
                ' homogeneous samples have at least ',
                min.sample.for.prps,
                ' samples in at least two batches of ',
                main.uv.variable, '.'),
            color = 'blue',
            verbose = verbose )
        # check connection between PRPS sets ####
        if(check.prps.connectedness){
            printColoredMessage(
                message = paste0(
                    '-- Check the connection between possible PRPS sets across batches of ',
                    main.uv.variable,
                    '.'),
                color = 'magenta',
                verbose = verbose )
            samples.dis <- checkPRPSconnectedness(
                data.input = samples.dis,
                min.samples = min.sample.for.prps,
                batch.name = main.uv.variable)
        }
        # check samples abundance ####
        printColoredMessage(
            message = '-- Check samples abundance of groups before create PRPS:',
            color = 'magenta',
            verbose = verbose)
        samples.dis <- samples.dis[rowSums(samples.dis >= min.sample.for.prps) > 1 , ]
        mean.samples <- round(mean(samples.dis[samples.dis >= min.sample.for.prps]), digits = 0)
        printColoredMessage(
            message = paste0(
                'The average sample size of groups to select ',
                min.sample.for.prps,
                ' samples for PRPS is ',
                mean.samples,
                '.'),
            color = 'red',
            verbose = verbose)
        printColoredMessage(
            message = 'Please note, the high average sample size, may result in a contamination of PRPS with some other variation.',
            color = 'red',
            verbose = verbose)
        printColoredMessage(
            message = 'To avoid this, you may consider increase number of clustering of continuous biological and unwanted variables.',
            color = 'red',
            verbose = verbose)
        # create PRPS sets ####
        prps.sets <- lapply(
            1:nrow(samples.dis),
            function(x) {
                selected.batches <- colnames(samples.dis)[samples.dis[x , ] >= min.sample.for.prps]
                ps.matrix <- sapply(
                    selected.batches,
                    function(y) {
                        index.samples <- all.groups$bio.batch == row.names(samples.dis)[x] & all.groups$uv.group == y
                        rowMeans(expre.data[, index.samples])
                    })
                colnames(ps.matrix) <- rep(
                    paste(row.names(samples.dis)[x], main.uv.variable, sep = '||'),
                    ncol(ps.matrix))
                ps.matrix
            })
        prps.sets <- do.call(cbind, prps.sets)
        printColoredMessage(
            message = paste0(
                length(unique(colnames(prps.sets))),
                ' PRPS sets with the total number of ',
                ncol(prps.sets),
                ' pseudo-samples are created.'),
            color = 'blue',
            verbose = verbose)
    }
    # PRPS within just homogeneous biological populations ####
    if(is.null(other.uv.variables)) {
        all.groups <- data.frame(
            uv.group = se.obj[[main.uv.variable]],
            bio.groups = homo.bio.groups
            )
        samples.dis <- table(all.groups$bio.groups, all.groups$uv.group)

        ## check connection between PRPS sets ####
        if(check.prps.connectedness){
            printColoredMessage(
                message = paste0(
                    '-- Check the connection between possible PRPS sets across batches of ',
                    main.uv.variable,
                    '.'),
                color = 'magenta',
                verbose = verbose )
            samples.dis <- checkPRPSconnectedness(
                data.input = samples.dis,
                min.samples = min.sample.for.prps,
                batch.name = main.uv.variable)
        }
        # check samples abundance ####
        printColoredMessage(
            message = '-- Check samples abundance of groups before create PRPS:',
            color = 'magenta',
            verbose = verbose)
        samples.dis <- samples.dis[rowSums(samples.dis >= min.sample.for.prps) > 1 , ]
        samples.dis <- samples.dis[rowSums(samples.dis >= min.sample.for.prps) > 1 , ]
        mean.samples <- round(mean(samples.dis[samples.dis >= min.sample.for.prps]), digits = 0)
        printColoredMessage(
            message = paste0('The average sample size of groups to select ', min.sample.for.prps, ' samples for PRPS is ', mean.samples, '.'),
            color = 'red',
            verbose = verbose)
        printColoredMessage(
            message = 'Please note, the high average sample size, may result in a contamination of PRPS with some other variation.',
            color = 'red',
            verbose = verbose)
        printColoredMessage(
            message = 'To avoid this, you may consider increase the number of clusters of continuous biological and unwanted variables.',
            color = 'red',
            verbose = verbose)

        # create PRPS sets ####
        prps.sets <- lapply(
            1:nrow(samples.dis),
            function(x) {
                selected.batches <- colnames(samples.dis)[samples.dis[x , ] >= min.sample.for.prps]
                ps.matrix <- sapply(
                    selected.batches,
                    function(y) {
                        index.samples <- all.groups$bio.groups == row.names(samples.dis)[x] &
                            all.groups$uv.group == y
                        rowMeans(expre.data[, index.samples])
                    })
                colnames(ps.matrix) <- rep(
                    paste(row.names(samples.dis)[x], main.uv.variable, sep = '||'),
                    ncol(ps.matrix))
                ps.matrix
            })
        prps.sets <- do.call(cbind, prps.sets)
        printColoredMessage(
            message = paste0(
            length(unique(colnames(prps.sets))),
            ' PRPS sets with the total number of ',
            ncol(prps.sets),
            ' pseudo-samples are created.'),
            color = 'blue',
            verbose = verbose)
    }
    # plot output #####
    if (plot.output) {
        ## PRPS map plot
        printColoredMessage(message = '-- Plot PRPS map:',
                            color = 'magenta',
                            verbose = verbose)
        if(!is.null(other.uv.variables)){
            groups <- all.groups$bio.batch
        } else groups <- homo.bio.groups
        catvar <- use <- n <- NULL
        info <- as.data.frame(colData(se.obj))
        info$catvar <- as.factor(paste0(info[, main.uv.variable]))
        info$groups <- as.factor(groups)
        df.count <- info %>% dplyr::count(catvar, groups)
        df.count$use <- 'unselected'
        df.count$use[df.count$n >= min.sample.for.prps] <- 'Selected'
        p <- ggplot(df.count, aes(x = catvar, y = groups)) +
            geom_count(aes(color = use)) +
            geom_text(aes(
                label = n,
                hjust = 0.5,
                vjust = 0.5
            )) +
            xlab(main.uv.variable) +
            ylab('Homogeneous groups') +
            theme_bw() +
            theme(
                axis.line = element_line(colour = 'black', size = .85),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                axis.text.x = element_text(
                    size = 10,
                    angle = 45,
                    hjust = 1
                ),
                axis.text.y = element_text(
                    size = 12,
                    angle = 45,
                    hjust = 1
                ),
                legend.position = 'none'
            )
        if(verbose) print(p)
    }
    # saving the output ####
    if(!is.null(other.uv.variables)) {
        out.put.name <- paste0(
            main.uv.variable,
            '|:',
            paste0(bio.variables, collapse = '&'),
            '|:',
            paste0(other.uv.variables, collapse = '&'),
            '|',
            assay.name)
    } else{
        out.put.name <- paste0(
            main.uv.variable,
            '|:',
            paste0(bio.variables, collapse = '&'),
            '|',
            assay.name)
    }
    if (save.se.obj) {
        ## check if metadata PRPS already exists
        if (length(se.obj@metadata) == 0) {
            se.obj@metadata[['PRPS']] <- list()
        }
        ## check if metadata PRPS already exists
        if (!'PRPS' %in% names(se.obj@metadata)) {
            se.obj@metadata[['PRPS']] <- list()
        }
        ## check if metadata PRPS already exist for supervised
        if (!'supervised' %in% names(se.obj@metadata[['PRPS']])) {
            se.obj@metadata[['PRPS']][['supervised']] <- list()
        }
        ## Check if metadata PRPS already exist for supervised
        se.obj@metadata[['PRPS']][['supervised']][[out.put.name]] <- prps.sets
        printColoredMessage(
            message = paste0(
                'The PRPS are saved to metadata@PRPS$supervised: ',
                out.put.name,
                '.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The prpsForCategoricalUV function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else {
        printColoredMessage(message = '------------The prpsForCategoricalUV function finished.',
                            color = 'white',
                            verbose = verbose)
        return(prps.sets)
    }
}


#
#     ### Add plots to SummarizedExperiment object
#     printColoredMessage(message = '### Saving the PRPS map plot to the metadata of the SummarizedExperiment object.',
#                         color = 'magenta',
#                         verbose = verbose)
#     ## Check if metadata plot already exist
#     if (length(se.obj@metadata) == 0) {
#         se.obj@metadata[['plot']] <- list()
#     }
#     ## Check if metadata plot already exist
#     if (!'plot' %in% names(se.obj@metadata)) {
#         se.obj@metadata[['plot']] <- list()
#     }
#     ## Check if metadata plot already exist for this metric
#     if (!'PRPS' %in% names(se.obj@metadata[['plot']])) {
#         se.obj@metadata[['plot']][['PRPS']] <- list()
#     }
#
#     for (v in categorical.uv) {
#         ## Save the new plot
#         se.obj@metadata[['plot']][['PRPS']][[v]] <-
#             plot.all[[v]]
#     }
#
# }
