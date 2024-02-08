#' create PRPS sets across batches using integration anchors.

#' @author Ramyar Molania

#' @description
#' This functions employs the FindIntegrationAnchors function from the Seurat R package to create PRPS sets for the
#' RUV-III normalization of RNA-seq data. This function can be used in situations in which biological variation are unknown.


#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. A symbol indicating the assay name within the SummarizedExperiment object for the creation
#' of PRPS data. The assay must be the one that will be used as data input for the RUV-III-PRPS normalization.
#' @param uv.variable Symbol. A symbol specifying the name of columns in the sample annotation of the SummarizedExperiment
#' object. This variable can to be categorical or continuous variable. If a continuous variable is provide, this will be
#' divided into groups using the clusteintf method.
#' @param clustering.method Symbol.A symbol indicating the choice of clustering method for grouping the 'uv.variable'
#' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
#' continuous variable. The default is 3.
#' @param hvg Vector. A vector containing the names of highly variable genes. These genes will be utilized to identify
#' anchor samples across different batches. The default value is set to 'NULL'.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is 'TRUE'.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay before applying
#' log transformation to avoid -Inf for raw counts that are equal to 0. The default is 1.
#' @param anchor.features Numeric. A numeric value indicating the provided number of features to be used in anchor finding.
#' The default is 2000. We refer to the FindIntegrationAnchors R function for more details.
#' @param scale Logical. Whether or not to scale the features provided. Only set to FALSE if you have previously scaled
#' the features you want to use for each object in the object.list.
#' @param min.prps.samples Numeric. The minimum number of samples to be averaged to create a pseudo-sample. The default
#' is 3. The minimum value is 2.
#' @param max.prps.samples Numeric value indicating the maximum number of samples to be averaged for creating a pseudo-sample.
#' The default is 'inf'. Please note that averaging a high number of samples may lead to inaccurate PRPS estimates.
#' @param max.prps.sets Numeric. The maximum number of PRPS sets across batches. The default in 10.
#' @param min.score Numeric. A cut off to filter anchors based the scores that the FindIntegrationAnchors calculates.
#' The default is NULL, indicating no filtration will be applied based on the scores. The value is between 0 or 1.
#' @param normalization.method Symbol. Indicate which normalization methods should be used before finding the anchors.
#' The options are "LogNormalize" or "SCT". The default is "LogNormalize".
#' @param sct.clip.range Numeric. Numeric of length two specifying the min and max values the Pearson residual will be
#' clipped to. The default is 'NULL'. We refer to the FindIntegrationAnchors R function for more details.
#' @param reduction Symbol. Indicates which dimensional reduction to perform when finding anchors. The options are "cca":
#' canonical correlation analysis, "rpca": reciprocal PCA and "rlsi": Reciprocal LSI. The default is "cca".
#' @param l2.norm Logical. Indicates whether to perform L2 normalization on the CCA samples embeddings after dimensional
#' reduction or not. The default is 'TRUE'.
#' @param dims Numeric. Indicates which dimensions to use from the CCA to specify the neighbor search space. Th default
#' is 10.
#' @param k.anchor Numeric. How many neighbors (k) to use when picking anchors. Th default is 3.
#' @param k.filter Numeric. How many neighbors (k) to use when filtering anchors. Th default is 200.
#' @param k.score Numeric. How many neighbors (k) to use when scoring anchors. Th default is 30.
#' @param max.features Numeric. The maximum number of features to use when specifying the neighborhood search space in
#' the anchor filtering.Th default is 30.
#' @param nn.method Symbol.Method for nearest neighbor finding. Options include: "rann", "annoy". The defauly is "annoy".
#' @param n.trees Numeric. More trees gives higher precision when using annoy approximate nearest neighbor search
#' @param eps Numeric. Error bound on the neighbor finding algorithm (from RANN/Annoy)
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment
#' object or to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom Seurat VariableFeatures FindIntegrationAnchors
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom stats setNames
#' @importFrom purrr map_df
#' @importFrom Matrix rowMeans
#' @export

createUnSupervisedPRPSbyAnchors <- function(
        se.obj,
        assay.name,
        uv.variable,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        hvg = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        anchor.features = 2000,
        scale = TRUE,
        min.prps.samples = 3,
        max.prps.samples = 'inf',
        max.prps.sets = 10,
        min.score = NULL,
        normalization.method = "LogNormalize",
        sct.clip.range = NULL,
        reduction = "cca",
        l2.norm = TRUE,
        dims = 1:10,
        k.anchor = 3,
        k.filter = 200,
        k.score = 30,
        max.features = 200,
        nn.method = "annoy",
        n.trees = 50,
        eps = 0,
        assess.se.obj = TRUE,
        remove.na = 'assays',
        save.se.obj = TRUE,
        verbose = TRUE) {
    printColoredMessage(message = '------------The createUnSupervisedPRPSbyAnchors function starts:',
                        color = 'white',
                        verbose = verbose)
    # check input #####
    if (is.list(assay.name)) {
        stop('The "assay.name" must be the name of an assay in the SummarizedExperiment object.')
    } else if (length(assay.name) > 1) {
        stop('The "assay.name" must be a single name of an assay in the SummarizedExperiment object..')
    } else if (is.null(uv.variable)) {
        stop('The "uv.variable" cannot be empty. Note, unknown sources of UV can be found by the "identifyUnknownUV" function.')
    } else if (length(uv.variable) > 1) {
        stop('The "uv.variable" must be a single variable name in the SummarizedExperiment object.')
    } else if (!uv.variable %in% colnames(colData(se.obj))) {
        stop('The "uv.variable" cannot be found in the SummarizedExperiment object')
    } else if (k.anchor == 0) {
        stop('The k.anchor cannot be 0.')
    } else if (min.score > 1 | min.score < 0){
        stop('The value of "min.score" must be between 0 and 1.')
    } else if (!normalization.method %in% c("LogNormalize", "SCT")) {
        stop('The "normalization.method" must be one of the "LogNormalize" or "SCT".')
    } else if (!reduction %in% c("cca", "rpca", 'rlsi')) {
        stop('The "reduction" should be one of the "cca", "rpca" or "rlsi"')
    } else if (sum(hvg %in% row.names(se.obj)) != length(hvg)) {
        stop('All the "hvg" genes are not found in the SummarizedExperiment object.')
    }

    # check the SummarizedExperiment ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(uv.variable),
            remove.na = 'assays',
            verbose = verbose)
    }
    ini.variable <- se.obj[[uv.variable]]
    # check the main uv variable ####
    if (class(se.obj[[uv.variable]]) %in% c('integer', 'numeric')) {
        printColoredMessage(
            message = paste0('Then, each source will be divided into ', nb.clusters, ' groups using ', clustering.method, '.'),
            color = 'blue',
            verbose = verbose
        )
        if (clustering.method == 'kmeans') {
            set.seed(3456)
            uv.cont.clusters <- kmeans(
                x = colData(se.obj)[[uv.variable]],
                centers = nb.clusters,
                iter.max = 1000)
            se.obj[[uv.variable]] <- factor(x = paste0(uv.variable, '_group', uv.cont.clusters$cluster))
        } else if (clustering.method == 'cut') {
            uv.cont.clusters <- as.numeric(cut(
                x = colData(se.obj)[[uv.variable]],
                breaks = nb.clusters,
                include.lowest = TRUE))
            se.obj[[uv.variable]] <- factor(x = paste0(uv.variable, '_group', uv.cont.clusters))
        } else if (clustering.method == 'quantile') {
            quantiles <- quantile(
                x = colData(se.obj)[[uv.variable]],
                probs = seq(0, 1, 1 / nb.clusters))
            uv.cont.clusters <- as.numeric(cut(
                x = colData(se.obj)[[uv.variable]],
                breaks = quantiles,
                include.lowest = TRUE))
            se.obj[[uv.variable]] <- factor(x = paste0(uv.variable, '_group', uv.cont.clusters))
        }
    } else if (is.factor(se.obj[[uv.variable]])) {
            se.obj[[uv.variable]] <- factor(x = se.obj[[uv.variable]])
    }

    # find anchors ####
    ## split the data into groups ####
    printColoredMessage(
        message = paste0('- Split the SummarizedExperiment object into ', length(levels(se.obj[[uv.variable]])), ' groups.'),
        color = 'blue',
        verbose = verbose
    )
    groups <- levels(se.obj[[uv.variable]])
    all.seurat.objects <- lapply(
        groups,
        function(x) {
            samples.index <- se.obj[[uv.variable]] == x
            seu.obj <- SeuratObject::CreateSeuratObject(
                counts = assay(x = se.obj[, samples.index], i = assay.name),
                project = x)
            if (!is.null(hvg))
                Seurat::VariableFeatures(seu.obj) <-hvg
            return(seu.obj)
        })
    names(all.seurat.objects) <- groups
    all.samples.index <- c(1:ncol(se.obj))

    ## find anchors  ####
    printColoredMessage(message = '- Apply the "FindIntegrationAnchors" function.',
                        color = 'blue',
                        verbose = verbose)
    all.anchors <- Seurat::FindIntegrationAnchors(
        object.list = all.seurat.objects,
        anchor.features = anchor.features,
        reference = NULL,
        scale = scale,
        normalization.method = normalization.method,
        sct.clip.range = sct.clip.range,
        reduction = reduction,
        l2.norm = l2.norm,
        dims = dims,
        k.anchor = k.anchor,
        k.filter = k.filter,
        k.score = k.score,
        max.features = max.features,
        nn.method = nn.method,
        n.trees = n.trees,
        eps = eps,
        verbose = verbose)
    score <- NULL
    all.anchors <- all.anchors@anchors
    all.anchors$score <- round(x = all.anchors$score, digits = 2)
    colnames(all.anchors)[1:2] <- c('sample1', 'sample2')
    all.anchors$dataset1.name <- groups[all.anchors$dataset1]
    all.anchors$dataset2.name <- groups[all.anchors$dataset2]
    printColoredMessage(
        message = paste0(nrow(all.anchors), ' sample pairs (anchors) are found across the subgroups of the ', uv.variable, ' variable.'),
        color = 'blue',
        verbose = verbose)

    ## sanity check  ####
    printColoredMessage(message = '- Apply a sanity check on the anchors.',
                        color = 'blue',
                        verbose = verbose)
    sanity.check <- lapply(
        1:length(groups),
        function(x) {
            max.index <- max(all.anchors$sample1[all.anchors$dataset1 == x])
            sample.size <- sum(se.obj[[uv.variable]] == groups[x])
            if (max.index > sample.size) {
                stop('There are someting wrong with the order of anchros.')
            }
            max.index <- max(all.anchors$sample2[all.anchors$dataset2 == x])
            sample.size <-
                sum(se.obj[[uv.variable]] == groups[x])
            if (max.index > sample.size) {
                stop('There are someting wrong with the order of anchros.')
            }
        })

    ## add overall sample index  ####
    printColoredMessage(message = '- Add overall sample number to the anchors.',
                        color = 'blue',
                        verbose = verbose)
    for (x in 1:length(groups)) {
        index.a <- all.anchors$sample1[all.anchors$dataset1 == x]
        all.anchors$sample.index1[all.anchors$dataset1 == x] <- all.samples.index[se.obj[[uv.variable]] == groups[x]][index.a]
        index.b <- all.anchors$sample2[all.anchors$dataset2 == x]
        all.anchors$sample.index2[all.anchors$dataset2 == x] <- all.samples.index[se.obj[[uv.variable]] == groups[x]][index.b]
    }
    rm(all.seurat.objects)
    gc()

    ## all possible sets of PRPS ####
    printColoredMessage(message = '- Find all possible sets of PRPS.',
                        color = 'blue',
                        verbose = verbose)
    all.prps.sets <- split(
        x = all.anchors,
        f = all.anchors$sample.index1)
    all.prps.sets <- lapply(
        all.prps.sets,
        function(x) {
            temp.anchors <- x
            temp.anchors <-
                rbind(temp.anchors, do.call(
                    rbind,
                    lapply(temp.anchors$sample.index2, function(j)
                        all.anchors[all.anchors$sample.index2 == j,]
                        )))
            average.scores <- mean(temp.anchors$score)
            all.datasets <- sort(unique(c( temp.anchors$dataset1, temp.anchors$dataset2)))
            anchor.sets <-
                lapply(all.datasets, function(x) {
                    unique(c(
                        temp.anchors$sample.index1[temp.anchors$dataset1 == x],
                        temp.anchors$sample.index2[temp.anchors$dataset2 == x]
                    ))
                })
            anchor.datasets <-
                unlist(lapply(all.datasets, function(x) {
                    unique(c(
                        temp.anchors$dataset1.name[temp.anchors$dataset1 == x],
                        temp.anchors$dataset2.name[temp.anchors$dataset2 == x]
                    ))
                }))
            length.sets <- sapply(anchor.sets, length)
            return(
                list(average.scores = average.scores,
                    anchor.sets = setNames(anchor.sets, anchor.datasets),
                    length.sets = length.sets))
        })
    names(all.prps.sets) <- paste0('Anchor', levels(as.factor(all.anchors$sample.index1)))
    remove.anchor.sets <- sapply(
        all.prps.sets,
        function(x)
        sum(x$length.sets == 1) >= c(length(x$length.sets) - 1))
    all.prps.sets <- all.prps.sets[!remove.anchor.sets]
    printColoredMessage(
        message = paste0(length(all.prps.sets), ' possible PRPS stes are found.'),
        color = 'blue',
        verbose = verbose)

    ## check initial coverage ####
    printColoredMessage(message = '- Check the distribution of the PRPS sets across the batches.',
                        color = 'blue',
                        verbose = verbose)
    prps.coverage <-matrix(0, nrow = length(all.prps.sets), ncol = length(groups))
    colnames(prps.coverage) <- groups
    prps.coverage <- lapply(seq_along(all.prps.sets), function(i) {
        index <- match(names(prps.coverage[i, ]), names(all.prps.sets[[i]]$anchor.sets))
        all.prps.sets[[i]]$length.sets[index]
    })
    prps.coverage <- do.call(rbind, prps.coverage)
    colnames(prps.coverage) <- groups
    prps.coverage[is.na(prps.coverage)] <- 0

    if (sum(colSums(prps.coverage) == 0)) {
        printColoredMessage(
            message = paste(
                paste0(colnames(prps.coverage)[colSums(prps.coverage) == 0], collapse = ' & '), 'are not covered by any PRPS set.'),
            color = 'blue',
            verbose = verbose)
    }
    if (sum(rowSums(prps.coverage > 0) == length(groups)) > 0) {
        printColoredMessage(
            message = paste0('There are ',sum(rowSums(prps.coverage >= 0) == length(groups)),
                ' anchor sets across all subgroups of ', uv.variable,'.'),
            color = 'blue',
            verbose = verbose)

        printColoredMessage(
            message = paste0('There are ', sum(rowSums(prps.coverage >= min.prps.samples) > 1),
                ' PRPS sets with at least ', min.prps.samples,
                ' samples within at least two batches across all subgroups of ', uv.variable, '.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = paste0('There are ', sum(rowSums(prps.coverage >= min.prps.samples) == length(groups)),
                ' PRPS sets with at least ', min.prps.samples, ' samples within each batch across all subgroups of ',
                uv.variable, '.'),
            color = 'blue',
            verbose = verbose)
    } else {
        # check connection
        printColoredMessage(
            message = 'There is no any acnhor or PRPS sets that cover all batches. We assess the connection between differet sets:',
            color = 'blue',
            verbose = verbose)
        prps.connection <- lapply(
            1:nrow(prps.coverage),
            function(y) {
                batch.names.a <- names(which(prps.coverage[y, ] > 0))
                con.prps <- lapply(c(1:nrow(prps.coverage))[-y],
                           function(z) {
                               batch.names.b <- names(which(prps.coverage[z, ] > 0))
                               inter.samples <- intersect(batch.names.a, batch.names.b)
                               if (length(inter.samples) > 0) {
                                   sort(unique(c( batch.names.a, batch.names.b)), decreasing = FALSE)
                               } else sort(batch.names.a, decreasing = FALSE)
                           })
                all <- sort(unique(unlist(Filter(Negate(is.null), con.prps))), decreasing = FALSE)
                if (all.equal(all, groups))
                    break
            })
        all.covered.batches <- Filter(Negate(is.null), prps.connection)
        all.covered.batches <- unique(all.covered.batches)
        all.not.covered.batches <- unique(unlist(all.covered.batches))[unique(unlist(all.covered.batches)) %in% groups]
        if (length(unique(unlist(all.covered.batches))) == length(groups)) {
            printColoredMessage(
                message = 'The connections between PRPS sets can cover all the batches.',
                color = 'white',
                verbose = verbose)
        } else {
            printColoredMessage(
                message = 'All batches are not covered by PRPS.',
                color = 'red',
                verbose = verbose)
            printColoredMessage(
                message = paste0('The PRPS sets donot cover ', paste0(all.not.covered.batches, collapse = ' & '),' batches.'),
                color = 'white',
                verbose = verbose)
        }
    }

    ## filter anchors based on the scores ####
    if (!is.null(min.score)) {
        printColoredMessage(message = '- Filter the anchors based on the scores.',
                            color = 'blue',
                            verbose = verbose)
        if (isTRUE(verbose)) {
            ggplot(all.anchors, aes(x = score)) +
                geom_histogram() +
                ggtitle('The distribution of the anchor scores') +
                xlab('Frequency') +
                ylab('Scores') +
                theme(
                    panel.background = element_blank(),
                    plot.title = element_text(size = 12),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 12),
                    axis.title.y = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    legend.position = 'bottom'
                )
        }
        if (sum(all.anchors$score < min.score) == nrow(all.anchors)) {
            stop(
                paste0(
                    'All anchor pairs have the score less than "min.score":',
                    min.score,
                    '.'
                )
            )
        }
        printColoredMessage(
            message = paste0(
                sum(all.anchors$score < min.score),
                ' sample pairs (anchors) have the scores lower than min.score:',
                min.score,
                '. These anchors are removed.'
            ),
            color = 'blue',
            verbose = verbose
        )
        all.anchors <-
            all.anchors[all.anchors$score >= min.score, ]
    }

    ## filter PRPS sets based on sample size ####
    printColoredMessage(message = '- Filter the PRPS sets based on sample size.',
                        color = 'blue',
                        verbose = verbose)
    selected.anchors.sets <- unlist(lapply(names(all.prps.sets),
                                           function(x) {
                                               if (max.prps.samples != 'inf') {
                                                   index <- all.prps.sets[[x]]$length.sets >= min.prps.samples &
                                                       all.prps.sets[[x]]$length.sets <= max.prps.samples
                                               } else {
                                                   index <- all.prps.sets[[x]]$length.sets >= min.prps.samples
                                               }
                                               sum(index)
                                           }))
    selected.anchors.sets <- selected.anchors.sets > 1
    all.prps.sets <- all.prps.sets[selected.anchors.sets]

    printColoredMessage(
        message = paste0(
            'There are ',
            length(all.prps.sets),
            ' PRPS sets with at least ',
            min.prps.samples,
            ' samples and maximum = ',
            max.prps.samples,
            ' within at least two batches across all subgroups of ',
            uv.variable,
            '.'
        ),
        color = 'blue',
        verbose = verbose
    )

    # based on PRPS size
    score <- NULL
    selected.anchors.sets <- purrr::map_df(names(all.prps.sets), ~ {
        anchor.names <- sort(names(all.prps.sets[[.x]]$anchor.sets))
        data.frame(
            prps = paste0(anchor.names, collapse = '_'),
            score = all.prps.sets[[.x]]$average.scores,
            anchor = .x
        )
    })
    selected.anchors.sets <-
        unlist(lapply(unique(selected.anchors.sets$prps),
                      function(x) {
                          temp.data <-
                              selected.anchors.sets[selected.anchors.sets$prps == x, ]
                          if (nrow(temp.data) > max.prps.sets) {
                              temp.data <- temp.data[order(temp.data$score, decreasing = T), ]
                              temp.data$anchor[1:max.prps.sets]
                          } else{
                              temp.data$anchor
                          }
                      }))

    all.prps.sets <- all.prps.sets[selected.anchors.sets]

    # final check coverage ####
    prps.coverage <-
        matrix(0,
               nrow = length(all.prps.sets),
               ncol = length(groups))
    colnames(prps.coverage) <- groups
    prps.coverage <- lapply(seq_along(all.prps.sets), function(i) {
        index <-
            match(names(prps.coverage[i, ]),
                  names(all.prps.sets[[i]]$anchor.sets))
        all.prps.sets[[i]]$length.sets[index]
    })
    prps.coverage <- do.call(rbind, prps.coverage)
    colnames(prps.coverage) <- groups
    prps.coverage[is.na(prps.coverage)] <- 0

    if (sum(colSums(prps.coverage) == 0)) {
        printColoredMessage(
            message = paste(
                paste0(colnames(prps.coverage)[colSums(prps.coverage) == 0], collapse = ' & '),
                'are not covered by any PRPS set.'
            ),
            color = 'blue',
            verbose = verbose
        )
    }
    if (sum(rowSums(prps.coverage > 0) == length(groups)) > 0) {
        printColoredMessage(
            message = paste0(
                'There are ',
                sum(rowSums(prps.coverage >= 0) == length(groups)),
                ' anchor sets across all subgroups of ',
                uv.variable,
                '.'
            ),
            color = 'blue',
            verbose = verbose
        )

        printColoredMessage(
            message = paste0(
                'There are ',
                sum(rowSums(
                    prps.coverage >= min.prps.samples
                ) > 1),
                ' PRPS sets with at least ',
                min.prps.samples,
                ' samples within at least two batches across all subgroups of ',
                uv.variable,
                '.'
            ),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0(
                'There are ',
                sum(
                    rowSums(prps.coverage >= min.prps.samples) == length(groups)
                ),
                ' PRPS sets with at least ',
                min.prps.samples,
                ' samples within each batch across all subgroups of ',
                uv.variable,
                '.'
            ),
            color = 'blue',
            verbose = verbose
        )
    } else {
        # check connection
        printColoredMessage(message = 'There is no any acnhor or PRPS sets that cover all batches. We assess the connection between differet sets:',
                            color = 'blue',
                            verbose = verbose)
        prps.connection <- lapply(1:nrow(prps.coverage),
                                  function(y) {
                                      batch.names.a <- names(which(prps.coverage[y,] > 0))
                                      con.prps <-
                                          lapply(c(1:nrow(prps.coverage))[-y],
                                                 function(z) {
                                                     batch.names.b <- names(which(prps.coverage[z,] > 0))
                                                     inter.samples <-
                                                         intersect(batch.names.a, batch.names.b)
                                                     if (length(inter.samples) > 0) {
                                                         sort(unique(c(
                                                             batch.names.a, batch.names.b
                                                         )), decreasing = FALSE)
                                                     } else {
                                                         sort(batch.names.a, decreasing = FALSE)
                                                     }
                                                 })
                                      all <-
                                          sort(unique(unlist(Filter(
                                              Negate(is.null), con.prps
                                          ))), decreasing = FALSE)
                                      if (all.equal(all, groups))
                                          break
                                  })
        all.covered.batches <-
            Filter(Negate(is.null), prps.connection)
        all.covered.batches <- unique(all.covered.batches)
        all.not.covered.batches <-
            unique(unlist(all.covered.batches))[unique(unlist(all.covered.batches)) %in% groups]
        if (length(unique(unlist(all.covered.batches))) == length(groups)) {
            printColoredMessage(
                message = 'The connections between PRPS sets can cover all the batches.',
                color = 'white',
                verbose = verbose
            )
        } else {
            printColoredMessage(
                message = 'All batches are not covered by PRPS.',
                color = 'red',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    'The PRPS sets donot cover ',
                    paste0(all.not.covered.batches, collapse = ' & '),
                    ' batches.'
                ),
                color = 'white',
                verbose = verbose
            )
        }
    }

    ## create PRPS data ####
    # data transformation and normalization ####
    printColoredMessage(message = '-- Data transformation and normalization:',
                        color = 'magenta',
                        verbose = verbose)
    ## apply log ####
    if (isTRUE(apply.log) & !is.null(pseudo.count)) {
        printColoredMessage(
            message = paste0(
                'Applying log2 + ',
                pseudo.count,
                ' (pseudo.count) on the ',
                assay.name,
                ' data.'
            ),
            color = 'blue',
            verbose = verbose
        )
        expr.data <-
            log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else if (isTRUE(apply.log) & is.null(pseudo.count)) {
        printColoredMessage(
            message = paste0('Applying log2 on the ',
                             assay.name,
                             ' data.'),
            color = 'blue',
            verbose = verbose
        )
        expr.data <- log2(assay(x = se.obj, i = assay.name))
    } else if (isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0(
                'The ',
                assay.name,
                ' data will be used without any log transformation.'
            ),
            color = 'blue',
            verbose = verbose
        )
        expr.data <- assay(x = se.obj, i = assay.name)
    }
    prps.data <- lapply(names(all.prps.sets),
                        function(x) {
                            prps.sets <- all.prps.sets[[x]]$anchor.sets
                            temp.data <- sapply(1:length(prps.sets),
                                                function(y) {
                                                    rowMeans(expr.data[, prps.sets[[y]], drop = FALSE])
                                                })
                            colnames(temp.data) <-
                                rep(x = paste0(uv.variable, 'UnSupprps', x), length(prps.sets))
                            temp.data
                        })
    prps.data <- do.call(cbind, prps.data)
    if (sum(table(colnames(prps.data)) == 1)) {
        stop(
            'There are something wrong with the prps.data. All the column names of the prps.data are the same.'
        )
    }
    # saving the output ####
    out.put.name <- paste0(uv.variable, '|', assay.name)

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
        if (!'Supervised' %in% names(se.obj@metadata[['PRPS']])) {
            se.obj@metadata[['PRPS']][['supervised']] <- list()
        }
        ## Check if metadata PRPS already exist for supervised
        se.obj@metadata[['PRPS']][['UnSupervised']][[out.put.name]] <- prps.data

        printColoredMessage(message = '------------The createUnSupervisedPRPSbyAnchors function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else{
        printColoredMessage(message = '------------The createUnSupervisedPRPSbyAnchors function finished.',
                            color = 'white',
                            verbose = verbose)
        return(list(prps.data = prps.data))
    }
}
