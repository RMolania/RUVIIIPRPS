#' is used to perform differential gene expression analysis of the assays in a SummarizedExperiment object.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Optional string or list of strings for the selection of the name(s) of the assay(s) of the
#' SummarizedExperiment class object to compute the ANOVA. By default all the assays of the SummarizedExperiment object
#' will be selected.
#' @param variable String of the label of a categorical variable such as sample types or batches from colData(se.obj).
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data. By default
#' the log transformation will be selected.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na TO BE DEFINED.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param plot.output Logical. Indicates whether to plot the boxplot of the F-test statistics, by default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#'
#' @return SummarizedExperiment A SummarizedExperiment object containing the log2 F-statistics of ANOVA on the continuous variable
#' and if requested the associated boxplot.


#' @importFrom SummarizedExperiment assays assay
#' @importFrom matrixTests row_wilcoxon_twosample
#' @export

genesDEA <- function(
        se.obj,
        assay.names = 'All',
        variable,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        plot.output = TRUE,
        verbose = TRUE
){
    printColoredMessage(message = '------------The genesDEA function starts:',
                        color = 'white',
                        verbose = verbose)
    # check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be null.')
    }
    if(length(assay.names) == 1 & assay.names!= 'All'){
        if(!assay.names %in% names(assays(se.obj)) )
            stop('The assay name cannot be found in the SummarizedExperiment object.')
    }
    if(length(assay.names) > 1){
        if(sum(!assay.names %in% names(assays(se.obj))) > 0 )
            stop('The assay names cannot be found in the SummarizedExperiment object.')
    }
    if (is.null(variable)) {
        stop('The variable cannot be empty.')
    } else if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop(paste0('The ', variable, ' should be a categorical variable.'))
    } else if (length(unique(se.obj@colData[, variable])) < 2) {
        stop(paste0('The ', variable, ', contains only one level. ANOVA cannot be performed.'))
    }
    if (apply.log){
        if (pseudo.count < 0)
            stop('The value of pseudo.count cannot be negative.')
    }

    # assays ####
    if (length(assay.names) == 1 && assay.names == 'All') {
        assay.names <- as.factor(names(assays(se.obj)))
    } else assay.names <- as.factor(unlist(assay.names))

    # assess the SummarizedExperiment object ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variable,
            remove.na = remove.na,
            verbose = verbose)
    }
    all.contrasts <- combn(
        x = unique(colData(se.obj)[[variable]]),
        m = 2)
    # anova ####
    all.wilcoxon.test <- lapply(
        levels(assay.names),
        function(x){
            # data transformation ####
            printColoredMessage(message = '-- Data transformation:',
                                color = 'magenta',
                                verbose = verbose)
            if (apply.log) {
                printColoredMessage(message = paste0(
                    'Perform log2 + ',
                    pseudo.count,
                    '(pseudo.count) on the ',
                    x,
                    ' data.'),
                    color = 'blue',
                    verbose = verbose)
                temp.data <- log2(assay(x = se.obj, i = x) +  pseudo.count)
            } else{
                printColoredMessage(
                    message = paste0('Please make sure that the ', x, ' data is log transformed.'),
                    color = 'blue',
                    verbose = verbose)
                temp.data <- assay(x = se.obj, i = x)
            }
            printColoredMessage(
                message = paste0('-- Perform Wilcoxon test.'),
                color = 'blue',
                verbose = verbose)
            de.results <- lapply(
                1:ncol(all.contrasts),
                function(i){
                    data1 <- temp.data[ , colData(se.obj)[['Call']] == all.contrasts[1 , i] ]
                    data2 <- temp.data[ , colData(se.obj)[['Call']] == all.contrasts[2 , i] ]
                    de.table <- row_wilcoxon_twosample(data1, data2)[ , c('obs.x', 'obs.y', 'pvalue')]
                })
            names(de.results) <- sapply(
                1:ncol(all.contrasts),
                function (x)
                    paste(all.contrasts[1 , x], all.contrasts[2 , x], sep = '&'))
        })
    names(all.wilcoxon.test) <- levels(assay.names)
    # save the results ####
    printColoredMessage(
        message = '-- Save the Wilcoxon test results:',
        color = 'magenta',
        verbose = verbose)
    ## add results to the SummarizedExperiment object ####
    if (save.se.obj == TRUE) {
        for (x in levels(assay.names)) {
            ## check if metadata metric already exist
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!'metric' %in% names(se.obj@metadata)) {
                se.obj@metadata[['metric']] <- list()
            }
            ## check if metadata metric already exist for this assay
            if (!x %in% names(se.obj@metadata[['metric']])) {
                se.obj@metadata[['metric']][[x]] <- list()
            }
            ## check if metadata metric already exist for this assay and this metric
            if (!paste0('gene.wilcoxon') %in% names(se.obj@metadata[['metric']][[x]])) {
                se.obj@metadata[['metric']][[x]][['gene.wilcoxon']] <- list()
            }
            ## check if metadata metric already exist for this assay, this metric and this variable
            se.obj@metadata[['metric']][[x]][['gene.wilcoxon']][[variable]] <- all.wilcoxon.test[[x]]
        }
        printColoredMessage(
            message = 'The Wilcoxon results for indiviaul assay are saved to metadata@metric',
            color = 'blue',
            verbose = verbose)

        ## Plot and save the plot into se.obj@metadata$plot
        if (plot.output == TRUE) {
            printColoredMessage(
                message = '-- Plot the Wilcoxon results:',
                color = 'magenta',
                verbose = verbose
            )
            printColoredMessage(
                message = 'A histogram of pvalues of the Wilcoxon results are saved to metadata@plot.',
                color = 'blue',
                verbose = verbose)
            se.obj <- plotMetric(
                se.obj,
                assay.names = assay.names,
                metric = 'gene.wilcoxon',
                variable = variable,
                verbose = verbose)
        }
        printColoredMessage(message = '------------The genesDEA function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
        ## return the results as a list ####
    } else if (save.se.obj == FALSE) {
        printColoredMessage(
            message = 'The Wilcoxon results for indiviaul assay are saved as list.',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The genesDEA function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.wilcoxon.test = all.wilcoxon.test)
    }

}
