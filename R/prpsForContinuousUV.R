#' is used to create PRPS for a continuous variable as source of unwanted variation.
#'
#'
#' @param se.obj A summarized experiment object.
#' @param assay.name A name of the assays in summarized experiment object.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data.
#' @param pseudo.count Numeric. A pseudo count to add to each gene before applying log.
#' @param bio.variable String of the label of a categorical variable that specifies major biological groups
#' such as samples types from colData(se).
#' @param uv.variable String of the label of a continuous or categorical variable.
#' @param min.sample.prps Numeric. The minimum number of homogeneous biological groups to create pseudo-sample.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' By default it is set to TRUE.
#' @param remove.na TO BE DEFINED.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result, by default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or
#' messages displayed during the execution of the functions, by default it is set to TRUE.

#' @return SummarizedExperiment A SummarizedExperiment object containing the PRPS data or just PRPS data.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr group_by arrange slice desc
#' @importFrom tidyr %>%
#' @import ggplot2
#' @export

prpsForContinuousUV <- function(se.obj,
                                assay.name,
                                apply.log = TRUE,
                                pseudo.count = 1,
                                uv.variable,
                                bio.variable,
                                min.sample.prps = 3,
                                assess.se.obj = TRUE,
                                remove.na = 'both',
                                save.se.obj = TRUE,
                                verbose = TRUE) {
    printColoredMessage(message = '------------The prpsForContinuousUV function starts.',
                        color = 'white',
                        verbose = verbose)
    ### Checking the function inputs
    if (length(assay.name) > 1) {
        stop('The function can only take a single assay name.')
    } else if (length(uv.variable) > 1) {
        stop('The function can only take a single uv.variable variable.')
    } else if (length(bio.variable) > 1) {
        stop('The function can only take a single bio.variable.')
    } else if (min.sample.prps < 1) {
        stop('The minimum value for the min_sample_prps is 1.')
    } else if(var(se.obj[[uv.variable]]) == 0 ){
        stop(paste('The variable ', uv.variable, ' has no variation. Then, there is no need to create PRPS.'))
    } else if(!class(se.obj[[uv.variable]]) %in% c('numeric', 'integer')){
        stop(paste('The variable ', uv.variable, ' should be a numeric or integer variable.'))
    } else if(!class(se.obj[[bio.variable]]) %in% c('factor', 'character')){
        stop(paste('The variable ', bio.variable, ' should be a factor or character variable.'))
    }
    # check the SummarizedExperiment ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(uv.variable, bio.variable),
            remove.na = remove.na,
            verbose = verbose
        )
    }
    # finding possible biological groups for PRPS ####
    printColoredMessage(
        message = paste0(
            '### Finding biological groups with at least (2* min.sample.prps) ',
            2 * min.sample.prps,
            ' samples.'
        ),
        color = 'magenta',
        verbose = verbose
    )
    bio.cont.prps <- findRepeatingPatterns(vector = colData(se.obj)[[bio.variable]],
                                           n = 2 * min.sample.prps)
    if (length(bio.cont.prps) == 1) {
        printColoredMessage(
            message = paste0(
                'There are a ',
                length(bio.cont.prps) ,
                ' homogeneous biological group that contains at least (2* min.sample.prps) ',
                2 * min.sample.prps,
                ' samples.'
            ),
            color = 'blue',
            verbose = verbose
        )
    } else if (length(bio.cont.prps) > 0) {
        printColoredMessage(
            message = paste0(
                'There are ',
                length(bio.cont.prps) ,
                ' homogeneous biological groups that each contains at least (2*min.sample.prps) ',
                2 * min.sample.prps,
                ' samples.'
            ),
            color = 'white',
            verbose = verbose
        )
    } else{
        stop(
            'There is not enough samples to create PRPS for contenious sources of unwanted variation.'
        )
    }
    printColoredMessage(
        message = '### Data preparation before creating PRPS.',
        color = 'magenta',
        verbose = verbose
    )
    # data transformation ####
    if (apply.log) {
        printColoredMessage(
            message = 'Applying log2 transformation on the data before creating PRPS.',
            color = 'blue',
            verbose = verbose
        )
        expre.data <- log2(assay(se.obj, assay.name) + pseudo.count)
    } else{
        printColoredMessage(
            message = 'The assay data will be used wihtout any transformation.',
            color = 'blue',
            verbose = verbose
        )
        expre.data <- assay(se.obj, assay.name)
    }
    # creating PRPS ####
    printColoredMessage(message = '### Creating a PRPS set with two PS for each individual homogeneous biological group.',
                        color = 'magenta',
                        verbose = verbose)
    se.obj <- se.obj[, se.obj[[bio.variable]] %in% bio.cont.prps]
    se.obj$sOrder <- c(1:ncol(se.obj))
    sample.annot <- as.data.frame(colData(se.obj))
    sample.annot.temp <-
        sample.annot[, c('sOrder', uv.variable, bio.variable)]
    top <- sample.annot.temp %>%
        arrange(desc(!!sym(uv.variable))) %>%
        group_by(!!sym(bio.variable)) %>%
        slice(1:min.sample.prps)
    bot <- sample.annot.temp %>%
        arrange(!!sym(uv.variable)) %>%
        group_by(!!sym(bio.variable)) %>%
        slice(1:min.sample.prps)
    prps.sets <- vector('list', length = ceiling(nrow(bot) / min.sample.prps))
    prps.sets <- lapply(
        seq(1, nrow(bot), min.sample.prps),
        function(y) {
            index <- y:((min.sample.prps + y) - 1)
            ps.all <- cbind(rowMeans(expre.data[, top[['sOrder']][index]]),
                      rowMeans(expre.data[, bot[['sOrder']][index]]))
            colnames(ps.all) <- rep(paste0(
                                           uv.variable,
                                           '||',
                                           unique(top[[bio.variable]][y])), 2)
            ps.all
        })
    prps.sets <- do.call(cbind, prps.sets)
    printColoredMessage(
        message = paste0(
            'In total ',
            length(unique(colnames(prps.sets))),
            ' PRPS sets with ',
            ncol(prps.sets),
            ' the total number of pseudo-samples are created for the removal of  ',
            uv.variable,
            ' effects.'
        ),
        color = 'blue',
        verbose = verbose
    )
    # saving the output ####
    if (save.se.obj) {
        ## Check if metadata PRPS already exists
        if(length(se.obj@metadata)==0 ) {
            se.obj@metadata[['PRPS']] <- list()
        }
        ## Check if metadata PRPS already exists
        if(!'PRPS' %in% names(se.obj@metadata) ) {
            se.obj@metadata[['PRPS']] <- list()
        }
        ## Check if metadata PRPS already exist for supervised
        if(!'supervised' %in% names(se.obj@metadata[['PRPS']]) ) {
            se.obj@metadata[['PRPS']][['supervised']] <- list()
        }

        ## Check if metadata PRPS already exist for supervised
        if(!paste0('bio:', bio.variable,'||','uv:',uv.variable,'||','data:',assay.name) %in% names(se.obj@metadata[['PRPS']][['supervised']])) {
            se.obj@metadata[['PRPS']][['supervised']][[paste0('bio:', bio.variable,'||','uv:',uv.variable,'||','data:',assay.name)]]<- list()
        } else {
            se.obj@metadata[['PRPS']][['supervised']][[paste0('bio:',
                                                          bio.variable,
                                                          '||',
                                                          'uv:',
                                                          uv.variable,
                                                          '||',
                                                          'data:',
                                                          assay.name)]] <- prps.sets
        }

    printColoredMessage(message= paste0(
        'The PRPS are saved to metadata@PRPS$supervised',
        paste0('$bio:', bio.variable,'||','uv:',uv.variable,'||','data:',assay.name),
        '.'),
        color = 'blue',
        verbose = verbose)
    return(se.obj)

    }else{
        return(prps.sets)
    }
    printColoredMessage(message = '------------The prpsForContinuousUV function finished.',
                        color = 'white',
                        verbose = verbose)
}