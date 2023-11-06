#' is used to create PRPS for a categorical variable as source of unwanted variation.
#'
#'
#' @param se.obj A summarized experiment object.
#' @param assay.name A name of the assays in summarized experiment object.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data.
#' @param pseudo.count Numeric. A pseudo count to add to each gene before applying log.
#' @param bio.variable String of the label of a categorical variable that specifies major biological groups
#' such as samples types from colData(se).
#' @param uv.variable String of the label of a continuous or categorical variable.
#' such as samples types, batch or library size from colData(se) that will be used to define PRPS.
#' @param min.sample.prps Numeric. The minimum number of homogeneous biological groups to create pseudo-sample.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' By default it is set to TRUE.
#' @param remove.na TO BE DEFINED.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class
#' object 'se.obj' or to output the result, by default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution of the functions, by default it is set to TRUE.

#' @return SummarizedExperiment A SummarizedExperiment object containing the PRPS data or just PRPS data.
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom knitr kable
#' @export

prpsForCategoricalUV <- function(se.obj,
                                 assay.name,
                                 uv.variable,
                                 bio.variable,
                                 min.sample.prps = 3,
                                 apply.log = TRUE,
                                 pseudo.count = 1,
                                 remove.na = 'both',
                                 assess.se.obj = TRUE,
                                 save.se.obj = TRUE,
                                 verbose = TRUE) {
    ## Check the input
    printColoredMessage(message = '------------The prpsForCategoricalUV function starts:',
                        color = 'white',
                        verbose = verbose)
    # assess se obj ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(uv.variable, bio.variable),
            remove.na = remove.na)
    }
    # checking the input of the function ####
    if (length(assay.name) > 1) {
        stop('The function can only take a single assay.name.')
    } else if (is.null(assay)) {
        stop('The assay.name cannot be empty.')
    } else if (length(uv.variable) > 1) {
        stop('The function can only take a single categorical variable for the uv.variable argument.')
    } else if (is.null(uv.variable)) {
        stop('The uv.variable cannot be empty.')
    } else if (!class(se.obj[[uv.variable]]) %in% c('character', 'factor')) {
        stop('The uv.variable should be a categorical source of unwanted variation, e.g., time effects')
    } else if (length(unique(se.obj[[uv.variable]])) == 1) {
        stop('The uv.variable should have at least two levels.')
    } else if(!is.null(bio.variable)){
        if (length(bio.variable) > 1) {
            stop('The function can only take a single biological variable.')
        } else if (!class(se.obj[[bio.variable]]) %in% c('character', 'factor')) {
            stop('The bio.variable should be a categorical biological variable, e.g., cancer subtypes')
        } else if (length(unique(se.obj[[bio.variable]])) == 1) {
            stop('The bio.variable should have at least two groups.')
        }
    } else if (min.sample.prps <= 1) {
        stop('The minimum value for the min.sample.prps is 2.')
    }
    ### replace "_" with "-"
    if(!is.null(bio.variable)){
        se.obj[[bio.variable]] <- gsub('_', '-', se.obj[[bio.variable]])
    }
    se.obj[[uv.variable]] <- gsub('_', '-', se.obj[[uv.variable]])
    # data transformation ####
    if (apply.log) {
        printColoredMessage(message = 'Applying log2 transformation on the data before creating PRPS.',
                            color = 'blue',
                            verbose = verbose)
        expre.data <- log2(assay(se.obj, assay.name) + pseudo.count)
    } else{
        printColoredMessage(message = 'The assay data will be used for PRPS without any transformation.',
                            color = 'blue',
                            verbose = verbose)
        expre.data <- assay(se.obj, assay.name)
    }

    #if(!is.null(bio.variable)){
        ### Table of biological variable and unwanted variable
        bio.batch.table <- table(
            colData(se.obj)[[bio.variable]],
            colData(se.obj)[[uv.variable]]
            )
        bio.dist <- rowSums(bio.batch.table >= min.sample.prps)
        bio.batch.table <- bio.batch.table[bio.dist > 1 , ]
        bio.batch.table.to.plot <- bio.batch.table
        if(nrow(bio.batch.table) == 0){
            stop('There are not enough samples to create PRPS across the batches.')
        }
        ### create PRPS across uv_variables
        printColoredMessage(
            message = paste0(
                '### Creating PRPS across all batches (ones that have at least ',
                min.sample.prps,
                ' samples of a homogenous biological population) of ',
                uv.variable,
                '.'
            ),
            color = 'magenta',
            verbose = verbose
        )
        prps.sets <- lapply(
            1:nrow(bio.batch.table),
            function(y) {
                selected.cat.var <- colnames(bio.batch.table)[bio.batch.table[y , ] >= min.sample.prps]
                ps.matrix <- sapply(
                    selected.cat.var,
                    function(z) {
                        index.sample <- colData(se.obj)[[bio.variable]] == row.names(bio.batch.table)[y] &
                            colData(se.obj)[[uv.variable]] == z
                        rowMeans(expre.data[, index.sample])
                    })
                colnames(ps.matrix) <- paste(uv.variable, row.names(bio.batch.table)[y] ,colnames(ps.matrix),sep = '||')
                ps.matrix
            })
        prps.sets <- do.call(cbind, prps.sets)

        # ######### THIS PARt NEED to BE tEStED ############
    # } else {
    #     batches <- findRepeatingPatterns(
    #         vector = colData(se.obj)[[uv.variable]],
    #         n = min.sample.prps
    #         )
    #     if(length(batches) > 1){
    #         prps.sets <- sapply(
    #             batches,
    #             function(x){
    #                 index.sample <- colData(se.obj)[[uv.variable]] == x
    #                 if(sum(index.sample) >= min.sample.prps){
    #                     ps.matrix <- as.matrix(rowMeans(expre.data[, index.sample]))
    #                     colnames(ps.matrix) <- paste('batch', x , sep = '||')
    #                     return(ps.matrix)
    #                 }
    #             })
    #     } else{
    #         stop('There are not enough samples to create PRPS across the batches.')
    #     }
    # }
    printColoredMessage(
        message = paste0(
            length(unique(colnames(
                prps.sets
            ))),
            ' sets of PRPS with the total number of ',
            ncol(prps.sets),
            ' pseudo-samples are created between different batches of the ',
            uv.variable,
            ' variable.'
        ),
        color = 'blue',
        verbose = verbose
    )
    if(verbose){
        bio.batch.table.to.plot[bio.batch.table.to.plot >= min.sample.prps] <- 'PS'
        bio.batch.table.to.plot[bio.batch.table.to.plot < min.sample.prps] <- 'No'
        print(kable(bio.batch.table.to.plot))
    }

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
        }
        se.obj@metadata[['PRPS']][['supervised']][[paste0('bio:',
                                                              bio.variable,
                                                              '||',
                                                              'uv:',
                                                              uv.variable,
                                                              '||',
                                                              'data:',
                                                              assay.name)]] <- prps.sets


        printColoredMessage(message= paste0(
            'The PRPS are saved to metadata@PRPS$supervised',
            paste0('$bio:', bio.variable,'||','uv:',uv.variable,'||','data:',assay.name),
            '.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(message = '------------The prpsForCategoricalUV function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    } else{
        printColoredMessage(message = '------------The prpsForCategoricalUV function finished.',
                            color = 'white',
                            verbose = verbose)
        return(prps.sets)
    }

}



