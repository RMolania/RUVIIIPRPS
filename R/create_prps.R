
#' is used to create pseudo-replicates of pseudo samples of a SummarizedExperiment class object
#' using RUVIII-PRPS method.
#'
#' We will create distinct group of pseudo-replicates for each source of unwanted variation.
#' To correct for batch effect defined in the 'batch' argument, several group of pseudo-samples
#' will be created by averaging the samples of the same biological subtype in each batch. Then those
#' pseudo-samples will be defined as pseudo-replicates.
#' To correct for library size defined in the 'librarySize' argument, several group of pseudo-samples
#' will be created by averaging the top and bottom-ranked samples by library size of the same biological subtype in each batch.
#' Then those pseudo-samples will be defined as pseudo-replicates.
#' Similarly to correct for purity defined in the 'purity' argument, several group of pseudo-samples
#' will be created by averaging the top and bottom-ranked samples by purity of the same biological subtype in each batch.
#' Then those pseudo-samples will be defined as pseudo-replicates.
#'
#'
#' @param se A SummarizedExperiment object that will be used to create PRPS.
#' @param raw_data_assay_label String for selection of the name of the raw data assay of the SummarizedExperiment class object.
#' @param biology Vector containing the biological subtypes of the samples.
#' @param batch Vector containing the batch variable of the samples.
#' @param minSamplesPerBatchPS Minimum number of samples per batch to create a pseudo-sample.
#' @param librarySize Vector containing the library size of the samples.
#' @param include.ls Create PRPS for library size. By default is set to FALSE.
#' @param purity Vector containing the purity variable of the samples.
#' @param include.purity Create PRPS for purity. By default is set to FALSE.
#' @param minSamplesForPurityPS Minimum number of samples per batch to create a pseudo-sample to correct
#' for purity.
#' @param minSamplesForPurityPerBiology Minimum number of samples with the same biology to create a
#' pseudo-sample for purity.
#' @param minSamplesForLibrarySizePerBatch Minimum number of samples per batch to create pseudo-samples
#' to correct for library size
#' @param minSamplesForLibrarySizePS Minimum number of samples with similar library size to create a
#' pseudo-sample to correct for library size.
#'
#' @return list List containing the PRPS created to correct for batch, for library size and for purity.

#' @importFrom Matrix rowMeans
#' @importFrom SummarizedExperiment assay colData
#' @export

create_prps <- function(
        se,
        raw_data_assay_label,
        librarySize,
        biology,
        batch,
        purity,
        include.ls = FALSE,
        include.purity = FALSE,
        minSamplesPerBatchPS = 3,
        minSamplesForPurityPerBiology = 12,
        minSamplesForPurityPS = 3,
        minSamplesForLibrarySizePerBatch = 10,
        minSamplesForLibrarySizePS = 3
){
    ### Check se
    if (!class(se)[1] == 'SummarizedExperiment') {
        stop('Please provide a summarized experiment object.\n')
    }
    ### Check minSamples
    if(include.purity & minSamplesForPurityPS > minSamplesForPurityPerBiology){
        stop('error: minSamplesForPurityPS can not be smaller than minSamplesForPurityPerBiology')
    } else if(include.purity & minSamplesForPurityPerBiology < 2*minSamplesForPurityPS){
        stop('error: minSamplesForPurityPerBiology should be at least two times larger than minSamplesForPurityPS')
    } else if(include.ls & minSamplesForLibrarySizePS > minSamplesForLibrarySizePerBatch) {
        stop('error: minSamplesForLibrarySizePerBatch can not be smaller than minSamplesForLibrarySizePS')
    } else if(include.ls & minSamplesForLibrarySizePerBatch < 2*minSamplesForLibrarySizePS ){
        stop('error: minSamplesForLibrarySizePerBatch should be at least two times larger than minSamplesForLibrarySizePS')
    }
    ### Biology
    expr.data = as.data.frame(assay(se, raw_data_assay_label))
    sample.info = droplevels(as.data.frame(colData(se)))
    row.names(sample.info) <- colnames(expr.data)
    sample.info$biology <- apply(
        sample.info[ , biology, drop = FALSE],
        1,
        paste,
        collapse = "-"
    )
    ### Biology - Batch
    sample.info$biology.batch <- apply(
        sample.info[, c(biology, batch)],
        1,
        paste,
        collapse = "_"
    )
    ### removing batch effects
    # create PS per biology/batch
    selected.biology.ps.batch <- unlist(lapply(
        unique(sample.info$biology),
        function(x){
            index <- sample.info$biology == x
            if(sum( table(sample.info$biology.batch[index] ) >= minSamplesPerBatchPS) > 1 ){
                x
            }
        }))
    if(length(selected.biology.ps.batch) > 0){
        message('PRPS are generated for batch effects')
    }else{
        message('error: there are not enough samples to create pseudo-samples for batch effects removal, you may want to lower minSamplesPerBatchPS')
    }
    sample.info.ps.batch <- sample.info[sample.info$biology %in% selected.biology.ps.batch , ]
    expr.data.ps.batch <- expr.data[, row.names(sample.info.ps.batch)]
    ### sort samples
    selected.batches <- names(which(table(sample.info.ps.batch$biology.batch) >= minSamplesPerBatchPS))
    ps.batch <- sapply(
        selected.batches,
        function(x) {
            index <- sample.info.ps.batch$biology.batch == x
            rowMeans(expr.data.ps.batch[, index])
        })

    if(include.ls){
        selected.batches.ls <- names(
            which(table(sample.info$biology.batch) >= minSamplesForLibrarySizePerBatch)
        )
        if(length(selected.batches.ls) > 0){
            message('PRPS are generated for library size effects')
            sample.info <- sample.info[
                with(sample.info,
                     order(sample.info[, 'biology.batch'],
                           sample.info[, librarySize])), ]
            expr.data <- expr.data[, row.names(sample.info)]
            ps.ls <- lapply(
                selected.batches.ls,
                function(x){
                    index <- sample.info$biology.batch == x
                    ls.data <- expr.data[ , index]
                    low.ls <- Matrix::rowMeans(ls.data[ , 1:minSamplesForLibrarySizePS])
                    high.ls <- rowMeans(ls.data[ , c(ncol(ls.data)-(minSamplesForLibrarySizePS - 1)):ncol(ls.data) ])
                    all <- cbind(low.ls, high.ls)
                    colnames(all) <- rep(paste(x, 'LS', sep = '-'), 2)
                    all
                })
            ps.ls <- do.call(cbind, ps.ls)

        }else{
            message('error: there are not enough samples to create pseudo-samples for removal of library size effects,
              you may want to lower minSamplesForLibrarySizePerBatch')
        }
    }else if (! include.ls){
        print('PRPS is not generated for librray size effects')
        ps.ls = list()
    }
    if(include.purity ){
        selected.biology.purity <- names(
            which(table(sample.info$biology) >= minSamplesForPurityPerBiology)
        )
        if(length(selected.biology.purity) > 0){
            message('PRPS are generated for purity effects')
            sample.info <- sample.info[
                with(sample.info,
                     order(sample.info[, 'biology'],
                     #order(sample.info[, 'biology.batch'],
                           sample.info[, purity])),]
            expr.data <- expr.data[, row.names(sample.info)]
            ps.purity <- lapply(
                selected.biology.purity,
                function(x) {
                    index <- sample.info$biology == x
                    purity.data <- expr.data[, index]
                    low.pur <- rowMeans(purity.data[, 1:minSamplesForPurityPS])
                    high.pur <- rowMeans(purity.data[, c(ncol(purity.data) - (minSamplesForPurityPS - 1)):ncol(purity.data)])
                    all <- cbind(low.pur, high.pur)
                    colnames(all) <- rep(paste(x, 'purity', sep = '-'), 2)
                    all
                })
            ps.purity <- do.call(cbind, ps.purity)
        }else{
            message('error: there are not enough samples to make pseudo-samples for purity variation,
              you may want to lower minSamplesForPurityPerBiology')
        }
    } else if (!include.purity){
        print('PRPS is not generated for purity effects')
        ps.purity = list()
    }
    return(list(ps.batch = ps.batch, ps.ls = ps.ls, ps.purity = ps.purity))
}
