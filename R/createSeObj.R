#' is used to create a SummarizedExperiment object.

#' @description
#' This function creates a SummarizedExperiment object from tabular expression data and sample annotation. In addition,
#' the function can remove lowly expressed genes, add a range of annotations for genes, and provide several sets of
#' housekeeping genes and a immunStroma gene signature.

#' @param assays A list of assays or expression data. The genes should be in row and samples in the column. The row names
#' of the assays should be in the same order.
#' @param remove.lowly.expressed.genes 	Logical. If TRUE the function removes lowly expressed genes from the assay that
#' is provide in the raw.count.assay.name argument.
#' @param raw.count.assay.name Symbol. The name of raw counts data in the assays.
#' @param count.cutoff Numeric. Minimum count required for at least some sample groups. If the "biological.group" argument
#' is equal to NULL, all samples will be considered as one group. Otherwise, the smallest subgroups of the "biological.group"
#' will be considered.
#' @param biological.group Symbol. Indicates a column name in the sample annotation that specifies the biological groups.
#' The smallest population is considered for removing lowly expressed genes.
#' @param minimum.proportion Numeric. In large sample situations, the minimum proportion of samples in a group that a
#' gene needs to be expressed in.
#' @param calculate.library.size Logical. If TRUE then library size is calculated using the raw.count.assay.name. The
#' library size should be calculated after removing lowly expressed genes.
#' @param sample.annotation A data frame, contains information for individual samples in the assay(s). The order of row
#' names of the sample annotation should be the same as in the assay(s).
#' @param create.sample.annotation Logical. If TRUE then a sample annotation the initially contains the column names of
#' the assay(s) will be created.
#' @param gene.annotation A data frame, that contains details e.g. chromosome names,GC content,... for individual genes.
#' @param create.gene.annotation Logical. If TRUE then a gene annotation that initially contains row names of the assay(s)
#' will be created.
#' @param add.gene.details Logical. If TRUE then a pr-set or provided gene details in the gene.details argument will be
#' added to the gene annotation.
#' @param gene.group A name of a gene class in the row names of the assay(s). This must be one of the 'entrezgene_id',
#' 'hgnc_symbol', 'ensembl_gene_id'.
#' @param gene.details A vector of gene details to be added to the gene annotation.
#' @param add.housekeeping.genes Logical. if TRUE then several sets of publicly available "housekeeping" will be added to
#' the gene annotation. The housekeeping could be potentially used as negative control genes for the RUV normalization.
#' @param add.immunStroma.genes Logical. If TRUE, the immune and stromal genes signature from Kosuke Yoshihara et.al will
#' be added to the gene annotation. These gene signatures, can be used to estimate tumor purity.
#' @param metaData Any metadata data. The metadata can be in any format and dimensions.
#' @param verbose Logical. If TRUE shows the messages.

#' @return A summarizedExperiment that contains assays, gene annotation, samples annotation and metadata.

#' @author Ramyar Molania

#' @importFrom SummarizedExperiment SummarizedExperiment assay rowData
#' @importFrom Matrix colSums
#' @importFrom dplyr left_join
#' @importFrom biomaRt getBM useMart useDataset
#' @importFrom S4Vectors DataFrame
#' @importFrom knitr kable
#' @export

createSeObj <- function(
        assays ,
        raw.count.assay.name = NULL,
        remove.lowly.expressed.genes = FALSE,
        count.cutoff = 10,
        biological.group = NULL,
        minimum.proportion = .5,
        calculate.library.size = FALSE,
        sample.annotation = NULL,
        create.sample.annotation = FALSE,
        gene.annotation = NULL,
        create.gene.annotation = FALSE,
        add.gene.details = FALSE,
        gene.group = NULL,
        gene.details = NULL,
        add.housekeeping.genes = FALSE,
        add.immunStroma.genes = FALSE,
        metaData = NULL,
        verbose = TRUE)
{
    printColoredMessage(message = '------------The createSeObj function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if(!inherits(assays, what = 'list')){
        stop('The assays should be a list of expression datastes (genes in rows and samples in columns.).')
    }
    ## dimensions and orders of the assays ####
    if(length(assays) > 1){
        dim.data <- unlist(lapply(
            names(assays),
            function(x) c(nrow(assays[[x]]), ncol(assays[[x]])))
        )
        if(length(unique(dim.data)) !=2 ){
            stop('The dimensions of the assays in the assays must be the same.')
        }
        all.assays <- combn(1:length(assays), 2)
        m.out <- lapply(
            1:ncol(all.assays),
            function(x){
                if(!all.equal(row.names(assays[[all.assays[ 1, x]]]),row.names(assays[[all.assays[ 2, x]]]))){
                    stop('The row names of the assays should be in the same order.')
                }
            })
        m.out <- lapply(
            1:ncol(all.assays),
            function(x){
                if(!all.equal(colnames(assays[[all.assays[ 1, x]]]),colnames(assays[[all.assays[ 2, x]]]))){
                    stop('The column names of the assays should be in the same order.')
                }
            })
    }
    ## check raw.count.assay.name ####
    if(!is.null(raw.count.assay.name)){
        if(length(raw.count.assay.name) > 1){
            stop('The raw.count.assay.name should contain a single assay name.')
        }
        if(!raw.count.assay.name %in% names(assays)){
            stop('The raw.count.assay.name should be in the assay list.')
        }
    }
    ## remove lowly expressed genes ####
    if(remove.lowly.expressed.genes){
        if(is.null(raw.count.assay.name)){
            stop('To remove lowly expressed genes, the raw.count.assay.name should be provided.')
        }
        if(!is.null(biological.group)){
            if(is.null(sample.annotation)){
                stop('To find biological.group, please provide a sample.annotation that contains the biological.group variable.')
            }
            if(!biological.group %in% colnames(sample.annotation)){
                stop('The biological.group cannot be found in the sample.annotation.')
            }
        }
        if(minimum.proportion > 1 | minimum.proportion < 0){
            stop('The minimum.proportion should between 0 to 1.')
        }
        if(is.null(count.cutoff)){
            stop('The count.cutoff cannot be empty.')
        }
        if(count.cutoff < 0){
            stop('The value of the count.cutoff cannot be negative.')
        }
    }
    ## library size calculation ####
    if(is.null(raw.count.assay.name) & calculate.library.size){
        stop('To calculate library size, a single raw count assay name should be provided.')
    }
    ## sample annotation ####
    if(!is.null(sample.annotation)){
        if(nrow(sample.annotation) != ncol(assays[[1]])){
            stop('The number of rows in sample.annotation and the number of columns in the assays should be the same.')
        } else if (!all.equal(row.names(sample.annotation), colnames(assays[[1]])) ){
            stop('The order and lables of the row names of sample.annotation and colummn names of the assays should be the same.')
        }
        if(create.sample.annotation){
            stop('A sample.annotation has been provided, then create.sample.annotation must be FALSE.')
        }
    }
    ## gene annotation ####
    if(!is.null(gene.annotation) & create.gene.annotation){
        stop('A gene annotattion is provided, then create.gene.annotation must be FALSE.')
    }
    if(add.gene.details & is.null(gene.annotation) & !create.gene.annotation){
        stop('To add add.gene.details , the create.gene.annotation should be TRUE.')
    }
    if(add.housekeeping.genes & is.null(gene.annotation) & !create.gene.annotation){
        stop('To add add.housekeeping.genes a gene annotattion should be provided or create.gene.annotation should be TRUE.')
    }
    if(add.immunStroma.genes & is.null(gene.annotation) & !create.gene.annotation){
        stop('To add add.immunStroma.genes a gene annotattion should be provided or create.gene.annotation should be TRUE.')
    }
    if(add.gene.details & is.null(gene.group)){
        stop('To add gene details, the gene.group should be specified (entrezgene_id, hgnc_symbol and ensembl_gene_id).')
    }
    if(add.gene.details &!is.null(gene.group)){
        if(!gene.group %in% c('entrezgene_id', 'hgnc_symbol', 'ensembl_gene_id')){
            stop('The gene group should be: entrezgene_id, hgnc_symbol or ensembl_gene_id.')
        }
    }
    if (!is.null(gene.annotation)) {
        if (nrow(gene.annotation) != nrow(assays[[1]])) {
            stop('The number of rows in gene.annotation and the number of rows in the assays should be the same.')
        } else if (!sum(row.names(gene.annotation) == row.names(assays[[1]])) == nrow(gene.annotation))  {
            stop('The row names in gene.annotation and the row names datastes should be identical.')
        } else if (!gene.group %in% colnames(gene.annotation)) {
            stop('The gene.group should be a column name in the gene.annotation.')
        }
    }
    if(add.housekeeping.genes & is.null(gene.group)){
        stop('To add housekeeping genes, the gene.group should be specified.')
    }
    if(add.housekeeping.genes & is.null(gene.annotation) & !create.gene.annotation){
        stop('To add housekeeping genes, the create.gene.annotation should be TRUE.')
    }
    # grammar
    if(length(assays) == 1){
        assay.n <- 'assay'
    } else assay.n <- 'assays'
    # remove lowly expressed genes ####
    if(remove.lowly.expressed.genes){
        printColoredMessage(message = paste0('-- Remove lowly expressed genes from the ', raw.count.assay.name, ' assay.'),
                            color = 'magenta',
                            verbose = verbose)
        library.size <- colSums(assays[[raw.count.assay.name]])
        cpm.data <- cpm(y = assays[[raw.count.assay.name]], lib.size = NULL)
        cpm.cutoff <- round(count.cutoff/median(library.size) * 1e6, digits = 2)
        if(!is.null(biological.group)){
            sample.size <- min(table(sample.annotation[[biological.group]]))
        }
        if(!is.null(minimum.proportion)){
            sample.size <- round(ncol(cpm.data) * minimum.proportion, digits = 0)
        }
        if(is.null(minimum.proportion) & is.null(biological.group)){
            sample.size <- ncol(cpm.data)
        }
        keep.genes <- rowSums(cpm.data >= cpm.cutoff) >= sample.size
        printColoredMessage(
            message = paste0(
                sum(keep.genes),
                ' of ',
                nrow(cpm.data),
                ' genes with expression cpm cutoff => ',
                cpm.cutoff,
                ' in at least ',
                sample.size,
                ' samples are kept as highly expressed genes.'),
            color = 'blue',
            verbose = verbose
        )
        names.assays <-  names(assays)
        assays <- lapply(
            names(assays),
            function(x) assays[[x]][keep.genes ,])
        names(assays) <- names.assays
        rm(cpm.data, library.size)
        if(!is.null(gene.annotation)){
            gene.annotation <- gene.annotation[keep.genes , ]
        }
    }
    ### calculate library size
    if(calculate.library.size){
        printColoredMessage(
            message = ' -- Calculating library size.',
            color = 'magenta',
            verbose = verbose
        )
        printColoredMessage(
            message = 'Note, to calculate library size, a raw count data without any transformation should be provided.',
            color = 'red',
            verbose = verbose
        )
        library.size <- colSums(assays[[raw.count.assay.name]])
        printColoredMessage(
            message = 'The library size is calculated, with summaries (in millions):',
            color = 'blue',
            verbose = verbose
        )
        if(verbose) print(summary(library.size/10^6), color = 'blue')
    }
    # sample annotation ####
    printColoredMessage(
        message = '-- Sample annotation:',
        color = 'magenta',
        verbose = verbose
    )
    if(!is.null(sample.annotation)){
        printColoredMessage(
            message = 'The sample.annotation will be added to the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        if(calculate.library.size){
            sample.annotation[['library.size']] <- library.size
        }
    }
    if(is.null(sample.annotation) & create.sample.annotation){
        printColoredMessage(
            message = 'A sample.annotation will be created and added to the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        if(calculate.library.size){
            sample.annotation <- data.frame(
                sample.ids = colnames(assays[[1]]),
                library.size = library.size
            )
        } else{ sample.annotation <- data.frame(
                sample.ids = colnames(assays[[1]]))
        }
    }
    if (is.null(sample.annotation) & !create.sample.annotation){
        printColoredMessage(
            message = 'The SummarizedExperiment object will not contain sample annotation.',
            color = 'blue',
            verbose = verbose)
    }
    # gene annotation ####
    printColoredMessage(
        message = '-- Gene annotation:',
        color = 'magenta',
        verbose = verbose
    )
    if(!is.null(gene.annotation)){
        printColoredMessage(
            message = 'The gene.annotation will be added to the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
    } else if(is.null(gene.annotation) & create.gene.annotation){
        printColoredMessage(
            message = 'A gene.annotation (rowData) will be created and added to the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        gene.annotation <- data.frame(gene.ids = row.names(assays[[1]]))
        colnames(gene.annotation) <- gene.group
    } else if (is.null(gene.annotation)){
        printColoredMessage(
            message = 'The SummarizedExperiment object will not contain gene annotation.',
            color = 'blue',
            verbose = verbose)
    }
    # add gene details ####
    printColoredMessage(
        message = '-- Gene details:',
        color = 'magenta',
        verbose = verbose
    )
    if(add.gene.details){
        if(is.null(gene.details)){
            printColoredMessage(
                message = 'The gene.details is not specified, some pre-set details will be added to the gene annotation.',
                color = 'blue',
                verbose = verbose)
            printColoredMessage(
                message = 'Obtain the pre-set gene details from the bioMart R package, this may take a few minutes.',
                color = 'blue',
                verbose = verbose)
            ensembl <- useMart('ensembl') # 24 March 2022
            ensembl <- useDataset(
                mart = ensembl,
                'hsapiens_gene_ensembl'
            )
            bioMart.geneAnnot <- biomaRt::getBM(
                attributes = c(
                    'entrezgene_id',
                    'hgnc_symbol',
                    'gene_biotype',
                    'ensembl_gene_id',
                    'description',
                    'chromosome_name'),
                mart = ensembl
            )
            bioMart.geneAnnot <- bioMart.geneAnnot[!duplicated(bioMart.geneAnnot[[gene.group]]), ]
            gene.annotation <- left_join(
                x = gene.annotation,
                y = bioMart.geneAnnot,
                by = gene.group,
                multiple = 'first')
        } else if(!is.null(gene.details)){
            printColoredMessage(
                message = 'Obtain gene details from the bioMart R package.',
                color = 'blue',
                verbose = verbose)
            ensembl <- useMart('ensembl') # 24 March 2022
            ensembl <- useDataset(
                mart = ensembl,
                'hsapiens_gene_ensembl'
            )
            attributes.list <- biomaRt::listAttributes(mart = ensembl)
            if(sum(gene.details %in% attributes.list$name) == 0){
                stop('Non of the provided gene.details are found in the attributes list (biomaRt::listAttributes) in the biomaRt.')
            } else {
                printColoredMessage(
                    message = paste0(
                        sum(gene.details %in% attributes.list$name),
                        ' of ',
                        length(gene.details),
                        'gene.details are found.'),
                        color = 'blue',
                        verbose = verbose)
            }
            gene.details <- unique(gene.details, gene.group)
            bioMart.geneAnnot <- biomaRt::getBM(
                attributes = gene.details,
                mart = ensembl
            )
            bioMart.geneAnnot <- bioMart.geneAnnot[!duplicated(bioMart.geneAnnot[[gene.group]]), ]
            gene.annotation <- left_join(
                x = gene.annotation,
                y = bioMart.geneAnnot,
                multiple = 'first',
                by = gene.group)
        }
    } else{
        printColoredMessage(
            message = 'Any extra gene details are not specified.',
            color = 'blue',
            verbose = verbose
        )
    }
    # add housekeeping genes list ####
    if(add.housekeeping.genes){
        printColoredMessage(
            message = '-- Add several lists of housekeeping genes to the gene annotation:',
            color = 'magenta',
            verbose = verbose
        )
        kh.im.genes <- hk_immunStroma
        keep.cols <- c(
            which(colnames(kh.im.genes) %in% gene.group),
            4:9)
        gene.annotation <- left_join(
            x = gene.annotation,
            y = kh.im.genes[ , keep.cols],
            by = gene.group,
            multiple = 'first'
            )
        printColoredMessage(
            message = 'Seven different lists of housekeeping genes are added to the gene annotation.',
            color = 'blue',
            verbose = verbose
        )
        nb.hk.genes <- lapply(colnames(kh.im.genes)[4:9], function(x) sum(gene.annotation[[x]] == 'yes'))
        names(nb.hk.genes) <- colnames(kh.im.genes)[4:9]
        if(verbose) print(kable(unlist(nb.hk.genes),
                                caption = 'Number of genes in each list of housekeeping genes:',
                                col.names = 'nb.genes'))
    }
    # add immune and stromal genes signatures ####
    if(add.immunStroma.genes){
        printColoredMessage(
            message = '-- Add immune and stromal genes signature to the gene annotation:',
            color = 'magenta',
            verbose = verbose
        )
        kh.im.genes <- hk_immunStroma
        keep.cols <- c(
            which(colnames(kh.im.genes) %in% gene.group),
            10:ncol(kh.im.genes))
        gene.annotation <- as.data.frame(left_join(
            x = gene.annotation,
            y = kh.im.genes[ , keep.cols],
            by = gene.group,
            multiple = 'first'
            ))
        printColoredMessage(
            message = 'The immune and stromal genes signature from Kosuke Yoshihara et.al are added.',
            color = 'blue',
            verbose = verbose
        )
        nb.genes <- lapply(colnames(kh.im.genes)[10:11], function(x) sum(gene.annotation[[x]] == 'yes'))
        names(nb.genes) <- colnames(kh.im.genes)[10:11]
        if(verbose) print(
            kable(unlist(nb.genes),
                  caption = 'Number of genes in the immune and stromal gene signatures:',
                  col.names = 'nb.genes'))
    }
    # outputs ####
    printColoredMessage(
        message = '-- Create a SummarizedExperiment object:',
        color = 'magenta',
        verbose = verbose
    )
    if (is.null(gene.annotation) & is.null(sample.annotation) & is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(assays = assays)
    } else if (!is.null(gene.annotation) & is.null(sample.annotation) & is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(
            assays = assays,
            rowData = gene.annotation)
    } else if (is.null(gene.annotation) & !is.null(sample.annotation) & is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(
            assays = assays,
            colData = DataFrame(sample.annotation),
            metadata = metaData
        )
    } else if (is.null(gene.annotation) & is.null(sample.annotation) & !is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(
            assays = assays,
            metadata = metaData)
    } else if (!is.null(gene.annotation) & !is.null(sample.annotation) & is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(
            assays = assays,
            rowData = gene.annotation,
            colData = sample.annotation)
    } else if (is.null(gene.annotation) & !is.null(sample.annotation) & is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(
            assays = assays,
            rowData = gene.annotation,
            colData = DataFrame(sample.annotation)
        )
    } else if (!is.null(gene.annotation) & !is.null(sample.annotation) & !is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(
            assays = assays,
            rowData = gene.annotation,
            colData = sample.annotation,
            metadata = metaData
        )
    }
    printColoredMessage(
        message = paste0('A summarizedExperiment object is created with:'),
        color = 'blue',
        verbose = verbose
    )
    printColoredMessage(
        message = paste0(
            '-',
            nrow(se.obj),
            ' measurements (e.g. genes) and ',
            ncol(se.obj),
            ' assays (e.g.samples)'
        ),
        color = 'blue',
        verbose = verbose
    )
    printColoredMessage(
        message = paste0('-', length(assays(se.obj)), ' data sets (assays)'),
        color = 'blue',
        verbose = verbose
    )
    if (!is.null(sample.annotation))
        printColoredMessage(
            message = paste0('-', ncol(colData(se.obj)), ' annotations for the samples'),
            color = 'blue',
            verbose = verbose
        )
    if (!is.null(gene.annotation))
        printColoredMessage(
            message = paste0('-', ncol(rowData(se.obj)), ' annotations for the genes'),
            color = 'blue',
            verbose = verbose
        )
    if (!is.null(metaData))
        printColoredMessage(message = '- a metadata',
                            color = 'blue',
                            verbose = verbose)
    printColoredMessage(message = '------------The createSeObj function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}
