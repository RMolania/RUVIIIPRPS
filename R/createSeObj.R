#' is used to create a summarizedExperiment object.
#'
#' @param assays A list of assays.
#' @param remove.lowly.expressed.genes 	logical, if TRUE then remove lowly expressed genes.
#' @param raw.count.assay.name The name of raw counts data in the assays.
#' @param count.cutoff numeric. Minimum count required for at least some samples.
#' @param biological.group A column name in sample annotation that specify biological groups. The smallest population is considred for relvling lowly expressed genes.
#' @param minimum.proportion numeric. In large sample situations, the minimum proportion of samples in a group that a gene needs to be expressed in. See Details below for the exact formula.
#' @param calculate.library.size logical, if TRUE then library size is calculated on the raw.count.assay.name.
#' @param sample.annotation a data frame, contains information for samples.
#' @param create.sample.annotation logical, if TRUE then a sample annotation the initially contains column names of the assays.
#' @param gene.annotation a data frame, contains information for genes.
#' @param create.gene.annotation logical, if TRUE then a gene annotation that initially contains row names of the assays.
#' @param add.gene.details logical, if TRUE then a pr-set or provided gene details will be added to the gene annotation.
#' @param gene.group A name of gene class in the assays, this should be in 'entrezgene_id', 'hgnc_symbol', 'ensembl_gene_id'.
#' @param gene.details A list of gene details to be added to the gene annotation.
#' @param add.housekeeping.genes logical. if TRUE then several sets of publicly available sets of "housekeeping" will be added to the gene annotation.
#' @param add.immunStroma.genes Logical. If TRUE, the immune and stromal genes signature from Kosuke Yoshihara et.al will be added.
#' @param metaData Any metadata data for in the data.
#' @param verbose logical, if TRUE shows the messages.

#' @return a summarizedExperiment that contains assays, gene annotation, samples annotation and metadata.

#' @importFrom Matrix colSums
#' @importFrom dplyr left_join
#' @importFrom SummarizedExperiment assay SummarizedExperiment rowData
#' @importFrom biomaRt getBM useMart useDataset
#' @importFrom S4Vectors DataFrame
#' @importFrom utils read.table
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
    if(class(assays)!= 'list'){
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
        if(count.cutoff <0){
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
        stop('A gene annotattion has been provided, then create.gene.annotation must be FALSE.')
    }
    if(add.gene.details & is.null(gene.annotation) & !create.gene.annotation){
        stop('To add add.gene.details , the create.gene.annotation should be TRUE.')
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
    } else{
        assay.n <- 'assays'
    }
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
        hk.gene.lists <- read.table(
            file = '../Data/houskeeping.genes.lists.txt'
            )
        keep.cols <- c(
            which(colnames(hk.gene.lists) %in% gene.group),
            4:9)
        gene.annotation <- left_join(
            x = gene.annotation,
            y = hk.gene.lists[ , keep.cols],
            by = gene.group,
            multiple = 'first'
            )
        printColoredMessage(
            message = 'Seven different lists of housekeeping genes are added to the gene annotation.',
            color = 'blue',
            verbose = verbose
        )
        nb.hk.genes <- lapply(colnames(hk.gene.lists)[4:9], function(x) sum(gene.annotation[[x]]))
        names(nb.hk.genes) <- colnames(hk.gene.lists)[4:9]
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
        hk.gene.lists <- read.table(
            file = '../Data/houskeeping.genes.lists.txt'
            )
        keep.cols <- c(
            which(colnames(hk.gene.lists) %in% gene.group),
            10:ncol(hk.gene.lists))
        gene.annotation <- left_join(
            x = gene.annotation,
            y = hk.gene.lists[ , keep.cols],
            by = gene.group,
            multiple = 'first'
            )
        printColoredMessage(
            message = 'The immune and stromal genes signature from Kosuke Yoshihara et.al are added.',
            color = 'blue',
            verbose = verbose
        )
        nb.hk.genes <- lapply(colnames(hk.gene.lists)[10:11], function(x) sum(gene.annotation[[x]]))
        names(nb.hk.genes) <- colnames(hk.gene.lists)[10:11]
        if(verbose) print(
            kable(unlist(nb.hk.genes),
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
        message = paste0('-', length(assays(se.obj)), ' data sets (assay)'),
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


# red.cols <- c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol', 'gene_biotype')
# red.cols <- unlist(lapply(red.cols, function(x){
#     index <- grep(x, colnames(gene.annotation))
#     if(length(index) > 1) index[2:length(index)]
# }))
# if(is.numeric(red.cols)) gene.annotation <- gene.annotation[ , -c(red.cols)]
