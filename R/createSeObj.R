#' Create a SummarizedExperiment object.

#' @author Ramyar Molania

#' @description
#' This function creates a SummarizedExperiment object from tabular expression data set(s) and sample annotation (if
#' available). The function can identify and remove lowly expressed genes if a raw count data is provided, add a range
#' of annotations for individual genes e.g. biotype, chromosome names, ... . , estimate tumor purity, and provide several
#' sets of housekeeping genes and a immune and stromal gene signature. The housekeeping gene sets could be suitable
#' negative control genes for the RUV methods. The immune and stromal gene signature could be used to estimate tumor purity
#' variation in cancer RNA-seq data.

#' @details
#' The SummarizedExperiment object is a data structure used in the R for representing and manipulating high-dimensional
#' experimental data. Here are some key features and components of the SummarizedExperiment object:
#' Assays:
#' The SummarizedExperiment allows for the incorporation of multiple assays (data). Each assay is a separate matrix of data
#' associated with the same features and samples. Rows typically represent features (e.g., genes, transcripts), and columns
#' represent samples or experimental conditions.
#' Row and Column Metadata:
#' The SummarizedExperiment object includes metadata associated with both rows and columns. Row metadata can contain
#' information about the features, such as gene annotations or genomic coordinates. Column metadata may include sample
#' information, experimental conditions, or other relevant details.
#' Metadata
#' The Metadata of the SummarizedExperiment allows for flexibility in terms of data types and structures. This makes it
#' suitable for saving metrics and plots in the RUVIIIPRPS R package. We refer to the SummarizedExperiment R package for
#' more details.


#' @param data.sets List. A list containing assay(s) or expression data. For individual datasets, genes should be arranged
#' in rows and samples in columns. If multiple datasets are provided, ensure that the row names of the assays are in the
#' same order.
#' @param raw.count.assay.name Symbol. The name of raw counts RNA-seq data within the list of assay(s) or expression data.
#' Raw counts data must be included if 'remove.lowly.expressed.genes' or "calculate.library.size" is set to 'TRUE'.
#' @param remove.lowly.expressed.genes 	Logical. If 'TRUE', the function identifies and removes lowly expressed genes from
#' the assay specified in the 'raw.count.assay.name' argument. The default is 'FALSE'.
#' @param count.cutoff Numeric. Minimum count required for at least some sample groups. If the 'biological.group' argument
#' is set to 'NULL', all samples will be considered as single group. Otherwise, the smallest subgroups of the
#' 'biological.group' will be taken into account to remove lowly expressed genes. We refer to the filterByExpr function
#' from the edgeR R package.
#' @param biological.group Symbol. Indicates a column name in the sample annotation, specifying the biological groups.
#' The smallest biological groups will be taken into account for removing lowly expressed genes. If is 'NULL', all samples
#' will be considered as single group.
#' @param minimum.proportion Numeric. In large sample situations, the minimum proportion of samples within a group in which
#' a gene must be expressed. The default is 0.5.
#' @param calculate.library.size Logical. If 'TRUE', then library size (total counts) is calculated using the data provided
#' in 'raw.count.assay.name'. The library size will be calculated after removing lowly expressed genes.
#' @param estimate.tumor.purity Symbol. A symbol indicating which methods to be used to estimt the tumour purity. Opoptiona
#' include:'estimate', 'singscore', 'both' and 'NULL'. If 'estimate' is selected, the function will applied the ESTIMATE method to
#' estimate tumor purity. If 'singscore', the function will utilize the 'singscore' method and if 'both', the two methods
#' will be applied. The default is 'NULL'.
#' @param assay.name.to.estimate.purity Symbol.The name of an assay data within the list of assay(s) or expression data to
#' be used to estimate tumor purity.
#' @param sample.annotation A data frame containing information for individual samples in the assay(s). The order of row
#' names of the sample annotation should match with the columns names of the assay(s).
#' @param create.sample.annotation Logical. If 'TRUE', a sample annotation will be generated, initially contains the
#' column names of the assay(s). Additional variables such as library size will be included in the sample annotation.
#' @param gene.annotation A data frame, containing details e.g. chromosome names,GC content,... for individual genes. The
#' list of housekeeping genes, and immune and stromal genes will be included in the gene annotation.
#' @param create.gene.annotation Logical. If 'TRUE', a gene annotation will be generated, initially contains row names of
#' the assay(s). Additional gene lists and details will be included in the gene annotation.
#' @param add.gene.details Logical. If 'TRUE', preset or provided gene details in the 'gene.details' argument will be
#' added to the gene annotation.
#' @param gene.group Symbol. The name of a gene class of the row names of the assay(s). This must be one of the 'entrezgene_id',
#' 'hgnc_symbol', 'ensembl_gene_id'.
#' @param gene.details Symbol. A symbol or a vector symbols indicating the gene details to be included to the gene annotation.
#' @param add.housekeeping.genes Logical. if 'TRUE', several sets of publicly available "housekeeping" will be included in
#' the gene annotation. The housekeeping could be potentially used as negative control genes for the RUV normalization.
#' @param add.immun.stroma.genes Logical. If 'TRUE', the immune and stromal genes signature from Kosuke Yoshihara et.al will
#' be added to the gene annotation. These genes signatures, can be used to estimate tumor purity in cancer RNA-seq data.
#' @param metaData Any metadata data. The metadata can be in any format and dimensions.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return A SummarizedExperiment object containing assay(s) and also samples annotation, gene annotation and metadata
#' if there are specified.

#' @importFrom SummarizedExperiment SummarizedExperiment assay rowData
#' @importFrom tidyestimate filter_common_genes estimate_score
#' @importFrom biomaRt getBM useMart useDataset
#' @importFrom singscore rankGenes simpleScore
#' @importFrom Matrix colSums rowSums
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr left_join
#' @importFrom knitr kable
#' @importFrom edgeR cpm
#' @export

createSeObj <- function(
        data.sets ,
        raw.count.assay.name = NULL,
        remove.lowly.expressed.genes = FALSE,
        count.cutoff = 10,
        biological.group = NULL,
        minimum.proportion = 0.5,
        calculate.library.size = FALSE,
        estimate.tumor.purity = NULL,
        assay.name.to.estimate.purity = NULL,
        sample.annotation = NULL,
        create.sample.annotation = FALSE,
        gene.annotation = NULL,
        create.gene.annotation = FALSE,
        add.gene.details = FALSE,
        gene.group = NULL,
        gene.details = NULL,
        add.housekeeping.genes = FALSE,
        add.immun.stroma.genes = FALSE,
        metaData = NULL,
        verbose = TRUE)
{
    printColoredMessage(message = '------------The createSeObj function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if(!inherits(data.sets, what = 'list')){
        stop('The "data.sets" must be a list of expression dataste(s) (genes in rows and samples in columns.).')
    }
    if(isTRUE(calculate.library.size)){
        if(is.null(sample.annotation) & isFALSE(create.sample.annotation)){
            stop('To add the calculated library size, either a "sample.annotation" should be provided or set the create.sample.annotation=TRUE.')
        }
    }
    if(!is.null(estimate.tumor.purity)){
        if(is.null(sample.annotation) & isFALSE(create.sample.annotation)){
            stop('To add the calculated tumour purity, either a "sample.annotation" should be provided or set the create.sample.annotation=TRUE.')
        } else if(is.null(assay.name.to.estimate.purity)){
            stop('To estimate the tumour purity, the "assay.name.to.estimate.purity" must be provided.')
        }
    }
    if(isTRUE(add.gene.details)){
        if(is.null(gene.annotation) & isFALSE(create.gene.annotation)){
            stop('To add gene details, either a "gene.annotation" should be provided or set the "create.gene.annotation=TRUE".')
        }
    }
    if(isTRUE(add.housekeeping.genes)){
        if(is.null(gene.annotation) & isFALSE(create.gene.annotation)){
            stop('To add housekeeping genes, either a "gene.annotation" should be provided or set the "create.gene.annotation=TRUE".')
        }
    }
    if(isTRUE(add.immun.stroma.genes)){
        if(is.null(gene.annotation) & isFALSE(create.gene.annotation)){
            stop('To add immun stroma genes, either a "gene.annotation" should be provided or set the "create.gene.annotation=TRUE".')
        }
    }

    ## dimensions and orders of the assays ####
    if(length(data.sets) > 1){
        dim.data <- unlist(lapply(
            names(data.sets),
            function(x) c(nrow(data.sets[[x]]), ncol(data.sets[[x]])))
        )
        if(length(unique(dim.data)) !=2 ){
            stop('The "data.sets" must have matching dimensions.')
        }
        all.assays <- combn(1:length(data.sets), 2)
        m.out <- lapply(
            1:ncol(all.assays),
            function(x){
                if(!all.equal(row.names(data.sets[[all.assays[ 1, x]]]), row.names(data.sets[[all.assays[ 2, x]]])))
                    stop('The row names of the "data.sets" must be in the same order.')
            })
        m.out <- lapply(
            1:ncol(all.assays),
            function(x){
                if(!all.equal(colnames(data.sets[[all.assays[ 1, x]]]),colnames(data.sets[[all.assays[ 2, x]]])))
                    stop('The column names of the "data.sets" must be in the same order.')
            })
    }
    ## check raw.count.assay.name ####
    if(!is.null(raw.count.assay.name)){
        if(length(raw.count.assay.name) > 1){
            stop('The "raw.count.assay.name" must contain a single assay name.')
        }
        if(!raw.count.assay.name %in% names(data.sets)){
            stop('The "raw.count.assay.name" must be in the assay list.')
        }
    }
    ## remove lowly expressed genes ####
    if(isTRUE(remove.lowly.expressed.genes)){
        if(is.null(raw.count.assay.name)){
            stop('To find and remove lowly expressed genes, the "raw.count.assay.name" must be provided.')
        }
        if(!is.null(biological.group)){
            if(is.null(sample.annotation)){
                stop('To find "biological.group", please provide a sample annotation that contains the "biological.group" variable.')
            }
            if(!biological.group %in% colnames(sample.annotation)){
                stop('The "biological.group" cannot be found in the sample annotation.')
            }
        }
        if(minimum.proportion > 1 | minimum.proportion < 0){
            stop('The "minimum.proportion" must between 0 to 1.')
        }
        if(is.null(count.cutoff)){
            stop('The count.cutoff cannot be empty.')
        }
        if(count.cutoff < 0){
            stop('The value of the "count.cutoff" cannot be negative.')
        }
    }
    ## library size calculation ####
    if(is.null(raw.count.assay.name) & isTRUE(calculate.library.size)){
        stop('To calculate library size, raw count data must be provided.')
    }
    ## sample annotation ####
    if(!is.null(sample.annotation)){
        if(nrow(sample.annotation) != ncol(data.sets[[1]])){
            stop('The number of rows in "sample.annotation" and the number of columns in the "data.sets" must be the same.')
        } else if (!all.equal(row.names(sample.annotation), colnames(data.sets[[1]])) ){
            stop('The order and lables of the row names of "sample.annotation" and colummn names of the "data.sets" should be the same.')
        }
        if(create.sample.annotation){
            stop('A "sample.annotation" is provided, then "create.sample.annotation" must be set to "FALSE".')
        }
    }
    ## gene annotation ####
    if(!is.null(gene.annotation) & isTRUE(create.gene.annotation)){
        stop('A gene annotation is provided, then "create.gene.annotation" must be "FALSE".')
    }
    if(isTRUE(add.gene.details) & is.null(gene.annotation) & !create.gene.annotation){
        stop('To add "add.gene.details" , the "create.gene.annotation" must be set to "TRUE".')
    }
    if(isTRUE(add.housekeeping.genes) & is.null(gene.annotation) & !create.gene.annotation){
        stop('To add "add.housekeeping.genes", gene annotation must be provided or "create.gene.annotation" must be set to "TRUE" .')
    }
    if(isTRUE(add.immun.stroma.genes) & is.null(gene.annotation) & !create.gene.annotation){
        stop('To add "add.immun.stroma.genes", gene annotation must be provided or "create.gene.annotation" must set to "TRUE".')
    }
    if(isTRUE(add.gene.details) & is.null(gene.group)){
        stop('To add gene details, the "gene.group" must be specified (entrezgene_id, hgnc_symbol and ensembl_gene_id).')
    }
    if(isTRUE(add.gene.details) &!is.null(gene.group)){
        if(!gene.group %in% c('entrezgene_id', 'hgnc_symbol', 'ensembl_gene_id')){
            stop('The "gene.group" must be one of the "entrezgene_id", "hgnc_symbol" or "ensembl_gene_id".')
        }
    }
    if (!is.null(gene.annotation)) {
        if (nrow(gene.annotation) != nrow(data.sets[[1]])) {
            stop('The number of rows in "gene.annotation" and the number of rows in the "data.sets" must be the same.')
        } else if (!sum(row.names(gene.annotation) == row.names(data.sets[[1]])) == nrow(gene.annotation))  {
            stop('The row names in "gene.annotation" and the row names datastes should be identical.')
        } else if (!gene.group %in% colnames(gene.annotation)) {
            stop('The "gene.group" must be a column name in the "gene.annotation".')
        }
    }
    if(isTRUE(add.housekeeping.genes) & is.null(gene.group)){
        stop('To add housekeeping genes, the "gene.group" must be specified.')
    }
    if(isTRUE(add.housekeeping.genes) & is.null(gene.annotation) & !create.gene.annotation){
        stop('To add housekeeping genes, the "create.gene.annotation" must be "TRUE".')
    }
    # grammar
    if(length(data.sets) == 1){
        assay.n <- 'assay'
    } else assay.n <- 'assays'
    # remove lowly expressed genes ####
    if(remove.lowly.expressed.genes){
        printColoredMessage(message = paste0('-- Remove lowly expressed genes from the ', raw.count.assay.name, ' assay.'),
                            color = 'magenta',
                            verbose = verbose)
        library.size <- Matrix::colSums(data.sets[[raw.count.assay.name]])
        cpm.data <- edgeR::cpm(y = data.sets[[raw.count.assay.name]], lib.size = NULL)
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
        keep.genes <- Matrix::rowSums(cpm.data >= cpm.cutoff) >= sample.size
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
        names.assays <-  names(data.sets)
        data.sets <- lapply(
            names.assays,
            function(x) data.sets[[x]][keep.genes ,])
        names(data.sets) <- names.assays
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
        library.size <- Matrix::colSums(data.sets[[raw.count.assay.name]])
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
            message = 'The "sample.annotation" will be added to the SummarizedExperiment object.',
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
                sample.ids = colnames(data.sets[[1]]),
                library.size = library.size
            )
        } else sample.annotation <- data.frame(sample.ids = colnames(data.sets[[1]]))
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
            verbose = verbose)
    } else if(is.null(gene.annotation) & isTRUE(create.gene.annotation)){
        printColoredMessage(
            message = 'A gene.annotation (rowData) will be created and added to the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        gene.annotation <- data.frame(gene.ids = row.names(data.sets[[1]]))
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
    if(isTRUE(add.gene.details)){
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
            gene.annotation <- dplyr::left_join(
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
            gene.annotation <- dplyr::left_join(
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
        hk.im.genes <- hk_immunStroma
        keep.cols <- c(which(colnames(hk.im.genes) %in% gene.group), 4:9)
        gene.annotation <- dplyr::left_join(
            x = gene.annotation,
            y = hk.im.genes[ , keep.cols],
            by = gene.group,
            multiple = 'first'
            )
        printColoredMessage(
            message = 'Seven different lists of housekeeping genes are added to the gene annotation.',
            color = 'blue',
            verbose = verbose
        )
        nb.hk.genes <- lapply(
            colnames(hk.im.genes)[4:9],
            function(x) sum(gene.annotation[[x]]))
        names(nb.hk.genes) <- colnames(hk.im.genes)[4:9]
        if(verbose) print(kable(unlist(nb.hk.genes),
                                caption = 'Number of genes in each list of housekeeping genes:',
                                col.names = 'nb.genes'))
    }
    # add immune and stroma genes signatures ####
    if(add.immun.stroma.genes){
        printColoredMessage(
            message = '-- Add immune and stromal genes signature to the gene annotation:',
            color = 'magenta',
            verbose = verbose
        )
        hk.im.genes <- hk_immunStroma
        keep.cols <- c(
            which(colnames(hk.im.genes) %in% gene.group),
            10:ncol(hk.im.genes))
        gene.annotation <- as.data.frame(dplyr::left_join(
            x = gene.annotation,
            y = hk.im.genes[ , keep.cols],
            by = gene.group,
            multiple = 'first'
            ))
        printColoredMessage(
            message = 'The immune and stromal genes signature from Kosuke Yoshihara et.al are added.',
            color = 'blue',
            verbose = verbose
        )
        nb.genes <- lapply(
            colnames(hk.im.genes)[10:11],
            function(x) sum(gene.annotation[x]))
        names(nb.genes) <- colnames(hk.im.genes)[10:11]
        if(verbose) print(
            kable(unlist(nb.genes),
                  caption = 'Number of genes in the immune and stromal gene signatures:',
                  col.names = 'nb.genes'))
    }
    # estimate tumor purity ####
    if(!is.null(estimate.tumor.purity)){
        printColoredMessage(message = '-- Estimate tumour purity:',
            color = 'magenta',
            verbose = verbose)
        if(estimate.tumor.purity == 'estimate'){
            printColoredMessage(
                message = '-- Estimate tumour purity using the ESTIMATE method:',
                color = 'blue',
                verbose = verbose)
            tumour.purity <- tidyestimate::filter_common_genes(
                df = data.sets[[assay.name.to.estimate.purity]],
                id = "hgnc_symbol",
                tidy = FALSE,
                tell_missing = verbose,
                find_alias = TRUE)
            tumour.purity <- tidyestimate::estimate_score(
                df = tumour.purity,
                is_affymetrix = TRUE)
            tumour.purity <- tumour.purity$purity
            sample.annotation[['tumour.purity']] <- tumour.purity
        } else if (estimate.tumor.purity == 'singscore'){
            printColoredMessage(
                message = '-- Estimate tumour purity using the singscore method:',
                color = 'blue',
                verbose = verbose)
            im.str.gene.sig <- hk_immunStroma$immune.gene.signature == 'TRUE' |
                hk_immunStroma$stromal.gene.signature == 'TRUE'
            if(gene.group == "entrezgene_id"){
                im.str.gene.sig <- hk_immunStroma$entrezgene_id[im.str.gene.sig]
            } else if(gene.group == 'hgnc_symbol'){
                im.str.gene.sig <- hk_immunStroma$hgnc_symbol[im.str.gene.sig]
            } else if(gene.group == 'ensembl_gene_id')
                im.str.gene.sig <- hk_immunStroma$ensembl_gene_id[im.str.gene.sig]
            tumour.purity <- singscore::rankGenes(data.sets[[assay.name.to.estimate.purity]])
            tumour.purity <- singscore::simpleScore(
                rankData = tumour.purity,
                upSet = im.str.gene.sig)
            tumour.purity <- tumour.purity$TotalScore
            sample.annotation[['tumour.purity']] <- tumour.purity
        } else if (estimate.tumor.purity == 'both'){
            printColoredMessage(
                message = '-- Estimate tumour purity using both ESTIMATE and singscore methods:',
                color = 'blue',
                verbose = verbose)
            tumour.purity <- tidyestimate::filter_common_genes(
                df = data.sets[[assay.name.to.estimate.purity]],
                id = "hgnc_symbol",
                tidy = FALSE,
                tell_missing = verbose,
                find_alias = TRUE)
            tumour.purity <- tidyestimate::estimate_score(
                df = tumour.purity,
                is_affymetrix = TRUE)
            tumour.purity.estimate <- tumour.purity$purity
            im.str.gene.sig <- hk_immunStroma$immune.gene.signature == 'TRUE' |
                hk_immunStroma$stromal.gene.signature == 'TRUE'
            if(gene.group == "entrezgene_id"){
                im.str.gene.sig <- hk_immunStroma$entrezgene_id[im.str.gene.sig]
            } else if(gene.group == 'hgnc_symbol'){
                im.str.gene.sig <- hk_immunStroma$hgnc_symbol[im.str.gene.sig]
            } else if(gene.group == 'ensembl_gene_id')
                im.str.gene.sig <- hk_immunStroma$ensembl_gene_id[im.str.gene.sig]
            tumour.purity <- singscore::rankGenes(data.sets[[assay.name.to.estimate.purity]])
            tumour.purity <- singscore::simpleScore(
                rankData = tumour.purity,
                upSet = im.str.gene.sig)
            tumour.purity.singscore <- tumour.purity$TotalScore
            sample.annotation[['tumour.purity.estimate']] <- tumour.purity.estimate
            sample.annotation[['tumour.purity.singscore']] <- tumour.purity.singscore
        }
    }
    # outputs ####
    printColoredMessage(
        message = '-- Create a SummarizedExperiment object:',
        color = 'magenta',
        verbose = verbose
    )
    if (is.null(gene.annotation) & is.null(sample.annotation) & is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(assays = data.sets)
    } else if (!is.null(gene.annotation) & is.null(sample.annotation) & is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(
            assays = data.sets,
            rowData = gene.annotation)
    } else if (is.null(gene.annotation) & !is.null(sample.annotation) & is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(
            assays = data.sets,
            colData = DataFrame(sample.annotation),
            metadata = metaData
        )
    } else if (is.null(gene.annotation) & is.null(sample.annotation) & !is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(
            assays = data.sets,
            metadata = metaData)
    } else if (!is.null(gene.annotation) & !is.null(sample.annotation) & is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(
            assays = data.sets,
            rowData = gene.annotation,
            colData = sample.annotation)
    } else if (is.null(gene.annotation) & !is.null(sample.annotation) & is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(
            assays = data.sets,
            rowData = gene.annotation,
            colData = DataFrame(sample.annotation)
        )
    } else if (!is.null(gene.annotation) & !is.null(sample.annotation) & !is.null(metaData)) {
        se.obj <- SummarizedExperiment::SummarizedExperiment(
            assays = data.sets,
            rowData = gene.annotation,
            colData = sample.annotation,
            metadata = metaData
        )
    }
    printColoredMessage(
        message = paste0('A SummarizedExperiment object is created with:'),
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


