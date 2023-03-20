#' is used to assess the performance of the normalisation of a SummarizedExperiment class object.
#' Several assessment will be performed:
#' 1) PCA plot of each categorical variable.
#' 2) Silhouette and ARI computed on categorical variable.
#' 3) Combined Silhouette plot of the combined pair of categorical variable.
#' 4) Linear regression between the first cumulative PC and continuous variable.
#' 5) Spearman correlation between gene expression and continuous variable.
#'
#' @param se Dataset that will be used to assess the performance of the normalisation of the data.
#' @param assay_names Optional string or list of strings for selection of the names
#' of the assays of the SummarizedExperiment class object.
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#' @param cat_var_label String or vector of strings of the label of categorical variable(s) such as
#' sample types or batches from colData(se).
#' @param cont_var_label String or vector of strings of the label of continuous variable(s)
#' such as library size from colData(se).
#' @param output_file Path and name of the output file to save the assessments plots in a pdf format.
#' @param n.cores is the number of cpus used for mclapply parallelization. Default is set to 5.
#'
#'
#' @return plots List of assessments plots
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom gridExtra grid.arrange
#' @importFrom SummarizedExperiment assays colData
#' @export

norm_assessment = function(
        se,
        assay_names=NULL,
        apply.log = FALSE,
        cat_var_label,
        cont_var_label,
        #biological_subtypes,
        #library_size,
        #batch,
        output_file=NULL,
        n.cores=5
){
    ### Check se and assay names
    if (!class(se)[1] == 'SummarizedExperiment') {
        stop('Please provide a summarized experiment object.')
    } else if (class(se)[1] == 'SummarizedExperiment' & ncol(colData(se)) == 0){
        stop('The Summarized experiment object does not contain
           sample annotations, please provide colData of the summarized
           experiment object before running the norm_assessment() function.')
    } else if((!is.null(assay_names))&!(assay_names %in% names(assays(se)))){
        stop('The selected assay(s) is/are not in the assays names of the SummarizedExperiment class object.')
    }

    ### Check cat_var_label and cont_var_label
    sample.annot <- as.data.frame(colData(x = se))
    all.var.label<- c(cat_var_label, cont_var_label)
    exist.var.label <- colnames(sample.annot)[colnames(sample.annot) %in% all.var.label]
    if (!sum(all.var.label %in% exist.var.label) == length(all.var.label)){
        stop(
            'Provided variable label from cat_var_label and cont_var_label "',
            paste0(
                all.var.label[!all.var.label %in% exist.var.label], collapse = ' & '),
            '" are not in the colData of the summarized experiment object.')
    }

    ### Compute PCA
    data_pca=RUVPRPS::compute_pca(se,apply.log = apply.log)

    ## Get all the available assays (i.e. normalizations methods)
    if (!is.null(assay_names)){
        normalization=assay_names
    }else{
        normalization=names(
            SummarizedExperiment::assays(se))
    }
    ################# Categorical variable ################
    #nb_catvar=
    biological_subtypes=sample.annot[, cat_var_label]
    ## PCA Color Biology
    colfunc <- colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = 'Spectral')[-6])
    color.subtype<- colfunc(length(unique(biological_subtypes)))
    names(color.subtype) <- levels(biological_subtypes)
    message("PCA based on Biology")
    ### Compute PCA Biology
    PCA_BIO=RUVPRPS::plot_pca(data_pca,
                    variable= biological_subtypes,
                    variable.name =  'Biology',
                     color = color.subtype)

    ## Compute Silhouette based on biology
    message("Silhouette coefficient based on biology")
    silh_bio=RUVPRPS::compute_silhouette(data_pca,
                          cat_var=biological_subtypes)

    ## Compute ARI based on biology
    message("ARI based on biology")
    ari_bio=RUVPRPS::compute_ari(data_pca,
                          cat_var=biological_subtypes)

    ################# Assessment on the batch effect ##################
    batch=sample.annot[, cat_var_label]
    ## PCA Color Batch
    colfunc <- colorRampPalette(brewer.pal(n = 4, name = 'Set1')[-6])
    color.batch <- colfunc(length(unique(batch)))
    names(color.batch) <- levels(batch)
    message("PCA based on Batch")
    ### Compute PCA Batch
    PCA_BATCH=RUVPRPS::plot_pca(data_pca,
                    variable= batch,
                    variable.name =   'Batch',
                    color =color.batch)

    ## Compute Silhouette based on batch
    message("Silhouette coefficient based on batch")
    silh_batch=RUVPRPS::compute_silhouette(data_pca,cat_var=batch)

    ## Compute ARI based on biology
    message("ARI based on batch")
    ari_batch=RUVPRPS::compute_ari(data_pca,cat_var=batch)


    ################## Assessment on the library size ##################
    library_size=sample.annot[, cont_var_label]
    ## Compute regression between library size and PCs
    message("Linear regression between the first cumulative PC and library size")
    reg_lib_size= RUVPRPS::regression_pc_contvar(pca=data_pca,
                               cont_var = library_size)

    ## Compute Spearman correlation between gene expression and library size
    message("Spearman correlation between individual gene expression and library size")
    corr_lib_size=RUVPRPS::correlation_gene_exp_contvar(se,
                                                        library_size,
                                                        apply.log)

    ## Plot combined silhouette based on batch and biology
    message("Combined silhouette plot")
    combined_silh_plot=RUVPRPS::plot_combined_silh_batch_bio(silh_bio,silh_batch)

    ################## Generate pdf file to save the plots #####################
    if (!is.null(output_file)){
        pdf(output_file)
            do.call(grid.arrange,
                c(PCA_BIO,
                  ncol = 4))
            plot(silh_bio$plot)
            plot(ari_bio$plot)
            do.call(grid.arrange,
                c(PCA_BATCH,
                  ncol = 4))
            plot(silh_batch$plot)
            plot(ari_batch$plot)
            plot(reg_lib_size$plot)
            plot(corr_lib_size$plot)
            plot(combined_silh_plot)
        dev.off()
    }
        res=list(PCA_bio=PCA_BIO,
                 silh_bio=silh_bio,
                 PCA_batch=PCA_BATCH,
                 silh_batch=silh_batch,
                 plot_reg_lib_size=reg_lib_size$plot,
                 plot_cor_gen_exp_lib_size=corr_lib_size$plot,
                 combined_silh_plot)
    return(res)
}
