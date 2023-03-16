#' is used to assess the performance of the normalisation of a SummarizedExperiment class object.
#' It will first generate two PCA plots: one colored by biology and another one by batch,
#' each plot will display PCA plot for each assay in the following order: raw, fpkm, fpkm.uq and ruvprps.
#' It will also compute the regression between library size and PCs,
#'
#' @param se Dataset that will be used to assess the performance of the normalisation of the data.
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#' @param biological_subtypes Vector containing the biological subtype of each sample.
#' @param library_size Vector containing the library size of each sample.
#' @param batch Vector containing the batch (plates/years) of each sample.
#' @param output_file Path and name of the output file to save the assessments plots in a pdf format.
#' @param n.cores is the number of cpus used for mclapply parallelization. Default is set to 5.
#'
#'
#' @return plots List of assessments plots
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom gridExtra grid.arrange
#' @export

norm_assessment = function(
        se,
        apply.log = FALSE,
        biological_subtypes,
        library_size,
        batch,
        output_file=NULL,
        n.cores=5
){
    ### Compute PCA
    data_pca=RUVPRPS::compute_pca(se,apply.log = apply.log)

    ## Get all the available assays (i.e. normalizations methods)
    normalizations <- names(
        SummarizedExperiment::assays(se)
    )

    ################# Assessment on the biology ################
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
