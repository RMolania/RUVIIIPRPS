#' is used to assess the performance of the normalisation of the data
#'
#'
#' @param sce Dataset that will be used to assess the performance of the normalisation of the data.
#' It will first generate two PCA plots: one colored by biology and another one by batch,
#' each plot will display PCA plot for each assay in the following order: raw, fpkm, fpkm.uq and ruvprps.
#' It will also compute the regression between library size and PCs,
#' and if asked it will compute a differential analysis between sample with low and high library size (i.e early years vs late years)
#' @param apply.log Indicates whether to apply a log-transformation to the data
#' @param biological_subtypes Vector containing the biological subtypes of each sample
#' @param library_size Vector containing the library size of each sample
#' @param batch Vector containing the batch (plates/years) of each sample
#' @param output_file Path and name of the output file to save the assessments plots in a pdf format
#' @param catvar_da_library_size Vector containing the categorical variable of high vs low library size of each sample to compute differential analysis
#' @param n.cores is the number of cpus used for mclapply parallelization
#'
#'
#' @return plots List of assessments plots
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom gridExtra grid.arrange
#' @export

norm_assessment = function(
        sce,
        apply.log = FALSE,
        biological_subtypes,
        library_size,
        batch,
        output_file=NULL,
        catvar_da_library_size=NULL,
        n.cores=5
){
    ### Compute PCA
    data_pca=RUVPRPS::compute_pca(sce,apply.log = apply.log)

    ## Get all the available assays (i.e. normalizations methods)
    normalizations <- names(
        SummarizedExperiment::assays(sce)
    )

    ################# Assessment on the biology ################
    ## PCA Color Biology
    colfunc <- colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = 'Spectral')[-6])
    color.subtype<- colfunc(length(unique(biological_subtypes)))
    names(color.subtype) <- levels(biological_subtypes)
    message("PCA based on Biology")
    ### Compute PCA Biology
    # pp_bio <- lapply(
    #     normalizations,
    #     function(x){
    #         pcs <- data_pca[[x]]
    #         p1 <- RUVPRPS::pca_plot(
    #             pca = pcs,
    #             variable= biological_subtypes,
    #             variable.name =  'Biology',
    #             color = color.subtype)
    #         p1
    #     })
    # names(pp_bio) <- normalizations
    PCA_BIO=RUVPRPS::plot_pca(data_pca,
                    variable= biological_subtypes,
                    variable.name =  'Biology',
                     color = color.subtype)
    # PCA_BIO=c(pp_bio[[1]],
    #            pp_bio[[2]],
    #            pp_bio[[3]],
    #            pp_bio[[4]])

    ## Compute Silhouette based on biology
    message("Silhouette coefficient based on biology")
    silh_bio=RUVPRPS::compute_silhouette(data_pca,
                          normalizations,
                          cat_var=biological_subtypes)

    ## Compute ARI based on biology
    message("ARI based on biology")
    ari_bio=RUVPRPS::compute_ari(data_pca,
                          normalizations,
                          cat_var=biological_subtypes)

    ################# Assessment on the batch effect ##################
    ## PCA Color Batch
    colfunc <- colorRampPalette(brewer.pal(n = 4, name = 'Set1')[-6])
    color.batch <- colfunc(length(unique(batch)))
    names(color.batch) <- levels(batch)
    message("PCA based on Batch")
    ### Compute PCA Batch
    # pp_batch <- lapply(
    #     normalizations,
    #     function(x){
    #         pcs <- data_pca[[x]]
    #         p1 <- RUVPRPS::pca_plot(
    #             pca = pcs,
    #             variable= batch,
    #             variable.name =  'Batch',
    #             color = color.batch)
    #         p1
    #     })
    # names(pp_batch) <- normalizations
    PCA_BATCH=RUVPRPS::plot_pca(data_pca,
                    variable= batch,
                    variable.name =   'Batch',
                    color =color.batch)
    # PCA_BATCH=c(pp_batch[[1]],
    #             pp_batch[[2]],
    #             pp_batch[[3]],
    #             pp_batch[[4]])

    ## Compute Silhouette based on batch
    message("Silhouette coefficient based on batch")
    silh_batch=RUVPRPS::compute_silhouette(data_pca,
                            cat_var=batch)

    ## Compute ARI based on biology
    message("ARI based on batch")
    ari_batch=RUVPRPS::compute_ari(data_pca,cat_var=batch)

    ## Plot combined silhouette based on batch and biology
    #combined_silh=


    ################## Assessment on the library size ##################
    ## Compute regression between library size and PCs
    message("Linear regression between the first cumulative PC and library size")
    reg_lib_size= RUVPRPS::regression_pc_contvar(pca=data_pca,
                               normalization=normalizations,
                               cont_var = library_size)

    ## Compute Spearman correlation between gene expression and library size
    message("Spearman correlation between individual gene expression and library size")
    corr_lib_size=RUVPRPS::correlation_gene_exp_contvar(sce,
                                                        library_size,
                                                        apply.log)

    # ## DA between sample with low and high library size
    # if (!is.null(catvar_da_library_size)){
    #     message("Differential analysis using Wilcoxon test between samples with high vs low library size")
    #     da_analysis_lib_size=RUVPRPS::da_analysis_wilcoxon_gene_exp_catvar_all_assays(sce,
    #                                                                                   catvar_da_library_size,
    #                                                                                   apply.log)
    # }

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
            # if (!is.null(catvar_da_library_size)){
            #     plot(da_analysis_lib_size$plot)
            # }
        dev.off()
    }
    # if (!is.null(catvar_da_library_size)){
    #     res=list(PCA_bio=PCA_BIO,
    #              silh_bio=silh_bio,
    #              ari_bio=ari_bio,
    #              PCA_batch=PCA_BATCH,
    #              silh_batch=silh_batch,
    #              ari_batch=ari_batch,
    #              plot_reg_lib_size=reg_lib_size$plot,
    #              plot_cor_gen_exp_lib_size=corr_lib_size$plot,
    #              plot_da_analysis_lib_size=da_analysis_lib_size$plot)
    # }else{
        res=list(PCA_bio=PCA_BIO,
                 silh_bio=silh_bio,
                 PCA_batch=PCA_BATCH,
                 silh_batch=silh_batch,
                 plot_reg_lib_size=reg_lib_size$plot,
                 plot_cor_gen_exp_lib_size=corr_lib_size$plot)
        #}
    return(res)
}
