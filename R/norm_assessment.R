#' is used to assess the performance of the normalisation of a SummarizedExperiment class object.
#' Several assessment will be performed:
#' 1) PCA plot of each categorical variable.
#' 2) Silhouette and ARI computed on categorical variable.
#' 3) Combined Silhouette plot of the combined pair of categorical variable.
#' 4) Linear regression between the first cumulative PC and continuous variable.
#' 5) Spearman correlation between gene expression and continuous variable.
#'
#' @param se A SummarizedExperiment object that will be used to assess the performance of the normalisation of the data.
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
    } else if((!is.null(assay_names))&&!(assay_names %in% names(assays(se)))){
        stop('The selected assay(s) is/are not in the assays names of the SummarizedExperiment class object.')
    }

    ### Check cat_var_label and cont_var_label
    sample.annot <- as.data.frame(colData(x = se))
    all.var.label<- c(cat_var_label, cont_var_label)
    exist.var.label <- colnames(sample.annot)[colnames(sample.annot) %in% all.var.label]
    if (!sum(all.var.label %in% exist.var.label) == length(all.var.label)){
        stop(
            'Provided variable label from cat_var_label and cont_var_label "',
            paste0(all.var.label[!all.var.label %in% exist.var.label], collapse = ' & '),
            '" are not in the colData of the Summarized Experiment object.')
    }

    ### Compute PCA
    data_pca=RUVPRPS::compute_pca(se,apply.log = apply.log,assay_names = assay_names)

    ## Get all the available assays (i.e. normalizations methods)
    if (!is.null(assay_names)){
        normalization=assay_names
    }else{
        normalization=names(assays(se))
    }
    ################# Categorical variable ################
    nb_cat_var=length(cat_var_label)
    cat.var.assessment<- lapply(
        cat_var_label,
        function(x){
            group=as.factor(sample.annot[ , x])
            ## PCA Color
            colfunc <- colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = 'Spectral')[-6])
            color.group<- colfunc(length(unique(group)))
            names(color.group) <- levels(group)
            message(paste("PCA based on: ",x,sep=""))
            ### Compute PCA
            PCA=RUVPRPS::plot_pca(data_pca,assay_names = assay_names,
                                  cat_var=group,
                                  cat_var_label = x,
                                  color = color.group)

            ## Compute Silhouette
            message(paste("Silhouette coefficient based on: ",x,sep=""))
            silh=RUVPRPS::compute_silhouette(data_pca,assay_names = assay_names,
                                             cat_var=group,
                                             cat_var_label = x)

            ## Compute ARI
            message(paste("ARI based on: ",x,sep=""))
            ari=RUVPRPS::compute_ari(data_pca,assay_names = assay_names,
                                     cat_var=group,
                                     cat_var_label = x)
            return(list(PCA=PCA,sil=silh,ari=ari))
        })
    names(cat.var.assessment)=cat_var_label


    ## Plot combined silhouette based on all pairs of cat var
    Combined_sil_plot<-NULL
    if (nb_cat_var>1){
        message("Combined silhouette plot of cat var")
        for (v in 1:(nb_cat_var-1)){
            for (v2 in ((v+1):nb_cat_var)){
                p=RUVPRPS::plot_combined_silh(
                        cat.var.assessment[[cat_var_label[v]]][['sil']],
                        cat.var.assessment[[cat_var_label[v+1]]][['sil']])
                Combined_sil_plot[[paste0(cat_var_label[v],"_",cat_var_label[v2])]]=p
            }
        }
    }



    ################## Assessment on the library size ##################
    library_size=sample.annot[, cont_var_label]
    ## Compute regression between library size and PCs
    message("Linear regression between the first cumulative PC and library size")
    reg_lib_size= RUVPRPS::regression_pc_contvar(pca=data_pca,assay_names = assay_names,
                               cont_var = library_size)

    ## Compute Spearman correlation between gene expression and library size
    # message("Spearman correlation between individual gene expression and library size")
    # corr_lib_size=RUVPRPS::correlation_gene_exp_contvar(se,assay_names = assay_names,
    #                                                     library_size,
    #                                                     apply.log)


    ################## Generate pdf file to save the plots #####################
    if (!is.null(output_file)){
        pdf(output_file)
        ## Plot PCA
        for (v in 1:(nb_cat_var)){
            plot(cat.var.assessment[[v]][['PCA']])
        }
        p <- lapply(names(Combined_sil_plot), function(x) {
            plot(Combined_sil_plot[[x]])
        })
            plot(reg_lib_size$plot)

            #plot(corr_lib_size$plot)
            #plot(combined_silh_plot)
        dev.off()
    }
        res=list(#PCA_bio=PCA_BIO,
                 #silh_bio=silh_bio,
                 #PCA_batch=PCA_BATCH,
                 #silh_batch=silh_batch,
                 plot_reg_lib_size=reg_lib_size$plot)#,
                 #plot_cor_gen_exp_lib_size=corr_lib_size$plot)#,
                 #combined_silh_plot)
    return(res)
}
