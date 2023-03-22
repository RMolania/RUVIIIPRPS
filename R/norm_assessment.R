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
        cat_var_label=NULL,
        cont_var_label=NULL,
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

    ### Check cat or cont var are selected
    if (is.null(cat_var_label) && is.null(cont_var_label)){
        stop('Please provide at least cat_var_label or cont_var_label')
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
    if (!is.null(cat_var_label)){
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
                    PCA=RUVPRPS::plot_pca(data_pca,
                                          assay_names = assay_names,
                                          cat_var=group,
                                          cat_var_label = x,
                                          color = color.group)

                    ## Compute Silhouette
                    message(paste("Silhouette coefficient based on: ",x,sep=""))
                    silh=RUVPRPS::compute_silhouette(data_pca,
                                                     assay_names = assay_names,
                                                     cat_var=group,
                                                     cat_var_label = x)

                    ## Compute ARI
                    message(paste("ARI based on: ",x,sep=""))
                    ari=RUVPRPS::compute_ari(data_pca,
                                             assay_names = assay_names,
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
    }

    ################# Continous variable ################
    if (!is.null(cont_var_label)){
        nb_cont_var=length(cont_var_label)
        cont.var.assessment<- lapply(
            cont_var_label,
            function(x){
                group=sample.annot[ , x]

                ## Compute regression between library size and PCs
                message("Linear regression between the first cumulative PC and library size")
                reg= RUVPRPS::regression_pc_contvar(pca=data_pca,
                                                    assay_names = assay_names,
                                                    cont_var = group,
                                                    cont_var_label=x)

                ## Compute Spearman correlation between gene expression and library size
                message("Spearman correlation between individual gene expression and library size")
                corr=RUVPRPS::correlation_gene_exp_contvar(se = se,
                                                           assay_names = assay_names,
                                                           cont_var = group,
                                                           cont_var_label=x,
                                                           apply.log)
                return(list(reg=reg,corr=corr))
            })
        names(cont.var.assessment)=cont_var_label
    }


    ################## Generate pdf file to save the plots #####################
    if (!is.null(output_file)){
        pdf(output_file)
        ## Categorical variable
            for (v in 1:(nb_cat_var)){
                plot(cat.var.assessment[[v]][['PCA']])
            }
            ## Combined silhouette
            p <- lapply(names(Combined_sil_plot),
                        function(x){
                        plot(Combined_sil_plot[[x]])
            })
        ## Continuous variable
        for (v in 1:(nb_cont_var)){
            plot(cont.var.assessment[[v]][['reg']][['plot']])
            plot(cont.var.assessment[[v]][['corr']][['plot']])
        }
        dev.off()
    }

        ### Results to return
        if (nb_cat_var>0){
        cat.var.ass=cat.var.assessment
        } else {
        cat.var.ass=NULL
        }
        if (nb_cont_var>0){
        cont.var.ass=cont.var.assessment
        } else {
        cont.var.ass=NULL
        }

        return(list(cat.var.ass=cat.var.ass,cont.var.ass=cont.var.ass))
}
