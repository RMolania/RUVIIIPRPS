#' is used to identify the unwanted variation of a SummarizedExperiment class object.
#' Several assessment will be performed:
#' 1) PCA plot of each categorical variable.
#' 2) Silhouette and ARI computed on categorical variable.
#' 3) Combined Silhouette plot of the combined pair of categorical variable.
#' 4) Linear regression between the first cumulative PC and continuous variable.
#' 5) Spearman correlation between gene expression and continuous variable.
#'
#' @param se A SummarizedExperiment object that will be used to assess the performance of the normalisation of the data.
#' @param assay_names String or list of strings for selection of the name
#' of the assay of the SummarizedExperiment class object.
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

identify_unwanted_variation = function(
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
        stop('Please provide a summarized experiment object.\n')
    } else if (class(se)[1] == 'SummarizedExperiment' & ncol(colData(se)) == 0){
        stop('The Summarized experiment object does not contain
           sample annotations, please provide colData of the summarized
           experiment object before running the norm_assessment() function.\n')
    } else if((!is.null(assay_names))&&!(assay_names %in% names(assays(se)))){
        stop('The selected assay is not in the assays names of the SummarizedExperiment class object.\n')
    }

    ### Check cat_var_label and cont_var_label
    sample.annot <- as.data.frame(colData(x = se))
    all.var.label<- c(cat_var_label, cont_var_label)
    exist.var.label <- colnames(sample.annot)[colnames(sample.annot) %in% all.var.label]
    if (!sum(all.var.label %in% exist.var.label) == length(all.var.label)){
        stop(
            'Provided variable label from cat_var_label and cont_var_label "',
            paste0(all.var.label[!all.var.label %in% exist.var.label], collapse = ' & '),
            '" are not in the colData of the Summarized Experiment object.\n')
    }

    ### Check cat or cont var are selected
    if (is.null(cat_var_label) && is.null(cont_var_label)){
        stop('Please provide at least cat_var_label or cont_var_label.\n')
    }


    assessment=RUVPRPS::norm_assessment(se=se,assay_names = assay_names,apply.log = apply.log,
                                        cat_var_label = cat_var_label,cont_var_label =cont_var_label,
                                        output_file =  output_file)
    return(assessment)
}
