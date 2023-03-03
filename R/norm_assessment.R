#' is used to assess the performance of the normalisation of the data
#'
#'
#' @param sce Dataset that will be used to assess the performance of the normalisation of the data.
#' It will generate two PCA plots: one colored by biology and another one by time, each plot will display PCA plot for each assay in the following order: raw, fpkm, fpkm.uq and ruvprps.
#' @param apply.log Indicates whether to apply a log-transformation to the data
#' @param biological_subtypes Vector containing the biological subtypes of each sample
#' @param library_size Vector containing the library size of each sample
#' @param time Vector containing the time (plates/years) of each sample
#' @param output_file Path and name of the output file to save the assessments plots in a pdf format
#' @param var_da_library_size Vector containing the categorical variable of high vs low library size of each sample to compute differential analysis
#' @param n.cores is the number of cpus used for mclapply parallelization
#'
#'
#' @return plots List of assessments plots
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom gridExtra grid.arrange
#' @export

## put as options the different variable
## batch instead of time
norm_assessment = function(
        sce,
        apply.log = FALSE,
        biological_subtypes,
        library_size,
        time,
        output_file=NULL,
        var_da_library_size=NULL,
        n.cores=5
){
    ### Compute PCA
    data_pca=RUVPRPS::compute_pca(sce,apply.log = apply.log)

    ## Get all the available assays (i.e. normalizations methods)
    normalizations <- names(
        SummarizedExperiment::assays(sce)
    )

    ### Assessment on the biology ####
    # Color Biology
    colfunc <- colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = 'Spectral')[-6])
    color.subtype<- colfunc(length(unique(biological_subtypes)))
    names(color.subtype) <- levels(biological_subtypes)
    message("PCA based on Biology")
    ### Compute PCA Biology
    pp_bio <- lapply(
        normalizations,
        function(x){
            pcs <- data_pca[[x]]
            p1 <- RUVPRPS::pca_plot_squared(
                pca = pcs,
                variable= biological_subtypes,
                variable.name =  'Biology',
                color = color.subtype)
            p1
        })
    names(pp_bio) <- normalizations
    plot_BIO=c(pp_bio[[1]],
               pp_bio[[2]],
               pp_bio[[3]],
               pp_bio[[4]])

    ### Assessment on the time effect ####
    # Color Time (years)
    colfunc <- colorRampPalette(brewer.pal(n = 4, name = 'Set1')[-6])
    color.time <- colfunc(length(unique(time)))
    names(color.time) <- levels(time)
    message("PCA based on Time")
    ### Compute PCA Time
    pp_time <- lapply(
        normalizations,
        function(x){
            pcs <- data_pca[[x]]
            p1 <- RUVPRPS::pca_plot_squared(
                pca = pcs,
                variable= time,
                variable.name =  'Time',
                color = color.time)
            p1
        })
    names(pp_time) <- normalizations
    plot_TIME=c(pp_time[[1]],
                pp_time[[2]],
                pp_time[[3]],
                pp_time[[4]])


    # ### Assessment on the library size ####
    ## Compute regression between library size and PCs
    ## Regression on the library size
    message("Regression based on Library size")
    reg_lib_size= RUVPRPS::regression_pc(pca=data_pca,
                               normalization=normalizations,
                               regression_var=library_size)

    ## DE between sample with low and high library size
    if (!is.null(var_da_library_size)){
    de_analysis_lib_size=de_analysis(sce,
        var_da_library_size,
        apply.log)
    }

    #### Generate pdf file to save the plots
    if (!is.null(output_file)){
        pdf(output_file)
            do.call(grid.arrange,
                c(plot_BIO,
                  ncol = 4))
            do.call(grid.arrange,
                c(plot_TIME,
                  ncol = 4))
            plot(reg_lib_size$plot)
            if (!is.null(var_da_library_size)){
                plot(de_analysis_lib_size$plot)
                }
        dev.off()
    }
    if (!is.null(var_da_library_size)){
        res=list(plot_bio=plot_BIO,
                 plot_time=plot_TIME,
                 plot_reg_lib_size=reg_lib_size$plot,
                 plot_de_analysis_lib_size=de_analysis_lib_size$plot)
    }else{
        res=list(plot_bio=plot_BIO,
                 plot_time=plot_TIME,
                 plot_reg_lib_size=reg_lib_size$plot)
        }
    return(res)
}
