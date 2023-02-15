#' is used to assess the performance of the normalisation of the data
#'
#'
#' @param sce Dataset that will be used to assess the performance of the normalisation of the data
#' @param apply.log Indicates whether to apply a log-transformation to the data
#' @param biological_subtypes Vector containing the biological subtypes of each sample
#' @param library_size Vector containing the library size of each sample
#' @param time Vector containing the time (plates/years) of each sample
#' @param output_file Path and name of the output file to save the assessments plots in a pdf format
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
        time,
        output_file=NULL
){
    ### Compute PCA
    data_pca=RUVPRPS::compute_pca(sce,apply.log = apply.log)

    ## Get all the available normalizations methods
    normalizations <- names(
        SummarizedExperiment::assays(sce)
    )

    ### Assessment on the biology ####
    # Color Biology
    colfunc <- colorRampPalette(brewer.pal(n = 11, name = 'Spectral')[-6])
    message("Test print")
    message(unique(biological_subtypes))
    color.subtype<- colfunc(length(unique(biological_subtypes)))
    print(color.subtype)
    message("Test print 2")
    message(color.subtype)
    names(color.subtype) <- levels(biological_subtypes)
    print(color.subtype)
    message("Test print 3")
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
    print(pp_bio[[1]])
    message(pp_bio[[1]])
    pp_bio[[1]]
    plot_BIO=do.call(
        grid.arrange,
        c(pp_bio[[1]],
          pp_bio[[2]],
          pp_bio[[3]],
          pp_bio[[4]],
          ncol = 4))

    # ### Assessment on the time effect ####
    # # Color Time (years)
    # colfunc <- colorRampPalette(brewer.pal(n = 4, name = 'Set1')[-6])
    # color.time <- colfunc(length(unique(time)))
    # print(color.time)
    # names(color.time) <- levels(time)
    # print(color.time)
    # ### Compute PCA Time
    # pp_time <- lapply(
    #     normalizations,
    #     function(x){
    #         pcs <- data_pca[[x]]
    #         p1 <- RUVPRPS::pca_plot_squared(
    #             pca = pcs,
    #             variable= time,
    #             variable.name =  'Time',
    #             color = color.time)
    #         p1
    #     })
    # names(pp_time) <- normalizations
    # plot_TIME=do.call(
    #     grid.arrange,
    #     c(pp_time[[1]],
    #       pp_time[[2]],
    #       pp_time[[3]],
    #       pp_time[[4]],
    #       ncol = 4))
    #
    #
    # ### Assessment on the library size ####
    #
    # if (!is.null(output_file)){
    #     pdf(output_file)
    #     plot_BIO
    #     plot_TIME
    #     dev.off()
    # }
    # return(list(plot_bio=plot_BIO,plot_time=plot_TIME))
    return(list(plot_bio=plot_BIO,plot_time="t2"))
}
