#' is used to plot the RLE of a SummarizedExperiment class object.
#'
#' @param se A SummarizedExperiment object that will be used to compute the correlation
#' @param assay_names Optional string or list of strings for selection of the names
#' of the assays of the SummarizedExperiment class object to compute the correlation.
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom SummarizedExperiment assays assay
#' @importFrom kunstomverse geom_boxplot2
#' @importFrom matrixStats rowMedians colMedians colIQRs
#' @import ggplot2
#' @export

plot_RLE<-function(
        se,
        assay_names=NULL,
        apply.log=FALSE
){

    ## Assays
    if (!is.null(assay_names)){
        normalization=assay_names
    }else{
        normalization=names(assays(se))
    }

    ### RLE
    message('RLE plot')


    ### across all samples
    rle.all <- lapply(
        normalization,
        function(x){
            rle.comp <- function(expr) {
                rle.data <- expr - rowMedians(expr)
                rle.med <- colMedians(rle.data)
                rle.iqr <- colIQRs(rle.data)
                return(list(
                    rle = rle.data,
                    rle.med = rle.med,
                    rle.iqr = rle.iqr
                ))
            }

            ### log transformation
            if(apply.log){
                #message('Performing log + 1 transformation on the data')
                expr <- log2(assay(x = se, x) + 1)
            }else{
                expr <- assay(x = se, x)
            }
            rle=rle.comp(expr=expr)
            rle
        })
        names(rle.all)=normalization

        ## Plot RLE
        samples<-rle<-everything<-NULL
        plot.rle <- lapply(
            normalization,
            function(x){
                tmp=rle.all[[x]]$rle
                rle_plot=as.data.frame(tmp) %>% pivot_longer(everything(),
                    names_to = 'samples',values_to = 'rle') %>%
                    mutate(samples = factor(samples))
                p=ggplot(rle_plot, aes(x = samples,y=rle)) +
                geom_boxplot2(width=0.01,lwd=0.5) +
                ylab('RLE') +
                xlab('') +
                theme(panel.background = element_blank(),
                axis.line = element_line(colour = 'black'),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 8))+
                ggtitle(paste(" RLE plot distribution of ",x,sep=""))
            p
            })

    names(plot.rle) <- normalization
    return(list(rle=rle.all,plot=plot.rle))

}
