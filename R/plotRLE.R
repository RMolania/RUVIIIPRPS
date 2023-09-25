#' is used to plot the RLE (Relative Log Expression) of a SummarizedExperiment class object.
#'
#' @param se.obj A SummarizedExperiment object that will be used to plot the RLE.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to plot RLE. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#' @param save.se.obj Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param assess.se.obj Indicates whether to assess the SummarizedExperiment class object.
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#'
#' @return list List of the computed RLE and the associated plot
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom SummarizedExperiment assays assay
#' @importFrom kunstomverse geom_boxplot2
#' @importFrom matrixStats rowMedians colMedians colIQRs
#' @importFrom stats median
#' @import ggplot2
#' @export

plotRLE<-function(
        se.obj,
        assay.names="All",
        apply.log=FALSE,
        pseudo.count = 1,
        save.se.obj = TRUE,
        assess.se.obj = TRUE,
        verbose = TRUE
){
    printColoredMessage(message = '------------The plotRLE function starts:',
                        color = 'white',
                        verbose = verbose)

    ### Assess the se.obj
    if(assess.se.obj){
        se.obj <- checkSeObj(se.obj = se.obj,
                             assay.names = assay.names,
                             variables = NULL,
                             remove.na = 'measurements',
                             verbose = verbose)
    }

    ## Assays
    if(length(assay.names) == 1 && assay.names=='All'){
        assay.names=as.factor(names(assays(se.obj)))
    }else {
        assay.names=as.factor(unlist(assay.names))
    }

    ### Compute RLE across all Assays
    rle.all <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(message = paste0(
                '### Computing RLE on the ',
                x,
                ' assay.'
            ),
            color = 'magenta',
            verbose = verbose)
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
                expr <- log2(assay(x = se.obj, x) + pseudo.count)
            }else{
                expr <- assay(x = se.obj, x)
            }
            rle=rle.comp(expr=expr)
            rle
        })
        names(rle.all)=levels(assay.names)

        ## Plot RLE
        samples<-rle<-everything<-NULL
        plot.rle <- lapply(
            levels(assay.names),
            function(x){
                tmp=rle.all[[x]]$rle
                rle_plot=as.data.frame(tmp) %>% pivot_longer(everything(),
                    names_to = 'samples',values_to = 'rle') %>%
                    mutate(samples = factor(samples))
                p=ggplot(rle_plot, aes(x = samples,y=rle)) +
                #geom_boxplot2(width.errorbar =0.01,outlier.alpha=0.2) +
                    geom_boxplot(width=0.01,lwd=0.5,alpha=0.2)
                ylab('RLE') +
                xlab('') +
                coord_cartesian(ylim=c(-6,6))+
                stat_summary(geom = "crossbar", width=5, fatten=8, color="darkgreen",
                                 fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})+
                theme(panel.background = element_blank(),
                axis.line = element_line(colour = 'black'),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 8))+
                ggtitle(paste(" RLE plot distribution of ",x,sep=""))+
                    geom_hline(yintercept=0)
            p
        })
    names(plot.rle) <- levels(assay.names)


    ### Add results to the SummarizedExperiment object
    if(save.se.obj == TRUE){
        printColoredMessage(message= '### Saving the RLE plot to the metadata of the SummarizedExperiment object.',
                            color = 'magenta',
                            verbose = verbose)
        ## For all assays
        for (x in levels(assay.names)){
            ## Check if metadata plot already exist
            if(length(se.obj@metadata)==0 ) {
                se.obj@metadata[['plot']] <- list()
            }
            ## Check if metadata plot already exist
            if(!'plot' %in% names(se.obj@metadata) ) {
                se.obj@metadata[['plot']] <- list()
            }
            ## Check if metadata plot already exist for this metric
            if(!'rle' %in% names(se.obj@metadata[['plot']]) ) {
                se.obj@metadata[['plot']][['rle']] <- list()
            }
            ## Check if metadata plot already exist for this metric
            if(!x %in% names(se.obj@metadata[['plot']]) ) {
                ## Save the new plot
                se.obj@metadata[['plot']][['rle']][[x]] <- plot.rle[[x]]
            }
        }
        printColoredMessage(message= paste0(
            'The RLE plots are saved to metadata@$plot$rle'),
            color = 'blue',
            verbose = verbose)
        return(se.obj = se.obj)
    } else if(save.se.obj == FALSE){
        return(plot=plot.rle)
    }

        printColoredMessage(message = '------------The plotRLE function finished.',
                            color = 'white',
                            verbose = verbose)
}
