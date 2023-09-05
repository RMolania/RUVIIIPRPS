
#' is used to plot the metric computed on a variable across the assays selected of a
#' SummarizedExperiment class object.
#'
#' @param se.obj A SummarizedExperiment object that will be used.
#' @param assay.names String or list of strings for the selection of the name
#' of the assays of the SummarizedExperiment class object.
#' @param metric Metric being assessed.
#' @param variable String of the label of a variable such as
#' sample types or batches from colData(se.obj).
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.

#' @return SummarizedExperiment A SummarizedExperiment object containing the associated plot.
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer %>%
#' @importFrom SummarizedExperiment assays assay colData
#' @importFrom fastDummies dummy_cols
#' @import ggplot2
#' @importFrom wesanderson wes_palette
#' @export

plotMetric <- function(
        se.obj,
        assay.names='All',
        metric=c('gene.pearson.corr','gene.spearman.corr','gene.aov.anova','gene.welch.correction.anova','ari','sil','pcs.vect.corr','pcs.lm'),
        variable,
        verbose=TRUE
        ) {

    printColoredMessage(message = '------------The plotMetric function starts:',
                        color = 'white',
                        verbose = verbose)

    if (!metric %in% c('gene.pearson.corr','gene.spearman.corr','gene.aov.anova','gene.welch.correction.anova','ari','sil','pcs.vect.corr','pcs.lm')) {
        stop(paste0('The ', metric,'is not suitable. It has to be selected from the gene.pearson.corr, gene.spearman.corr, gene.aov.anova,
                gene.welch.correction.anova,ari, sil, pcs.vect.corr, pcs.lm metrics.'))
    }


    ## Assays
    if(length(assay.names) == 1 && assay.names=='All'){
        assay.names=as.factor(names(assays(se.obj)))
    }else {
        assay.names=as.factor(unlist(assay.names))
    }

    all.assays.metric <- lapply(
        levels(assay.names),
        function(x){
            if (!variable %in% names(se.obj@metadata[['metric']][[x]][[metric]])) {
                stop(paste0('The ', metric,'has not been computed yet for the ',variable,' variable and the ',x, ' assay.'))
            }
            se.obj@metadata[['metric']][[x]][[metric]][[variable]]
        })
    names(all.assays.metric)=levels(assay.names)
    everything<-datasets<-corr.coeff<-pcs<-NULL
    all.assays.metric <- as.data.frame(all.assays.metric)

    if (metric %in% c('pcs.vect.corr')){
        nb.pca.comp=10
        all.assays.metric= all.assays.metric %>% mutate(pcs=c(1:nb.pca.comp)) %>% pivot_longer(
            -pcs,
            names_to = 'datasets',
            values_to = {{metric}}) %>% mutate(datasets = factor(
                datasets,levels=levels(assay.names)))
    } else if (metric %in% c('pcs.lm')){
        nb.assays=length(assay.names)
        all.assays.metric=cbind(all.assays.metric,pcs=(1:10))
        all.assays.metric= all.assays.metric %>% pivot_longer(
            -(nb.assays+1),
            names_to = 'datasets',
            values_to = {{metric}}) %>% mutate(datasets = factor(
                datasets,levels=levels(assay.names)))
    } else {
        all.assays.metric= all.assays.metric %>% pivot_longer(
        everything(),
        names_to = 'datasets',
        values_to = {{metric}}) %>% mutate(datasets = factor(
            datasets,levels=levels(assay.names)))
    }


    p=ggplot(all.assays.metric, aes_string(x = 'datasets',y=metric, fill = 'datasets'))

    if (metric %in% c('gene.pearson.corr','gene.spearman.corr')){
        p=p+ geom_boxplot()
        xlabel=''
    } else if (metric %in% c('gene.aov.anova','gene.welch.correction.anova')) {
        p=p+ geom_boxplot2()
        xlabel=''
    } else if (metric %in% c('ari','sil')){
        p=p+ geom_point(aes(colour=datasets))
        xlabel=''
    } else if (metric %in% c('pcs.vect.corr','pcs.lm')){
        p=p+geom_line(aes(color = datasets), size = 1) +
            geom_point(aes(color = datasets), size = 3)
        xlabel='PCs'
    }
    p= p +
        ylab(paste0(metric)) +
        xlab(xlabel) +
        scale_y_continuous()+
        geom_hline(yintercept=0)+
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', size = 1),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12))+
        ggtitle(paste0("Assessment: ",metric," computed on ",variable))
    if  (length(levels(assay.names)) <= 4){
        dataSets.colors <- wes_palette(
            n = length(levels(assay.names)),
            name = "GrandBudapest1")[seq(1:length(levels(assay.names)))]
        p=p+scale_fill_manual(values = dataSets.colors, guide = 'none')
    }


    ### Add plots to SummarizedExperiment object
    printColoredMessage(message= '### Saving the plot to the metadata of the SummarizedExperiment object.',
                        color = 'magenta',
                        verbose = verbose)
    ## Check if metadata plot already exist
    if(length(se.obj@metadata)==0 ) {
            se.obj@metadata[['plot']] <- list()
    }
    ## Check if metadata plot already exist
    if(!'plot' %in% names(se.obj@metadata) ) {
        se.obj@metadata[['plot']] <- list()
    }
    ## Check if metadata plot already exist for this metric
    if(!metric %in% names(se.obj@metadata[['plot']]) ) {
            se.obj@metadata[['plot']][[metric]] <- list()
    }
    ## Save the new plot
    se.obj@metadata[['plot']][[metric]][[variable]]<- list()
    se.obj@metadata[['plot']][[metric]][[variable]] <- p

    printColoredMessage(message = '------------The plotMetric function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}

