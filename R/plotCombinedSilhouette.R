#' is used to plot a combined plot of the silhouette computed on two categorical variables
#' of a SummarizedExperiment class object (i.e. biology vs batch).
#'
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the differential expression ANOVA analysis.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the correlation. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param variable1 String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param variable2 String of the label of a categorical variable such as
#' sample types or batches from colData(se.obj).
#' @param method A character string indicating which method
#' is to be used for the differential analysis: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.
#' By default 'euclidean' will be selected.
#' @param save.se.obj Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#' @param assess.se.obj Indicates whether to assess the SummarizedExperiment class object.
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.

#' @return plot Combined plot of the silhouette coefficient on the two categorical variables.
#' @importFrom wesanderson wes_palette
#' @import ggplot2
#' @export


plotCombinedSilhouette<-function(
                    se.obj,
                    assay.names='All',
                    variable1,
                    variable2,
                    method='euclidean',
                    save.se.obj = TRUE,
                    assess.se.obj = TRUE,
                    verbose=TRUE
){
    ## variables = c(variable1,variable2)
    if (is.null(assay.names)) {
        stop('Please provide at least an assay name.')
    } else if (class(se.obj@colData[, variable1]) %in% c('numeric', 'integer')) {
    stop(paste0('The ', variable1,', is a numeric, but this should a categorical variable'))
    }


    printColoredMessage(message = '------------The plotCombinedSilhouette function starts:',
                        color = 'white',
                        verbose = verbose)

    ### check the inputs
    if (is.null(assay.names)) {
        stop('Please provide at least an assay name.')
    }

    ### Assess the se.obj
    if (assess.se.obj){
        se.obj <- checkSeObj(se.obj = se.obj,
                             assay.names = assay.names,
                             variables = c(variable1,variable2),
                             remove.na = 'both',
                             verbose = verbose)
    }

    ## Assays
    if(length(assay.names) == 1 && assay.names=='All'){
        assay.names=as.factor(names(assays(se.obj)))
    }else {
        assay.names=as.factor(unlist(assay.names))
    }

    ## Categorical variables
    var1=se.obj@colData[, variable1]
    var1.label=paste0(variable1)
    var2=se.obj@colData[, variable2]
    var2.label=paste0(variable2)

    ## Merge two silh coeffs
    all.assays.var1 <- lapply(
        levels(assay.names),
        function(x){
            if (!variable1 %in% names(se.obj@metadata[['metric']][[x]][[paste0('sil.',method)]])) {
                stop(paste0('The Silhouette based on ',method,' has not been computed yet for the ',variable1,' variable and the ',x, ' assay.'))
            }
            se.obj@metadata[['metric']][[x]][[paste0('sil.',method)]][[variable1]]
        })
    names(all.assays.var1)=levels(assay.names)
    all.assays.var1 <- as.data.frame(all.assays.var1)
    all.assays.var1= all.assays.var1 %>% pivot_longer(
        everything(),
        names_to = 'datasets',
        values_to = 'sil') %>% mutate(datasets = factor(
            datasets,levels=levels(assay.names)))
    all.assays.var2 <- lapply(
        levels(assay.names),
        function(x){
            if (!variable2 %in% names(se.obj@metadata[['metric']][[x]][[paste0('sil.',method)]])) {
                stop(paste0('The Silhouette based on ',method,' has not been computed yet for the ',variable2,' variable and the ',x, ' assay.'))
            }
            se.obj@metadata[['metric']][[x]][[paste0('sil.',method)]][[variable2]]
        })
    names(all.assays.var2)=levels(assay.names)
    all.assays.var2 <- as.data.frame(all.assays.var2)
    all.assays.var2= all.assays.var2 %>% pivot_longer(
        everything(),
        names_to = 'datasets',
        values_to = 'sil') %>% mutate(datasets = factor(
            datasets,levels=levels(assay.names)))

    df=cbind(all.assays.var1,all.assays.var2)
    everything<-datasets<-Silh1<-Silh2<-NULL
    colnames(df)=c("datasets","Silh1","datasets","Silh2")
    df=df[,-3]

    ## Plot
    p=ggplot(df,aes(Silh1,Silh2,colour=datasets)) +
        geom_point(aes(color = datasets), size = 3) +
        xlab(paste('Silhouette based on ',method,' computed on ',var1.label,sep="")) +
        ylab (paste('Silhouette based on ',method,' computed on ',var2.label,sep="")) +
        theme(
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black', size = 1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16))+
        ggtitle("Combined silhouette coefficients")

    if  (length(levels(assay.names)) <= 4){
        dataSets.colors <- wes_palette(
        n = length(levels(assay.names)),
        name = "GrandBudapest1")[seq(1:length(levels(assay.names)))]
        p=p+scale_color_manual(
            values=c(dataSets.colors),
            name = 'Datasets')
    }

    printColoredMessage(message = '------------The plotCombinedSilhouette function finished.',
                    color = 'white',
                    verbose = verbose)
return(p)

}
