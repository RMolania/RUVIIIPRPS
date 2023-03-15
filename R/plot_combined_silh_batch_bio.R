#' is used to plot a combined plot of the silhouette computed on batch versus the one
#' computed on biology of a SummarizedExperiment class object.
#'
#'
#' @param silh_bio Vector of the Silhouette coefficient computed based on biology.
#' @param silh_batch Vector of the Silhouette coefficient based on batch.

#' @return Combined plot of the silhouette computed on batch versus the one
#' computed on biology
#' @importFrom wesanderson wes_palette
#' @import ggplot2
#' @export

plot_combined_silh_batch_bio<-function(
        silh_bio,
        silh_batch
){

## Merge two silh coeffs
df=cbind(silh_bio,silh_batch)
everything<-datasets<-SilhBio<-SilhBatch<-NULL
colnames(df)=c("datasets","SilhBio","datasets","SilhBatch")
df=df[,-3]

## Plot
dataSets.colors <- wes_palette(
        n = dim(df)[1],
        name = "GrandBudapest1")[c(1,2,4,3)]
p=ggplot(df,aes(SilhBio,SilhBatch,colour=datasets)) +
    geom_point(aes(color = datasets), size = 3) +
    xlab('Silhouette based on Biology') +
    ylab ('Silhouette based on Bacth') +
    scale_color_manual(
        values=c(dataSets.colors),
        name = 'Datasets') +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black', size = 1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16))
return(p)

}
