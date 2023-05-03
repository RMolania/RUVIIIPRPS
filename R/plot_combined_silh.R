#' is used to plot a combined plot of the silhouette computed on two categorical variable
#' of a SummarizedExperiment class object (i.e. biology vs batch).
#'
#'
#' @param silh1 Vector of the Silhouette coefficient computed based on a categorical variable
#' @param silh2 Vector of the Silhouette coefficient based on batch.

#' @return plot Combined plot of the silhouette coefficient on the two categorical variables.
#' @importFrom wesanderson wes_palette
#' @import ggplot2
#' @export

plot_combined_silh<-function(
        silh1,
        silh2
){

## Merge two silh coeffs
df=cbind(silh1$silh.coeff,silh2$silh.coeff)
everything<-datasets<-Silh1<-Silh2<-NULL
colnames(df)=c("datasets","Silh1","datasets","Silh2")
df=df[,-3]

## Plot
dataSets.colors <- wes_palette(
        n = dim(df)[1],
        name = "GrandBudapest1")[seq(1:dim(df)[1])]
p=ggplot(df,aes(Silh1,Silh2,colour=datasets)) +
    geom_point(aes(color = datasets), size = 3) +
    xlab(paste('Silhouette based on ',silh1$cat_var_label,sep="")) +
    ylab (paste('Silhouette based on ',silh2$cat_var_label,sep="")) +
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
        legend.title = element_text(size = 16))+
    ggtitle("Combined silhouette coefficients")

return(p)

}
