#' is used to assess the performance of the normalisation of a SummarizedExperiment class object.
#'
#' Several assessment will be performed:
#' For each categorical variable:
#' - PCA plot of the categorical variable.
#' - Silhouette and ARI computed on the categorical variable.
#' - Differential analysis based ANOVA between the gene expression and the categorical variable.
#' - Vector correlation between the first cumulative PCs of the gene expression and the categorical variable.
#' For each continous variable:
#' - Linear regression between the first cumulative PC and continuous variable.
#' - Correlation between gene expression and continuous variable.
#'
#' It will output the following plots:
#' - PCA plot of each categorical variable.
#' - Boxplot of the F-test distribution from ANOVA between the gene expression and each categorical variable.
#' - Vector correlation between the first cumulative PCs of the gene expression and each categorical variable.
#' - Combined Silhouette plot of the combined pair of all categorical variables.
#' - Linear regression between the first cumulative PC and continuous variable.
#' - Boxplot of the correlation between gene expression and continuous variable.
#' - It will also output the RLE plot distribution.
#'
#' @param se.obj A SummarizedExperiment object that will be used to compute the PCA.
#' @param assay.names Optional string or list of strings for the selection of the name(s)
#' of the assay(s) of the SummarizedExperiment class object to compute the correlation. By default
#  all the assays of the SummarizedExperiment class object will be selected.
#' @param bio.variable String of the label of a categorical variable, from colData(se.obj) that specify major biological groups.
#' @param uv.variables String or vector of strings of the label of continuous or categorical variable(s)
#' such as samples types, batch or library size from colData(se).
#' @param apply.log Indicates whether to apply a log-transformation to the data. By default
#' no transformation will be selected.
#' @param output_file Path and name of the output file to save the assessments plots in a pdf format.
#' @param fast.pca logical. Indicates whether to calculate a specific number of PCs instead of the full range
#' to speed up the process, by default is set to 'TRUE'.
#' @param nb.pcs Numeric. The number of first PCs to be calculated for the fast pca process, by default is set to 10.
#' @param assess.se.obj Indicates whether to assess the SummarizedExperiment class object.
#' @param verbose Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#'
#'
#' @return list List of assessments plots and metrics used for the assessment
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom gridExtra grid.arrange
#' @importFrom SummarizedExperiment assays colData
#' @export

## remove n.core and merge all variable together

normAssessment = function(
        se.obj,
        assay.names = 'All',
        apply.log = TRUE,
        bio.variable=NULL,
        uv.variables=NULL,
        output_file=NULL,
        fast.pca = TRUE,
        nb.pcs = 10,
        assess.se.obj = TRUE,
        verbose = TRUE,
        pseudo.count = 1
){
    printColoredMessage(message = '------------The normAssessment function starts.',
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

    ### check the inputs of PCA
    if (fast.pca & is.null(nb.pcs)) {
        stop('To perform fast PCA, the number of PCs (nb.pcs) must specified.')
    } else if (fast.pca & nb.pcs == 0) {
        stop('To perform fast PCA, the number of PCs (nb.pcs) must specified.')
    }


    ### Categorical, continuous, biological variables
    if (!is.null(uv.variables)){
        uv.class <- sapply(
            uv.variables,
            function(x) class(colData(se.obj)[[x]])
        )
        categorical.uv <- names(uv.class[which(uv.class %in% c('character', 'factor'))])
        continuous.uv <- uv.variables[!uv.variables %in% categorical.uv]
        if(!is.null(bio.variable)){
            categorical.uv <-c(bio.variable,categorical.uv)
        }
    } else if(!is.null(bio.variable)){
                categorical.uv <-bio.variable
                continuous.uv <-NULL
    } else {
        categorical.uv <-NULL
        continuous.uv <-NULL
    }

    ## Assays
    if(length(assay.names) == 1 && assay.names=='All'){
        assay.names=as.factor(names(assays(se.obj)))
    }else {
        assay.names=as.factor(unlist(assay.names))
    }


    ### Compute PCA
    printColoredMessage(message = paste0(
        '### Computing PCA.'
    ),
    color = 'magenta',
    verbose = verbose)
    se.obj=RUVPRPS::computePCA(se.obj=se.obj,
                                 assay.names = assay.names,
                                 apply.log = apply.log,
                                 pseudo.count = pseudo.count,
                                 fast.pca = fast.pca,
                                 nb.pcs = nb.pcs,
                                 assess.se.obj = assess.se.obj,
                                 verbose = verbose)

    ################# Categorical variable ################
    if (!is.null(categorical.uv)){
        nb.cat.var=length(categorical.uv)
        ## PCA plotting
        PCA.plots<- lapply(
            categorical.uv,
            function(x){
                ## PCA Color
                group=as.factor(se.obj@colData[, x])
                if (length(unique(group))<=11){
                    colfunc <- colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = 'Spectral')[-6])
                    color.group<- colfunc(length(unique(group)))
                    names(color.group) <- unique(group)
                } else {color.group=NULL}

                ### Plot PCA
                printColoredMessage(message = paste0(
                    '### Plotting PCA based on the ',
                    x,
                    ' variable.'
                ),
                color = 'magenta',
                verbose = verbose)
                PCA=RUVPRPS::plotPCA(se.obj=se.obj,
                                     assay.names = assay.names,
                                     variable=x,
                                     color = color.group,
                                     fast.pca=fast.pca,
                                     nb.pcs=nb.pcs,
                                     assess.se.obj = assess.se.obj,
                                     verbose = verbose)
                return(PCA)
        })
        names(PCA.plots)=categorical.uv

        ## Computing other metrics
        for (x in categorical.uv){
                ## Compute Silhouette
                printColoredMessage(message = paste0(
                    '### Computing Silhouette based on the ',
                    x,
                    ' variable.'
                ),
                color = 'magenta',
                verbose = verbose)
                se.obj=RUVPRPS::computeSilhouette(se.obj=se.obj,
                                                assay.names = assay.names,
                                                variable=x,
                                                fast.pca=fast.pca,
                                                assess.se.obj = assess.se.obj,
                                                verbose = verbose)

                ## Compute ARI
                printColoredMessage(message = paste0(
                    '### Computing ARI based on the ',
                    x,
                    ' variable.'
                ),
                color = 'magenta',
                verbose = verbose)
                se.obj=RUVPRPS::computeARI(se.obj=se.obj,
                                        assay.names = assay.names,
                                        variable=x,
                                        fast.pca=fast.pca,
                                        assess.se.obj = assess.se.obj,
                                        verbose = verbose)

                ## Compute ANOVA
                printColoredMessage(message = paste0(
                    '### Computing ANOVA based on the ',
                    x,
                    ' variable.'
                ),
                color = 'magenta',
                verbose = verbose)
                se.obj=RUVPRPS::genesVariableAnova(se.obj=se.obj,
                                                  assay.names = assay.names,
                                                  variable=x,
                                                  apply.log=apply.log,
                                                  pseudo.count = pseudo.count,
                                                  assess.se.obj = assess.se.obj,
                                                  verbose = verbose)

                ## Compute Vector correlation
                printColoredMessage(message = paste0(
                    '### Computing Vector Correlation between the first cumulative PCs and the ',
                    x,
                    ' variable.'
                ),
                color = 'magenta',
                verbose = verbose)
                se.obj=RUVPRPS::PCVariableCorrelation(se.obj=se.obj,
                                                    assay.names = assay.names,
                                                    variable=x,
                                                    fast.pca=fast.pca,
                                                    nb.pcs=nb.pcs,
                                                    assess.se.obj = assess.se.obj,
                                                    verbose = verbose)
        }

        ## Plot combined silhouette based on all pairs of cat var
        CombinedSilPlot<-NULL
        if (nb.cat.var>1){
            printColoredMessage(message = paste0(
                '### Plotting all combined silhouette plots of all categorical variables'),
            color = 'magenta',
            verbose = verbose)
            for (v in 1:(nb.cat.var-1)){
                for (v2 in ((v+1):nb.cat.var)){
                    p=RUVPRPS::plotCombinedSilhouette(se.obj=se.obj,
                                                      assay.names = assay.names,
                                                      variable1=categorical.uv[v],
                                                      variable2=categorical.uv[v2],
                                                      assess.se.obj = assess.se.obj,
                                                      verbose = verbose)
                    CombinedSilPlot[[paste0(categorical.uv[v],"_",categorical.uv[v2])]]=p
                }
            }
        }
    }

    ################# Continous variable ################
    if (!is.null(continuous.uv)){
        nb.cont.var=length(continuous.uv)
        ## Computing other metrics
        for (x in continuous.uv){
            ## Compute regression between library size and PCs
            printColoredMessage(message = paste0(
                '### Computing Linear regression between the first cumulative PCs and the ',
                x,
                ' variable.'
            ),
            color = 'magenta',
            verbose = verbose)
            se.obj=RUVPRPS::PCVariableRegression(se.obj,
                                                 assay.names = 'All',
                                                 variable=x,
                                                 fast.pca=fast.pca,
                                                 nb.pcs=nb.pcs,
                                                 assess.se.obj = assess.se.obj,
                                                 verbose = verbose)

            ## Compute Spearman correlation between gene expression and library size
            printColoredMessage(message = paste0(
                '### Computing Spearman correlation based on',
                x,
                ' variable.'
            ),
            color = 'magenta',
            verbose = verbose)
            se.obj=RUVPRPS::genesVariableCorrelation(se.obj,
                                                     assay.names = 'All',
                                                     variable=x,
                                                     apply.log=apply.log,
                                                     pseudo.count = pseudo.count,
                                                     assess.se.obj = assess.se.obj,
                                                     verbose = verbose)
        }
    }

    ########## RLE plot ############
    # RLE
    rle=RUVPRPS::plotRLE(se.obj=se.obj,
                          assay.names = assay.names,
                          apply.log=apply.log)

    ################## Generate pdf file to save the plots #####################
    if (!is.null(output_file)){
        printColoredMessage(message= paste0(
            'The plots are being saved into the output file'),
            color = 'blue',
            verbose = verbose)
        pdf(output_file)
        ## Categorical variable
        if (!is.null(categorical.uv)){
            for (v in 1:(nb.cat.var)){
                plot(PCA.plots[[v]])
                plot(se.obj@metadata[['plot']][['gene.aov.anova']][[categorical.uv[v]]])
                plot(se.obj@metadata[['plot']][['pcs.vect.corr']][[categorical.uv[v]]])
            }
            ## Combined silhouette
            p <- lapply(names(CombinedSilPlot),
                        function(x){
                            plot(CombinedSilPlot[[x]])
                        })
        }
        ## Continuous variable
        if (!is.null(continuous.uv)){
            for (v in 1:(nb.cont.var)){
                plot(se.obj@metadata[['plot']][['pcs.lm']][[categorical.uv[v]]])
                plot(se.obj@metadata[['plot']][['gene.spearman.corr']][[categorical.uv[v]]])
            }
        }
        ## RLE plot
        lreg.pcs<- lapply(
            levels(assay.names),
            function(x){
                plot(rle$plot[[x]])
            })
        dev.off()
    }
    printColoredMessage(message = '------------The normAssessment function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}
