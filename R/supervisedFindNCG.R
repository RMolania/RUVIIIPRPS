#' is used to find a set of negative control genes of a SummarizedExperiment class object
#' using RUVIII-PRPS method with a supervised approach using the 'bio.variables' and
#' 'uv.variables' given.
#'
#' @param se.obj A summarized experiment object.
#' @param assay.name String for the selection of the name of the assay of the SummarizedExperiment class object used to define NCG.
#' If you provide the raw data assay, we recommend to set apply.normalization to TRUE.
#' @param bio.variables String or vector of strings of the label of a categorical variable that specifies major biological groups
#' such as samples types from colData(se) that will be used to find the negative controls.
#' @param uv.variables String or vector of strings of the label of continuous or categorical variable(s)
#' such as samples types, batch or library size from colData(se) that will be used to find the negative controls.
#' @param no.ncg Logical, TO BE BETTER DEFINED. if TRUE then a sample annotation the initially contains column names of the assays.???
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param regress.out.uv.variables TO BE DEFINED.
#' @param regress.out.bio.variables TO BE DEFINED.
#' @param apply.normalization Logical Indicates whether to apply a normalization method when providing the raw data assay.
#' By default it is set to FALSE.
#' @param normalization String defining the normalization method to use from 'CPM', 'TMM', 'upper', 'full', 'median', 'VST',
#' and 'Quantile'. By default it is set to 'CPM'.
#' @param corr.method TO BE DEFINED.
#' @param a TO BE DEFINED.
#' @param rho TO BE DEFINED.
#' @param anova.method TO BE DEFINED.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' By default it is set to TRUE.
#' @param assess.variables Logical. TO BE DEFINED.
#' @param remove.na TO BE DEFINED.
#' @param plot.output Logical. Indicates whether to plot the PCA of the defined NCG colored by the categorical variables, by default it is set to TRUE.
#' @param fast.pca Logical. Indicates whether to calculate a specific number of PCs instead of the full range
#' to speed up the process, by default is set to 'TRUE'.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment
#' class object 'se.obj' or to output the result, by default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed
#' during the execution of the functions, by default it is set to TRUE.

#' @return SummarizedExperiment or a List A SummarizedExperiment object or a list containing
#' all the negative control genes defined.
#' @importFrom matrixTests row_oneway_equalvar
#' @importFrom ggpubr as_ggplot
#' @importFrom SummarizedExperiment assay
#' @importFrom stats as.formula
#' @importFrom Rfast correls
#' @importFrom graphics frame
#' @export

supervisedFindNGC <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        no.ncg = 1000,
        apply.log = TRUE,
        pseudo.count = 1,
        regress.out.uv.variables = FALSE,
        regress.out.bio.variables = FALSE,
        apply.normalization=FALSE,
        normalization = 'CPM',
        corr.method = "pearson",
        a = 0.05,
        rho = 0,
        anova.method = 'aov',
        assess.se.obj = TRUE,
        assess.variables = TRUE,
        remove.na = 'both',
        plot.output=TRUE,
        fast.pca = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The supervisedFindNGC function starts:',
                        color = 'white',
                        verbose = verbose)

    # check the SummarizedExperiment ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(bio.variables, uv.variables),
            remove.na = remove.na,
            verbose = verbose
        )
    }
    # check the variables ####
    if (assess.variables) {
        se.obj <- variablesCorrelation(
            se.obj = se.obj,
            bio.variables = bio.variables,
            uv.variables = uv.variables,
            cont.coef= c(0.95, 0.95),
            spearman.coef = c(0.95, 0.95),
            assess.se.obj = TRUE,
            remove.na = 'none',
            verbose = verbose
        )
        bio.variables <- se.obj$bio.variables
        uv.variables <- se.obj$uv.variables
        se.obj <- se.obj$se.obj
    }
    # data transformation ####
    if(apply.log){
        expr.data <- log2(assay(se.obj, assay.name) + pseudo.count)
    }else{
        expr.data <- assay(se.obj, assay.name)
    }
    ####
    if(isTRUE(regress.out.uv.variables)){
        expr.data.reg.uv <- transpose(expr.data)
        uv.variables.all <- paste('se.obj', uv.variables, sep = '$')
        y <- lm(as.formula(paste(
            'expr.data.reg.uv',
            paste0(uv.variables.all, collapse = '+') ,
            sep = '~'
        )))
        expr.data.reg.uv <- transpose(y$residuals)
        colnames(expr.data.reg.uv) <- colnames(se.obj)
        row.names(expr.data.reg.uv) <- row.names(se.obj)
    } else if(isTRUE(regress.out.bio.variables)){
        expr.data.reg.bio <- transpose(expr.data)
        bio.variables.all <- paste('se.obj', bio.variables, sep = '$')
        y <- lm(as.formula(paste(
            'expr.data.reg.bio',
            paste0(bio.variables.all, collapse = '+') ,
            sep = '~'
        )))
        expr.data.reg.bio <- transpose(y$residuals)
        colnames(expr.data.reg.bio) <- colnames(se.obj)
        row.names(expr.data.reg.bio) <- row.names(se.obj)
    }
    if(isTRUE(apply.normalization)){
        expr.data.nor <- applyOtherNormalizations(
            se.obj = se.obj,
            assay.name = assay.name,
            method = normalization,
            pseudo.count = pseudo.count,
            apply.log = apply.log,
            assess.se.obj = FALSE,
            save.se.obj=FALSE,
            remove.na = 'none',
            verbose = verbose
        )
    }

    # finding negative control genes ####
    ## step1: highly affected by unwanted variation ####
    printColoredMessage(message = '### Finding genes that are highly affected by unwanted variation:',
                        color = 'magenta',
                        verbose = verbose)
    uv.var.class <- unlist(lapply(uv.variables,
                                  function(x) {
                                      class(colData(se.obj)[[x]])
                                  }))
    categorical.uv <- uv.variables[uv.var.class %in% c('factor', 'character')]
    continuous.uv <- uv.variables[uv.var.class %in% c('numeric', 'integer')]
    all.tests <- list()
    ind=1
    ###
    if(!is.null(categorical.uv)){
        if(isTRUE(regress.out.bio.variables)){
            data.to.use <- expr.data.reg.bio
        } else{
            data.to.use <- expr.data
        }
        ### gene-batch anova
        printColoredMessage(
            message = paste0(
                'Performing ANOVA between individual gene ',
                'expression and each categorical source of unwanted variation:',
                paste0(categorical.uv, collapse = '&'),
                '.'
            ),
            color = 'blue',
            verbose = verbose
        )
        anova.gene.uv <- lapply(
            categorical.uv,
            function(x) {
                freq <- table(colData(se.obj)[[x]])
                bio.exc <- names(which(freq < 3))
                temp.data.a <- as.data.frame(data.to.use[,!colData(se.obj)[[x]] %in% bio.exc])
                anova.gene.batch <- as.data.frame(
                    row_oneway_equalvar(
                        x = temp.data.a,
                        g = se.obj@colData[, x]))
                anova.gene.batch$ranked.genes <- rank(anova.gene.batch[, 'statistic'])
                anova.gene.batch
            })
        names(anova.gene.uv) <- categorical.uv
        rm(data.to.use)
        all.tests[[ind]] <- anova.gene.uv
        names(all.tests) <- ('anova.gene.uv')
        ind=ind+1
    }
    if(!is.null(continuous.uv)){
        if(isTRUE(regress.out.bio.variables)){
            data.to.use <- expr.data.reg.bio
        } else{
            data.to.use <- expr.data
        }
        printColoredMessage(
            message = paste0(
                'Performing Spearman correlation between individual gene ',
                'expression and each continuous source of unwanted variation:',
                paste0(continuous.uv, collapse = '&'),
                '.'
            ),
            color = 'blue',
            verbose = verbose
        )
        corr.genes.uv <- lapply(
            continuous.uv,
            function(x) {
                corr.genes.var <- as.data.frame(correls(
                    y = se.obj@colData[, x],
                    x = transpose(data.to.use),
                    type = corr.method,
                    a = a ,
                    rho = rho
                ))
                corr.genes.var <- cbind(round(x = corr.genes.var[, 1:4], digits = 2),
                          corr.genes.var[, 5, drop = FALSE])
                corr.genes.var$ranked.genes <- rank(abs(corr.genes.var[, 'correlation']))
                row.names(corr.genes.var) <- row.names(data.to.use)
                corr.genes.var
            })
        names(corr.genes.uv) <- continuous.uv
        rm(data.to.use)
        all.tests[[ind]] <-corr.genes.uv
        names(all.tests)[[ind]] <- 'corr.genes.uv'
        ind=ind+1
    }

    ## step2: not highly affected by biology ####
    printColoredMessage(message = '### Finding genes that are not highly affected by biology',
                        color = 'magenta',
                        verbose = verbose)
    bio.var.class <- unlist(lapply(
        bio.variables,
        function(x) {
            class(colData(se.obj)[[x]])
        }))
    continuous.bio <- bio.variables[bio.var.class %in% c('numeric', 'integer')]
    categorical.bio <- bio.variables[bio.var.class %in% c('factor', 'character')]
    ###
    if(!is.null(continuous.bio)){
        if(isTRUE(apply.normalization)){
            data.to.use <- expr.data.nor
        } else if(isTRUE(regress.out.uv.variables)){
            data.to.use <- expr.data.reg.uv
        } else{
            data.to.use <- expr.data
        }
        ### gene-batch anova
        printColoredMessage(
            message = paste0(
                'Performing ANOVA between individual gene ',
                'expression and each categorical source of unwanted variation:',
                paste0(categorical.uv, collapse = '&'),
                '.'
            ),
            color = 'blue',
            verbose = verbose
        )
        corr.genes.bio <- lapply(
            continuous.bio,
            function(x) {
                corr.genes.var <- as.data.frame(correls(
                    y = se.obj@colData[, x],
                    x = transpose(data.to.use),
                    type = corr.method,
                    a = a ,
                    rho = rho
                ))
                corr.genes.var <- cbind(round(x = corr.genes.var[, 1:4], digits = 2),
                          corr.genes.var[, 5, drop = FALSE])
                row.names(corr.genes.var) <- row.names(data.to.use)
                corr.genes.var$ranked.genes <- rank(abs(corr.genes.var[, 'correlation']))
                corr.genes.var
            })
        names(corr.genes.bio) <- continuous.bio
        all.tests[[ind]] <-corr.genes.bio
        names(all.tests)[[ind]] <- 'corr.genes.bio'
        ind=ind+1
    }
    if(!is.null(continuous.uv)){
        if(isTRUE(normalization)){
            data.to.use <- expr.data.nor
        } else if(isTRUE(regress.out.uv.variables)){
            data.to.use <- expr.data.reg.uv
        } else{
            data.to.use <- expr.data
        }
        printColoredMessage(
            message = paste0(
                'Performing Spearman correlation between individual gene ',
                'expression and each continuous source of unwanted variation:',
                paste0(continuous.uv, collapse = '&'),
                '.'
            ),
            color = 'blue',
            verbose = verbose
        )
        anova.gene.bio <- lapply(
            categorical.bio,
            function(x) {
                freq <- table(colData(se.obj)[[x]])
                bio.exc <- names(which(freq < 3))
                se.obj.temp <- se.obj[,!colData(se.obj)[[x]] %in% bio.exc]
                temp.data <- assay(se.obj.temp, assay.name)
                anova.gene <- as.data.frame(
                    row_oneway_equalvar(x = data.to.use,
                                        g = se.obj@colData[, x]))
                anova.gene$ranked.genes <- rank(anova.gene[, 'statistic'])
                anova.gene
            })
        names(anova.gene.bio) <- categorical.bio
        all.tests[[ind]] <-anova.gene.bio
        names(all.tests)[[ind]] <- 'anova.gene.bio'
    }
    # step3: final selection ####
    #all.tests <- c('anova.gene.bio', 'anova.gene.uv', 'corr.genes.bio', 'corr.genes.uv')
    ncg.selected <- lapply(
        names(all.tests),
        function(x){
                temp <- all.tests[[x]]
                #if(!is.null(temp)){
                    ranks.data <- lapply(
                        names(temp),
                        function(y){
                            temp[[y]]$ranked.genes
                        })
                    do.call(cbind, ranks.data)
                #}
        })
    ncg.selected <- do.call(cbind, ncg.selected)
    ncg.selected <- rank(-apply(ncg.selected, 1, prod))
    ncg.selected <- ncg.selected < no.ncg + 1


    #### Ploting output
    cat.var=c(categorical.bio,categorical.uv)
    if (plot.output==TRUE && !is.null(cat.var)){
        se.ncg.selected=se[ncg.selected,]
        ### Compute PCA
        if (fast.pca) {
                se.ncg.selected=RUVPRPS::computePCA(se.obj=se.ncg.selected,
                                           assay.names = assay.name,
                                           apply.log = apply.log,
                                           pseudo.count = pseudo.count,
                                           fast.pca = fast.pca,
                                           nb.pcs = 10,
                                           assess.se.obj = FALSE,
                                           verbose = verbose)

        } else {
                se.ncg.selected=RUVPRPS::computePCA(se.obj=se.ncg.selected,
                                           assay.names = assay.name,
                                           apply.log = apply.log,
                                           pseudo.count = pseudo.count,
                                           fast.pca = fast.pca,
                                           nb.pcs = 10,
                                           assess.se.obj = FALSE,
                                           verbose = verbose)
        }

        ## PCA plotting
        printColoredMessage(message= '### Plotting PCA of the NCG colored by the categorical variables.',
                            color = 'magenta',
                            verbose = verbose)

        PCA.plots<- lapply(
            cat.var,
                function(x){
                    ## PCA Color
                    group=as.factor(se.ncg.selected@colData[, x])
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
                    PCA=RUVPRPS::plotPCA(se.obj=se.ncg.selected,
                                         assay.names = assay.name,
                                         variable=x,
                                         color = color.group,
                                         fast.pca=fast.pca,
                                         assess.se.obj = assess.se.obj,
                                         verbose = verbose)
                    return(PCA)
        })
        names(PCA.plots)=cat.var
        plot=do.call(grid.arrange,c(PCA.plots))

        ### Add plots to SummarizedExperiment object
        printColoredMessage(message= '### Saving the NCG plot to the metadata of the SummarizedExperiment object.',
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
        if(!'NCG' %in% names(se.obj@metadata[['plot']]) ) {
            se.obj@metadata[['plot']][['NCG']] <- list()
        }
        ## Save the new plot
        se.obj@metadata[['plot']][['NCG']]<-as_ggplot(plot)

    }


    ### Add results to the SummarizedExperiment object
    if(save.se.obj == TRUE){
        printColoredMessage(message= '### Saving the NCG vector to the metadata of the SummarizedExperiment object.',
                            color = 'magenta',
                            verbose = verbose)

        ## Check if metadata NCG already exists
        if(length(se.obj@metadata)==0 ) {
            se.obj@metadata[['NCG']] <- list()
        }
        ## Assign NCG
        se.obj@metadata[['NCG']] <- ncg.selected

        printColoredMessage(message=
            'The NCG are saved to metadata@NCG.',
            color = 'blue',
            verbose = verbose)
        return(se.obj)
    } else{
        return(ncg.selected)
    }

    printColoredMessage(message = '------------The supervisedFindNGC function finished.',
                        color = 'white',
                        verbose = verbose)
}
