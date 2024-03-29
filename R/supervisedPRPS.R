#' is used to create pseudo-replicates of pseudo samples of a SummarizedExperiment class object
#' using RUVIII-PRPS method with a supervised approach using the 'uv.variables', 'bio.variable' and 'batch.variable' given.
#'
#' We will create distinct group of pseudo-replicates for each source of unwanted variation defined in the 'uv.variables' argument.
#' For example to correct for batch effect if defined in the 'uv.variables' argument, several group of pseudo-samples
#' will be created by averaging the samples of the same biological subtype defined in 'bio.variable' in each batch. Then those
#' pseudo-samples will be defined as pseudo-replicates which constitutes a PRPS set.
#' For example to correct for library size if defined in the 'uv.variables' argument, several group of pseudo-samples
#' will be created by averaging the top and bottom-ranked samples by library size of the same biological subtype in each batch.
#' Then those pseudo-samples will be defined as pseudo-replicates which constitutes a PRPS set.
#' Similarly to correct for purity if defined in the 'uv.variables' argument, several group of pseudo-samples
#' will be created by averaging the top and bottom-ranked samples by purity of the same biological subtype in each batch.
#' Then those pseudo-samples will be defined as pseudo-replicates which constitutes a PRPS set.
#'
#'
#' @param se.obj A SummarizedExperiment object that will be used to create PRPS.
#' @param assay.name String for the selection of the name of the assay of the SummarizedExperiment class object.
#' @param bio.variable String of the label of a categorical variable that specifies major biological groups
#' such as samples types from colData(se).
#' @param uv.variables String or vector of strings of the label of continuous or categorical variable(s)
#' such as samples types, batch or library size from colData(se) that will be used to define PRPS.
#' @param batch.variable String of the label of a categorical variable that specifies major batch groups
#' such as plates from colData(se).
#' @param min.sample.for.prps Numeric. Indicates the minimum number of samples to create one pseudo-sample,
#' by default it is set to 3.
#' @param min.sample.per.batch Numeric. Indicates the minimum number of homogeneous biological samples within each batch
#' to create a PRPS set for each continuous variable. The minimum should be '2*min.sample.for.prps'. By default it is set to 6.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' By default it is set to TRUE.
#' @param assess.cor.variables Logical. Indicates whether to assess the association between pairs of categorical variables
#' and/or pairs of continuous variables in order to select only one of the variable of a pair of highly correlated
#' variables.
#' @param cat.cor.coef Vector of two numerical values. Indicates the cut-off of the correlation coefficient between each pair of
#' categorical variables. The first one is between each pair of 'uv.variables' and the second one is between each pair of 'bio.variables'.
#' The correlation is computed by the function ContCoef from the DescTools package. If the correlation of a pair of variable is higher than
#' the cut-off, then only the variable that has the highest number of factor will be kept and the other one will be excluded from the
#' remaining analysis. By default they are both set to 0.7.
#' @param cont.cor.coef Vector of two numerical values. Indicates the cut-off of the Spearman correlation coefficient between each pair of
#' continuous variables. The first one is between each pair of 'uv.variables' and the second one is between each pair of 'bio.variables'.
#' If the correlation of a pair of variable is higher than the cut-off, then only the variable that has the highest variance will
#' be kept and the other one will be excluded from the remaining analysis. By default they are both set to 0.7.
#' @param remove.na String. Indicates whether to remove NA or missing values from either the 'assays', the 'sample.annotation',
#' 'both' or 'none'. If 'assays' is selected, the genes that contains NA or missing values will be excluded. If 'sample.annotation' is selected, the
#' samples that contains NA or missing values for any 'bio.variables' and 'uv.variables' will be excluded. By default, it is set to
#' 'both'.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment
#' class object 'se.obj' or to output the result, by default it is set to TRUE.
#' @param plot.output Logical. Indicates whether to plot the PRPS map of the PRPS defined by the 'bio.variable'
#' and the categorical variables, by default it is set to TRUE.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed
#' during the execution of the functions, by default it is set to TRUE.

#' @return SummarizedExperiment or a List A SummarizedExperiment object or a list containing all the PRPS created.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr count
#' @importFrom tidyr %>%
#' @export

supervisedPRPS <- function(
        se.obj,
        assay.name,
        bio.variable,
        uv.variables,
        batch.variable,
        min.sample.for.prps = 3,
        min.sample.per.batch = 6,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        assess.cor.variables = FALSE,
        cat.cor.coef = c(0.7, 0.7),
        cont.cor.coef = c(0.7, 0.7),
        remove.na = 'both',
        save.se.obj = TRUE,
        plot.output=TRUE,
        verbose = TRUE
        ) {
    printColoredMessage(message = '------------The supervisedPRPS function starts.',
                          color = 'white',
                          verbose = verbose)
    ### Assess the input
    if (length(assay.name) > 1 || is.null(assay.name)) {
        stop('Please provide a single assay.name.')
    } else if (is.null(bio.variable)) {
        stop('The function requires some known biological groups.')
    } else if (is.null(uv.variables)) {
        stop('The function requires known sources of unwanted variation groups.')
    } else if (is.null(batch.variable)) {
        stop('The function requires a batch variable.')
    }

    ### Assess the correlation btw variables
    if(assess.cor.variables){
        se.obj <- variablesCorrelation(
            se.obj = se.obj,
            assay.name = assay.name,
            bio.variables = bio.variable,
            uv.variables = uv.variables,
            cat.cor.coef = cat.cor.coef,
            cont.cor.coef = cont.cor.coef,
            assess.se.obj = assess.se.obj,
            remove.na = remove.na,
            verbose = verbose)
        uv.variables <- se.obj$uv.variables
        bio.variable <- se.obj$bio.variable
        se.obj <- se.obj$se.obj
    }

    ### Define categorical and continuous variables
    uv.class <- sapply(
        uv.variables,
        function(x) class(colData(se.obj)[[x]])
         )
    categorical.uv <- names(uv.class[which(uv.class %in% c('character', 'factor'))])
    continuous.uv <- uv.variables[!uv.variables %in% categorical.uv]

    ### PRPS for categorical variables
    if(length(categorical.uv) > 0){
    printColoredMessage(
        message = '### Create PRPS for all categorical sources of unwanted variation.',
        color = 'blue',
        verbose = verbose
        )
        if(!save.se.obj){
            categorical.uv.prps <- lapply(
                categorical.uv,
                function(x){
                    prpsForCategoricalUV(
                        se.obj = se.obj,
                        assay.name = assay.name,
                        uv.variable = x,
                        bio.variable = bio.variable,
                        min.sample.for.prps = min.sample.for.prps,
                        apply.log = apply.log,
                        assess.se.obj = FALSE,
                        save.se.obj = save.se.obj,
                        pseudo.count = pseudo.count,
                        remove.na = remove.na,
                        verbose = verbose
                    )
                })
            names(categorical.uv.prps) <- categorical.uv
            categorical.uv.prps.all <- do.call(cbind, categorical.uv.prps)
        } else {
            for(x in categorical.uv){
                se.obj <- prpsForCategoricalUV(
                    se.obj = se.obj,
                    assay.name = assay.name,
                    uv.variable = x,
                    bio.variable = bio.variable,
                    min.sample.for.prps = min.sample.for.prps,
                    apply.log = apply.log,
                    assess.se.obj = FALSE,
                    save.se.obj = TRUE,
                    pseudo.count = pseudo.count,
                    remove.na = remove.na,
                    verbose = verbose
                )
            }
        }

    }
    if(length(continuous.uv) > 0){
    printColoredMessage(
        message = '### Create PRPS for all categorical sources of unwanted variation.',
        color = 'blue',
        verbose = verbose
        )
        if(!save.se.obj){
            continuous.uv.prps <- lapply(
                continuous.uv,
                function(x){
                    prpsForContinuousUV(
                        se.obj = se.obj,
                        assay.name = assay.name,
                        uv.variable = x,
                        bio.variable = bio.variable,
                        batch.variable = batch.variable,
                        min.sample.for.prps = min.sample.for.prps,
                        min.sample.per.batch=min.sample.per.batch,
                        apply.log = apply.log,
                        assess.se.obj = FALSE,
                        save.se.obj = save.se.obj,
                        pseudo.count = pseudo.count,
                        remove.na = remove.na,
                        verbose = verbose
                    )
                })
            names(continuous.uv.prps) <- continuous.uv
            continuous.uv.prps.all <- do.call(cbind, continuous.uv.prps)
        } else{
            for(x in continuous.uv){
                se.obj <- prpsForContinuousUV(
                    se.obj = se.obj,
                    assay.name = assay.name,
                    uv.variable = x,
                    bio.variable = bio.variable,
                    batch.variable = batch.variable,
                    min.sample.for.prps = min.sample.for.prps,
                    min.sample.per.batch=min.sample.per.batch,
                    apply.log = apply.log,
                    assess.se.obj = FALSE,
                    save.se.obj = TRUE,
                    pseudo.count = pseudo.count,
                    remove.na = remove.na,
                    verbose = verbose
                )
            }
        }
    }
    #### Ploting output
    if (plot.output==TRUE && (length(categorical.uv) > 0)){
        ## PRPS map plotting
        printColoredMessage(message= '### Plotting PRPS map for each uv categorical variable.',
                            color = 'magenta',
                            verbose = verbose)
        # Plot for each categorical variable
        catvar<-biology<-use<-n<-NULL
        plot.all<- lapply(
            categorical.uv,
            function(var){
                info <- as.data.frame(SummarizedExperiment::colData(se.obj))
                info$catvar <- as.factor(paste0(info[,var]))
                info$biology <- as.factor(paste0(info[,bio.variable]))
                df_count <- info %>%
                    count(catvar, biology)
                df_count$use <- 'unselected'
                df_count$use[df_count$n >= min.sample.for.prps] <- 'Selected'
                p=ggplot(df_count, aes(x = catvar, y = biology)) +
                    geom_count(aes(color = use)) +
                    geom_text(aes(
                        label = n,
                        hjust = 0.5,
                        vjust = 0.5
                    )) +
                    xlab(paste0(var)) +
                    ylab('Biological groups') +
                    theme_bw()+
                    theme(
                        axis.line = element_line(colour = 'black', size = .85),
                        axis.title.x = element_text(size = 18),
                        axis.title.y = element_text(size = 0),
                        axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
                        axis.text.y = element_text(size = 12, angle = 45, hjust = 1),
                        legend.position = 'none')
                p
            })
        names(plot.all)=categorical.uv

        ### Add plots to SummarizedExperiment object
        printColoredMessage(message= '### Saving the PRPS map plot to the metadata of the SummarizedExperiment object.',
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
        if(!'PRPS' %in% names(se.obj@metadata[['plot']]) ) {
            se.obj@metadata[['plot']][['PRPS']] <- list()
        }

        for (v in categorical.uv){
            ## Save the new plot
            se.obj@metadata[['plot']][['PRPS']][[v]]<-plot.all[[v]]
        }

    }

    ####### Save output
    if(save.se.obj){
        printColoredMessage(message = '------------The supervised.prps function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }else{
        printColoredMessage(message = '------------The supervised.prps function finished.',
                            color = 'white',
                            verbose = verbose)
        if(length(continuous.uv) > 0 & length(categorical.uv) > 0){
            return(list(
                categorical.uv.prps = categorical.uv.prps,
                continuous.uv.prps = continuous.uv.prps,
                all.prps = cbind(categorical.uv.prps.all, continuous.uv.prps.all))
            )
        } else if (length(continuous.uv) > 0 & length(categorical.uv) == 0){
            return(list(
                continuous.uv.prps = continuous.uv.prps,
                all.prps = continuous.uv.prps.all
            ))
        } else if (length(continuous.uv) == 0 & length(categorical.uv) > 0){
            return(list(
                categorical.uv.prps = categorical.uv.prps,
                all.prps = categorical.uv.prps.all
            ))
        }
    }

}





