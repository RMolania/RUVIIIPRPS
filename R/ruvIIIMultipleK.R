#' is used to run the RUVIII-PRPS method for a given assay and for multiple k values in one go.
#'
#'
#' @param se.obj A SummarizedExperiment object that will be used to computer fastRUV-III
#' @param assay.name String for the selection of the name of the assay data
#' of the SummarizedExperiment class object
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data, by default it is set to TRUE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements before log transformation,
#' by default it is set to 1.
#' @param replicate.data TO BE DEFINED
#' @param ctl Logical vector of length n of the negative control genes.
#' @param k A single value or a vector of values containing a single k or a range of k - the number of unwanted factors - to be tested.
#' @param eta A matrix with n columns, gene-wise (as opposed to sample-wise) covariates. By default is set to NULL.
#' @param include.intercept Logical. Add an intercept term to eta if it does not include one already. By default is set to TRUE.
#' @param apply.average.rep Average replicates after adjustment. By default is set to FALSE.
#' @param fullalpha Can be included to speed up execution. By default is set to NULL.
#' @param return.info If FALSE, only the adjusted data matrix is returned. If TRUE, additional information
#' is returned. By default is set to FALSE.
#' @param inputcheck Perform a basic sanity check on the inputs, and issue a warning if there is a problem.
#' By default is set to TRUE.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment class object.
#' @param remove.na TO BE DEFINED.
#' @param verbose Logical. Indicates whether to show or reduce the level of output or messages displayed during the execution
#' of the functions, by default it is set to TRUE.
#' @param save.se.obj Logical. Indicates whether to save the result in the metadata of the SummarizedExperiment class object 'se.obj' or
#' to output the result. By default it is set to TRUE.
#'
#'
#' @return list List containing the corrected gene expression matrix and the M, alpha and W for each k values tested.

#' @importFrom ruv RUV1 residop
#' @export

ruvIIIMultipleK <- function(
        se.obj,
        assay.name,
        apply.log=TRUE,
        pseudo.count = 1,
        replicate.data,
        ctl,
        k = NULL,
        eta = NULL,
        include.intercept = TRUE,
        apply.average.rep = FALSE,
        fullalpha = NULL,
        return.info = FALSE,
        inputcheck = TRUE,
        assess.se.obj = TRUE,
        remove.na = 'measurements',
        save.se.obj = TRUE,
        verbose = TRUE
){
    printColoredMessage(message = '------------The normalise function starts:',
                        color = 'white',
                        verbose = verbose)

    if(k == 0 || is.null(k)){
        stop('k cannot be 0. This means no adjustment will be made.')
    } else if(min(table(rownames(replicate.data))) == 1){
        stop('There are only replicated samples of a single sample in the replicate.data')
    }

    # check the SummarizedExperiment ####
    if (assess.se.obj) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose
        )
    }
    ## Run RUVIII_PRPS for all k value provided and save them into the se.obj
    if (save.se.obj) {
        for (x in k[1:length(k)]){
            se.obj<- ruvIII(
                se.obj = se.obj,
                assay.name=assay.name,
                apply.log=apply.log,
                pseudo.count = pseudo.count,
                replicate.data=replicate.data,
                ctl=ctl,
                k = x,
                eta = eta,
                include.intercept = include.intercept,
                apply.average.rep = apply.average.rep,
                fullalpha = fullalpha,
                return.info = return.info,
                inputcheck = inputcheck,
                assess.se.obj = assess.se.obj,
                remove.na = remove.na,
                save.se.obj = save.se.obj,
                verbose = verbose)
        }
        return(se.obj)

    ## Run RUVIII_PRPS for all k value provided and output the Y
    } else if (!return.info & !save.se.obj) {
        ## Run RUVIII_PRPS for the first k value provided
        Y.adj=ruvIII(
            se.obj = se.obj,
            assay.name=assay.name,
            apply.log=apply.log,
            pseudo.count = pseudo.count,
            replicate.data=replicate.data,
            k = k[1],
            eta = eta,
            include.intercept = include.intercept,
            apply.average.rep = apply.average.rep,
            fullalpha = fullalpha,
            return.info = return.info,
            inputcheck = inputcheck,
            assess.se.obj = assess.se.obj,
            remove.na = remove.na,
            save.se.obj = save.se.obj,
            verbose = verbose)
        ruv.adj=list(Y.adj)
        names(ruv.adj) <- paste0('RUVIII_K:', k[1], '_Data:', assay.name)

        ## if there are multiple k values provided
        if (length(k)>1) {
            ruv.adj.others.k <- lapply(
                k[2:length(k)],
                function(x){
                Y.adj.k<- ruvIII(
                    se.obj = se.obj,
                    assay.name=assay.name,
                    apply.log=apply.log,
                    pseudo.count = pseudo.count,
                    replicate.data=replicate.data,
                    k = x,
                    eta = eta,
                    include.intercept = include.intercept,
                    apply.average.rep = apply.average.rep,
                    fullalpha = fullalpha,
                    return.info = return.info,
                    inputcheck = inputcheck,
                    assess.se.obj = assess.se.obj,
                    remove.na = remove.na,
                    save.se.obj = save.se.obj,
                    verbose = verbose)
            })
            names(ruv.adj.others.k) <- paste0('RUV_K', k[2:length(k)], '_Data:', assay.name)
            ruv.adj.allk=append(ruv.adj.others.k, ruv.adj, after = 0)
            ruv.adj=ruv.adj.allk
        }
        ## Return a list containing the adjusted dataset(s) for single k or multiple k values
        return(ruv.adj)

    ## Run RUVIII_PRPS for all k value provided and output all the return info
    } else if (return.info & !save.se.obj) {
        return.info.k <- lapply(
            k,
            function(x){
                Y.adj.k<- ruvIII(
                    se.obj = se.obj,
                    assay.name=assay.name,
                    apply.log=apply.log,
                    pseudo.count = pseudo.count,
                    replicate.data=replicate.data,
                    k = x,
                    eta = eta,
                    include.intercept = include.intercept,
                    apply.average.rep = apply.average.rep,
                    fullalpha = fullalpha,
                    return.info = return.info,
                    inputcheck = inputcheck,
                    assess.se.obj = assess.se.obj,
                    remove.na = remove.na,
                    save.se.obj = save.se.obj,
                    verbose = verbose)
            })
            names(return.info.k) <- paste0('RUV_K', k[2:length(k)], '_Data:', assay.name)
            ## Return a list containing the info for the adjusted dataset(s) for single k or multiple k values
            return(return.info.k)
    }

    printColoredMessage(message = '------------The normalise function finished:',
                        color = 'white',
                        verbose = verbose)
}

