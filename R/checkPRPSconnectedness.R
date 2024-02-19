#' Assess the connectedness across PRPS sets.

#' @author Ramyar Molania

#' @description
#' A short description...

#' @param data.input TTTTT
#' @param min.samples TTTTTT
#' @param batch.name TTTTTT
#' @param verbose TTTTTT

#' @importFrom Matrix colSums rowSums


checkPRPSconnectedness <- function(
        data.input,
        min.samples,
        batch.name,
        verbose = TRUE) {
    colsum.data.input <- colSums(data.input >= min.samples)
    if (sum(colsum.data.input == 0) > 0) {
        printColoredMessage(
            message = paste0(
                sum(colsum.data.input == 0),
                ' batches of the ',
                batch.name ,
                ' do not have at least ',
                min.samples,
                ' samples across all the homogenous sample groups.'),
            color = 'red',
            verbose = verbose
        )
        stop(message = paste0(
                'So, no PRPS sets can be created for those batches. This may result in unsatisfactory removal of the ',
                batch.name,
                ' effects. Possibly, decreasing min.samples value can help.'))
    } else {
        selected.groups <- rowSums(data.input >= min.samples) > 1
        if (sum(selected.groups) == 0) {
            printColoredMessage(
                message =  paste0(
                    'There are not any homogenous sample groups that have at least ',
                    min.samples,
                    ' samples in more that one batch of.',
                    batch.name,
                    '.'),
                color = 'blue',
                verbose = verbose
            )
        } else{
            if (c(nrow(data.input) - sum(selected.groups)) > 0) {
                printColoredMessage(
                    message = paste0(
                        nrow(data.input) - sum(selected.groups),
                        ' homogenous sample groups are removed as they do not have at least ',
                        min.samples,
                        ' samples in more that one batch of ',
                        batch.name,
                        '.'
                    ),
                    color = 'blue',
                    verbose = verbose
                )
            }
            data.input <-
                data.input[selected.groups, , drop = FALSE]
            connection.check <- lapply(1:nrow(data.input),
                                       function(x) {
                                           groups.a <- names(which(data.input[x,] >= min.samples))
                                           sub.connections <-
                                               lapply(c(1:nrow(data.input))[-x],
                                                      function(y) {
                                                          groups.b <- names(which(data.input[y, ] >= min.samples))
                                                          inter.samples <-
                                                              intersect(groups.a, groups.b)
                                                          if (length(inter.samples) > 0) {
                                                              sort(unique(c(
                                                                  groups.a, groups.b
                                                              )), decreasing = FALSE)
                                                          } else
                                                              sort(groups.a, decreasing = FALSE)
                                                      })
                                           sort(unique(unlist(
                                               Filter(Negate(is.null), sub.connections)
                                           )), decreasing = FALSE)
                                       })
            covered.batches <-
                Filter(Negate(is.null), connection.check)
            covered.batches <-
                unlist(lapply(covered.batches, length))
            if (max(covered.batches) == ncol(data.input)) {
                printColoredMessage(
                    message = 'There is complete connection between the PRPS sets.',
                    color = 'blue',
                    verbose = verbose
                )
            } else {
                stop(
                    message = paste0(
                        'There is not any complete connection between possible PRPS sets. This may result in unsatifactoy removal of the ',
                        batch.name,
                        ' effects.'
                    )
                )
            }
        }
    }
    return(data.input)
}
