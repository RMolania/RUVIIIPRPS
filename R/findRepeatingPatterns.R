#' is used to find repeating patterns
#'
#'
#' @param vector TO BE DEFINED
#' @param n TO BE DEFINED


findRepeatingPatterns <- function(vector, n) {
    char.counts <- table(vector)
    repeated.chars <- subset(char.counts, char.counts >= n)
    return(names(repeated.chars))
}
