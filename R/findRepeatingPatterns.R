#' is used to find repeating factors ro characters in a vector.
#'
#' @param vec A vector of factors or characters.
#' @param n.repeat Numeric. Indicates the minimum repeat of individual factors ro characters in the vector.
#'
#'
#' @return A vec of factors ro characters that are repeated at least "n.repeat" times.

findRepeatingPatterns <- function(vec, n.repeat) {
    char.counts <- table(vec)
    repeated.chars <- subset(char.counts, char.counts >= n.repeat)
    return(names(repeated.chars))
}
