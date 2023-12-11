#' is used to print colored messages.
#'
#'
#' @param message Message to be printed.
#' @param color The color of the message.
#' @param verbose whether to show the messages or not.

############################################# print_colored_message
printColoredMessage <- function(message, color, verbose) {
    # ANSI escape sequence for color
    colors <- switch(
        tolower(color),
        "red"     = "\033[31m",
        "green"   = "\033[32m",
        "yellow"  = "\033[33m",
        "blue"    = "\033[34m",
        "magenta" = "\033[35m",
        "cyan"    = "\033[36m",
        "white"   = "\033[37m",
        "reset"   = "\033[0m",
        ""
    )
    # Print the colored message
    if(verbose){
        cat(paste0(colors, message, '\n'))
    } else {
        cat('')
    }
}
