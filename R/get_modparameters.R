#' Gives description of the parameters of modeling functions
#'
#' @import magrittr
#' @import dplyr
#'
#' @param fun character
#' @param parameter character or character vector
#'
#' @export
#'
#' @examples
#' get_modparameters(fun = "mod_popoccup")
#' get_modparameters(parameter = "z")
#' get_modparameters(parameter = c("z","N","B"))
get_modparameters <- function(fun = NULL, parameter = NULL) {
  if (!is.null(fun)) {
    modfunparameters <- filter(modfunparameters, Function %in% fun)
  }
  if (!is.null(parameter)) {
    modfunparameters <- filter(modfunparameters, Para %in% parameter)
  }
  modfunparameters <- select(modfunparameters, -Para, - Function) %>% distinct() %>% print(right = FALSE, row.names = FALSE)
}
