#' Converts quantitative variables into qualitative variables
#'
#' @import magrittr
#' @import dplyr
#' @importFrom rlang enquo
#'
#' @param df data frame
#' @param var_quali argument
#'
#' @return df data frame
#' @export
#'
#' @examples
#' data(riverfish)
#' mef_convertquali(riverfish, var_quali="flow")
mef_convertquali <- function(df, var_quali) {
  var_quali <- enquo(var_quali)
  name_quali <- colnames(select(df, !!var_quali))
  for (i in name_quali) {
    modality <- unique(pull(df, i)[!is.na(pull(df, i))])
    for (j in 2:length(modality)) {
      df <- df %>% mutate(across(paste(i), ~ifelse(is.na(.x), NA, ifelse(.x %in% modality[j], 1, 0)), .names = paste(i,modality[j],sep="_")))
      }
    }
  return(df)
}
