decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(abs(x - round(x)))-2
  } else {
    return(0)
  }
}
