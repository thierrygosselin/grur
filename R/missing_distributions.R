#' @title Missing data distributions
#' @description Graphical interface to choose parameters for missing data.
#' 
#' @return list of parameters 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} and 
#'   Thierry Gosselin \email{thierrygosselin@@icloud.com}
#'   
#' @export
#' 
missing_distributions <- function() {
  app.dir <- system.file("missing_distributions", package = "grur")
  shiny::runApp(app.dir, display.mode = "normal")
}
