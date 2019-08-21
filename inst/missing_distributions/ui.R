ui <- shiny::fluidPage(
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::h4("Select Pr(missing)"),
      shiny::h4("Loci"),
      shiny::fluidRow(
        column(6, shiny::numericInput("loc.shape1", "Beta v", min = 0, value = 0.2, step = 0.05)),
        column(6, shiny::numericInput("loc.shape2", "Beta w", min = 0, value = 1, step = 0.05))
      ),
      shiny::h4("Individuals"),
      shiny::fluidRow(
        column(6, shiny::numericInput("ind.shape1", "Beta v", min = 0, value = 0.2, step = 0.05)),
        column(6, shiny::numericInput("ind.shape2", "Beta w", min = 0, value = 1, step = 0.05))
      ),
      shiny::h4("Strata"),
      shiny::fluidRow(
        column(6, shiny::numericInput("strata.shape1", "Beta v", min = 0, value = 0.2, step = 0.05)),
        column(6, shiny::numericInput("strata.shape2", "Beta w", min = 0, value = 1, step = 0.05))
      ),
      shiny::h4("Genotype"),
      shiny::fluidRow(
        column(4, shiny::numericInput("pr.aa", "AA", min = 0, value = 0.00001, step = 0.00001)),
        column(4, shiny::numericInput("pr.ab", "Ab", min = 0, value = 0.00001, step = 0.00001)),
        column(4, shiny::numericInput("pr.bb", "bb", min = 0, value = 0.00001, step = 0.00001))
      ),
      shiny::sliderInput("locus.xlim", "Zoom", min = 0, max = 0.05, value = c(0, 0.05), dragRange = TRUE),
      shiny::actionButton("return", h5("Return parameters")),
      shiny::hr(),
      shiny::h4("Visualize joint distribution"),
      shiny::numericInput("num.loc", "Number of loci", min = 0, value = 10000, step = 1000),
      shiny::numericInput("num.ind", "Number of individuals", min = 0, value = 100, step = 10),
      shiny::actionButton("draw", "Draw from distributions"),
      shiny::hr(),
      shiny::numericInput("num.bins", "Number of bins", min = 1, value = 50)
    ),
    shiny::mainPanel(
      shiny::plotOutput("dist"),
      shiny::plotOutput("locus.hist"),
      shiny::tableOutput("smry")
    )
  )
)
