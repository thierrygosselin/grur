ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      h4("Select Pr(missing)"),
      h4("Loci"),
      fluidRow(
        column(6, numericInput("loc.shape1", "Beta v", min = 0, value = 0.2, step = 0.05)),
        column(6, numericInput("loc.shape2", "Beta w", min = 0, value = 1, step = 0.05))
      ),
      h4("Individuals"),
      fluidRow(
        column(6, numericInput("ind.shape1", "Beta v", min = 0, value = 0.2, step = 0.05)),
        column(6, numericInput("ind.shape2", "Beta w", min = 0, value = 1, step = 0.05))
      ),
      h4("Strata"),
      fluidRow(
        column(6, numericInput("strata.shape1", "Beta v", min = 0, value = 0.2, step = 0.05)),
        column(6, numericInput("strata.shape2", "Beta w", min = 0, value = 1, step = 0.05))
      ),
      h4("Genotype"),
      fluidRow(
        column(4, numericInput("pr.aa", "AA", min = 0, value = 0.00001, step = 0.00001)),
        column(4, numericInput("pr.ab", "Ab", min = 0, value = 0.00001, step = 0.00001)),
        column(4, numericInput("pr.bb", "bb", min = 0, value = 0.00001, step = 0.00001))
      ),
      sliderInput("locus.xlim", "Zoom", min = 0, max = 0.05, value = c(0, 0.05), dragRange = TRUE),
      actionButton("return", h5("Return parameters")),
      hr(),
      h4("Visualize joint distribution"),
      numericInput("num.loc", "Number of loci", min = 0, value = 10000, step = 1000),
      numericInput("num.ind", "Number of individuals", min = 0, value = 100, step = 10),
      actionButton("draw", "Draw from distributions"),
      hr(),
      numericInput("num.bins", "Number of bins", min = 1, value = 50)
    ),
    mainPanel(
      plotOutput("dist"),
      plotOutput("locus.hist"),
      tableOutput("smry")
    )
  )
)
