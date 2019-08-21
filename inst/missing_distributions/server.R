server <- function(input, output) {
  output$dist <- shiny::renderPlot({
    tibble::tibble(
      x = seq(0, 1, length.out = 10000),
      Loci = stats::dbeta(x, input$loc.shape1, input$loc.shape2),
      Individuals = stats::dbeta(x, input$ind.shape1, input$ind.shape2),
      Strata = stats::dbeta(x, input$strata.shape1, input$strata.shape2)
    ) %>%
      tidyr::gather(type, density, -x) %>%
      dplyr::mutate(type = factor(type, levels = c("Loci", "Individuals", "Strata"))) %>% 
      ggplot2::ggplot(ggplot2::aes(x, density)) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(type ~ ., scales = "free_y") +
      ggplot2::xlim(input$locus.xlim) +
      ggplot2::labs(x = "Pr(missing)", y = "Density", title = "Marginal PDFs")
  })

  pr.missing <- shiny::eventReactive(input$draw, ({
    loc <- stats::rbeta(input$num.loc, input$loc.shape1, input$loc.shape2)
    ind <- stats::rbeta(input$num.ind, input$ind.shape1, input$ind.shape2)
    loc %o% ind
  }))

  output$locus.hist <- shiny::renderPlot({
    ggplot2::ggplot(data.frame(pr = as.vector(pr.missing())), ggplot2::aes(pr)) +
      ggplot2::geom_histogram(bins = input$num.bins) +
      ggplot2::xlim(input$locus.xlim) +
      ggplot2::labs(
        x = "Pr(missing)", 
        y = "Number of genotypes", 
        title = "Joint loci x individual distribution"
      )
  })

  output$smry <- shiny::renderTable({
    data.frame(
      N = length(pr.missing()),
      Mean = mean(pr.missing()),
      Median = stats::median(pr.missing()),
      Max = max(pr.missing()),
      LCI = unname(stats::quantile(pr.missing(), 0.025)),
      UCI = unname(stats::quantile(pr.missing(), 0.975))
    )
  })
  
  shiny::observeEvent(input$return, {
    shiny::stopApp(list(
      loc = c(shape1 = input$loc.shape1, shape2 = input$loc.shape2),
      ind = c(shape1 = input$ind.shape1, shape2 = input$ind.shape2),
      strata = c(shape1 = input$strata.shape1, shape2 = input$strata.shape2),
      genotype = c(aa = input$pr.aa, ab = input$pr.ab, bb = input$pr.bb)
    ))
  })
}
