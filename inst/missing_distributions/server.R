server <- function(input, output) {
  output$dist <- renderPlot({
    tibble(
      x = seq(0, 1, length.out = 10000),
      Loci = dbeta(x, input$loc.shape1, input$loc.shape2),
      Individuals = dbeta(x, input$ind.shape1, input$ind.shape2),
      Strata = dbeta(x, input$strata.shape1, input$strata.shape2)
    ) %>%
      gather(type, density, -x) %>%
      mutate(type = factor(type, levels = c("Loci", "Individuals", "Strata"))) %>% 
      ggplot(aes(x, density)) +
      geom_line() +
      facet_grid(type ~ ., scales = "free_y") +
      xlim(input$locus.xlim) +
      labs(x = "Pr(missing)", y = "Density", title = "Marginal PDFs")
  })

  pr.missing <- eventReactive(input$draw, ({
    loc <- rbeta(input$num.loc, input$loc.shape1, input$loc.shape2)
    ind <- rbeta(input$num.ind, input$ind.shape1, input$ind.shape2)
    loc %o% ind
  }))

  output$locus.hist <- renderPlot({
    ggplot(data.frame(pr = as.vector(pr.missing())), aes(pr)) +
      geom_histogram(bins = input$num.bins) +
      xlim(input$locus.xlim) +
      labs(
        x = "Pr(missing)", 
        y = "Number of genotypes", 
        title = "Joint loci x individual distribution"
      )
  })

  output$smry <- renderTable({
    data.frame(
      N = length(pr.missing()),
      Mean = mean(pr.missing()),
      Median = median(pr.missing()),
      Max = max(pr.missing()),
      LCI = unname(quantile(pr.missing(), 0.025)),
      UCI = unname(quantile(pr.missing(), 0.975))
    )
  })
  
  observeEvent(input$return, {
    stopApp(list(
      loc = c(shape1 = input$loc.shape1, shape2 = input$loc.shape2),
      ind = c(shape1 = input$ind.shape1, shape2 = input$ind.shape2),
      strata = c(shape1 = input$strata.shape1, shape2 = input$strata.shape2),
      genotype = c(aa = input$pr.aa, ab = input$pr.ab, bb = input$pr.bb)
    ))
  })
}