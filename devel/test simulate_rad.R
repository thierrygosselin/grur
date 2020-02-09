rm(list = ls())
library(grur)

sim <- grur::simulate_rad(
  num.pops = 3,
  num.loci = 200,
  div.time = 1000,
  ne = c(50, 200),
  nm = c(0, 0.1),
  theta = 0.2,
  mig.type = c("island", "stepping.stone"),
  num.reps = 5,
  label = "sim.rad.test",
  parallel.core = 3
 )
