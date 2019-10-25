context("getGini")
library(SCOPE)

test_that("Gini Calculation works", {
  Gini <- get_gini(Y_sim)
  expect_equal(sum(Gini<=0.12), 2)
})
