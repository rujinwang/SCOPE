context("getGini")
library(SCOPE)

test_that("Gini Calculation works", {
  Gini = getGini(Y_sim)
  expect_equal(sum(Gini<=0.12), 2)
})
