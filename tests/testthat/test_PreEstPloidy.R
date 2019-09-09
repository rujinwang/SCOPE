context("PreEst_Ploidy")
library(SCOPE)

test_that("Ploidy Initialization works", {
  Gini = getGini(Y_sim)
  normObj.sim <- normalize_codex2_ns_noK(Y_qc =Y_sim,
                                         gc_qc = ref_sim$gc,
                                         norm_index = which(Gini<=0.12))

  Yhat.noK.sim=normObj.sim$Yhat
  ploidy.sim =  PreEst_ploidy(Y = Y_sim, Yhat = Yhat.noK.sim, ref = ref_sim)
  expect_equal(sum(Gini<=0.12), sum(ploidy.sim==2))
})
