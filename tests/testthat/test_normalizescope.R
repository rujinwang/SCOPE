context("normalizescope")
library(SCOPE)

test_that("basic argument errors thrown", {
  Gini = getGini(Y_sim)
  normObj.sim <- normalize_codex2_ns_noK(Y_qc =Y_sim,
                                         gc_qc = ref_sim$gc,
                                         norm_index = which(Gini<=0.12))

  Yhat.noK.sim=normObj.sim$Yhat
  beta.hat.noK.sim=normObj.sim$beta.hat
  ploidy.sim =  PreEst_ploidy(Y = Y_sim, Yhat = Yhat.noK.sim, ref = ref_sim)
  expect_error(normalize_scope_foreach(Y_qc = Y_sim, gc_qc = ref_sim$gc,
                                       K = 1:3, ploidyInt = ploidy.sim,
                                       norm_index = which(Gini<=0.12), T = 1:7,
                                       beta0 = beta.hat.noK.sim), "exceed the number")
})
