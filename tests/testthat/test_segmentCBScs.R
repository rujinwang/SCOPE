context("CrossSampleSegmentation")
library(SCOPE)

test_that("Cross-sample CBS segmentation works", {
  chrs = unique(as.character(seqnames(ref_sim)))
  segment_cs = vector('list',length = length(chrs))
  names(segment_cs) = chrs
  Yhat.sim = normObj.scopeDemo$Yhat[[which.max(normObj.scopeDemo$BIC)]]
  for (chri in chrs) {
    message('\n', chri, '\n')
    segment_cs[[chri]] = segmentCBScs(Y = Y_sim,
                                      Yhat = Yhat.sim,
                                      sampname = colnames(Y_sim),
                                      ref = ref_sim,
                                      chr = chri,
                                      mode = "integer", max.ns = 1)
  }
  iCN = do.call(rbind, lapply(segment_cs, function(z){z[["iCN"]]}))

  expect_equal(nrow(Y_sim), nrow(iCN))
  expect_equal(ncol(Y_sim), ncol(iCN))
})
