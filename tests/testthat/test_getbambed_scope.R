context("getbambed")
library(SCOPE)
library(WGSmapp)

test_that("File preparation works", {
  bamfolder <- system.file("extdata", package = "WGSmapp")
  bamFile <- list.files(bamfolder, pattern = '*.dedup.bam$')
  bamdir <- file.path(bamfolder, bamFile)
  sampname_raw = sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 1)
  bambedObj <- getbambed_scope(bamdir = bamdir, sampname = sampname_raw)
  bamdir <- bambedObj$bamdir
  sampname_raw <- bambedObj$sampname
  ref_raw <- bambedObj$ref

  expect_equal(length(bamFile), length(sampname_raw))
  expect_gte(length(as.character(unique(seqnames(ref_raw)))), 22)
})
