test_that("Enrichment returns expected output",{

  ## template signature and geneset
  data(COVID_urine)
  data(Reactome)

  hypeR_GEM_obj <- hypeR.GEM::signature2gene(signatures = COVID_urine,
                                             directional = FALSE,
                                             species = "Human",
                                             merge = TRUE,
                                             promiscuous_threshold = 500,
                                             ensemble_id = TRUE,
                                             reference_key = 'refmet_name',
                                             background = NULL)

  result <- hypeR.GEM::enrichment(hypeR_GEM_obj,
                                  genesets = Reactome$genesets,
                                  genesets_name = "REACTOME",
                                  method='unweighted',
                                  min_metabolite = 2,
                                  background=3068)

  testthat::expect_length(result, 2L)
  testthat::expect_type(result, "list")

  testthat::expect_length(result$up, 2L)
  testthat::expect_type(result$up, "list")

  testthat::expect_named(result$up,
                         c("info","data"))

  testthat::expect_named(result$up$info,
                         c("Test","Signature_size", "Genesets", "Background"))

  testthat::expect_named(result$up$data,
                         c("label","pval", "fdr", "signature",
                           "geneset", "overlap", "weighted_overlap",
                           "background", "gene_hits", "metabolite_hits",
                           "num_met_hits", "ratio_met_hits" ))



})
