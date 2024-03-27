test_that("signature2gene returns expected output",{

  ## template signature
  data(COVID_urine)

  result <- hypeR.GEM::signature2gene(signatures = COVID_urine,
                                      directional = FALSE,
                                      species = "Human",
                                      merge = TRUE,
                                      promiscuous_threshold = 500,
                                      ensemble_id = TRUE,
                                      reference_key = 'refmet_name',
                                      background = NULL)

  testthat::expect_length(result, 3L)
  testthat::expect_type(result, "list")
  testthat::expect_named(result,
                         c("promiscuous_met","mapped_metabolite_signatures","gene_tables"))


  testthat::expect_type(result$mapped_metabolite_signatures,"list")
  testthat::expect_length(result$mapped_metabolite_signatures, 2L)
  testthat::expect_type(result$gene_tables,"list")
  testthat::expect_length(result$gene_tables, 2L)

  testthat::expect_named(result$mapped_metabolite_signatures$up,
                         c("name","fullname","refmet_name","gene_association",
                           "formula", "exact_mass", "super_class", "main_class",
                           "sub_class"))

  testthat::expect_named(result$gene_tables$up,
                         c("name", "symbol", "associated_reactions",
                           "total_association", "signature_association", "signature_size" ,
                           "background", "p_value", "fdr", "one_minus_p_value",
                           "one_minus_fdr", "associated_metabolites"))
})
