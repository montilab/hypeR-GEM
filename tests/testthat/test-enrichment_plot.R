test_that("Enrichment returns expected output",{

  ## template signature and geneset
  data(COVID_urine)
  data(reactome)

  hypeR_GEM_obj <- hypeR.GEM::signature2gene(signatures = COVID_urine,
                                             directional = FALSE,
                                             species = "Human",
                                             merge = TRUE,
                                             promiscuous_threshold = 500,
                                             ensemble_id = TRUE,
                                             reference_key = 'refmet_name',
                                             background = NULL)

  enrichment_obj <- hypeR.GEM::enrichment(hypeR_GEM_obj,
                                  genesets = reactome,
                                  genesets_name = "REACTOME",
                                  method='unweighted',
                                  min_metabolite = 2,
                                  background=3068)


  result <- hypeR.GEM::enrichment_plot(enrichment_obj,
                                        top=40,
                                        abrv=50,
                                        size_by="significance",
                                        fdr_cutoff= 0.05,
                                        val='fdr')



  testthat::expect_named(result$data,
                         c("label","signature","significance", "geneset", "size"))










})
