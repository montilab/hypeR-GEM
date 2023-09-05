
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hypeR.GEM

## (1) Installation

-   Using `devtools` package

``` r
library(devtools)
devtools::install_github("montilab/hypeR-GEM")
```

## Usage

``` r
library(hypeR.GEM)
```

## Load signatures

``` r
data(EL_signature)

str(EL_signature)
```

    #> List of 2
    #>  $ EL_up  :'data.frame': 311 obs. of  5 variables:
    #>   ..$ CHEMICAL_NAME: chr [1:311] "c-glycosyltryptophan" "vanillactate" "hydroxyasparagine**" "cysteine s-sulfate" ...
    #>   ..$ refmet_name  : chr [1:311] "alpha-C-Mannosyltryptophan" "Vanillactic acid" "Hydroxyasparagine" "Cysteine-S-sulfate" ...
    #>   ..$ cent.Estimate: num [1:311] 0.824 1.049 0.858 0.882 1.44 ...
    #>   ..$ cent.p.value : num [1:311] 5.25e-14 3.92e-13 2.69e-12 3.77e-12 6.46e-12 ...
    #>   ..$ cent.q.value : num [1:311] 5.52e-11 2.06e-10 9.42e-10 9.92e-10 1.36e-09 ...
    #>  $ EL_down:'data.frame': 96 obs. of  5 variables:
    #>   ..$ CHEMICAL_NAME: chr [1:96] "epiandrosterone sulfate" "dehydroepiandrosterone sulfate (dhea-s)" "androstenediol (3beta,17beta) monosulfate (1)" "x-11378" ...
    #>   ..$ refmet_name  : chr [1:96] "Epiandrosterone" "Dehydroepiandrosterone sulfate" "4-Androsten-3beta,17beta-diol 17-sulfate" NA ...
    #>   ..$ cent.Estimate: num [1:96] -1.9 -1.44 -1.64 -1.12 -1.75 ...
    #>   ..$ cent.p.value : num [1:96] 7.73e-08 1.35e-07 1.95e-07 6.28e-07 1.07e-06 ...
    #>   ..$ cent.q.value : num [1:96] 1.77e-06 2.91e-06 3.51e-06 9.31e-06 1.40e-05 ...

## Map metabolite signtures to associated genes

-   `signatures` must be a named list
-   `species = c("Human", "Mouse", "Rat", "Zebrafish", "Worm", "Other")`
-   `merge`: merge metabolites from different department
-   `promiscuous_threshold`: gene association threshold of promiscuous
    metabolite
-   `ensemble_id`: for current version, if `species != 'Human`, use
    `ensemble_id = FALSE`
-   `reference_key = 'refmet_name` by default
-   `background` is used to compute gene-specific p-values, if
    `background = NULL`, then background = \# of metabolites associated
    with non-exchange reactions

``` r
hypeR_GEM_obj <- hypeR.GEM::signature2gene(signatures = EL_signature,
                                           species = "Human",
                                           merge = TRUE,
                                           promiscuous_threshold = 50,
                                           ensemble_id = TRUE,
                                           reference_key = 'refmet_name',
                                           background = NULL)


hypeR_GEM_obj$promiscuous_met
```

    #>        refmet_name gene_association
    #> 1    Phenylalanine               57
    #> 2          Proline               57
    #> 3 Arachidonic acid               77
    #> 4    Linoleic acid               62

``` r
library(hypeR)
kegg <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory = "CP:KEGG", clean = TRUE)
print(kegg)
```

    #> C2.CP:KEGG v7.5.1 
    #> Abc Transporters (44)
    #> Acute Myeloid Leukemia (57)
    #> Adherens Junction (73)
    #> Adipocytokine Signaling Pathway (67)
    #> Alanine Aspartate And Glutamate Metabolism (32)
    #> Aldosterone Regulated Sodium Reabsorption (42)

## Hypergeometric test

### KEGG (unweighted)

``` r
max_fdr <- 0.01
hypeR_GEM_enrichments_unweighted <- hypeR.GEM::enrichment(hypeR_GEM_obj,
                                         genesets = kegg$genesets,
                                         genesets_name = "KEGG",
                                         method='unweighted',
                                         background=3068)

hypeR_GEM_enrichments_unweighted$EL_down$info
```

    #> $Test
    #> [1] "Hypergeometric test"
    #> 
    #> $Signature_size
    #> [1] 182
    #> 
    #> $Genesets
    #> [1] "KEGG"
    #> 
    #> $Background
    #> [1] 3068

#### Plot

``` r
hypeR.GEM::enrichment_plot(hypeR_GEM_enrichments_unweighted,
                           top=40,
                           abrv=50,
                           size_by="significance",
                           fdr_cutoff= max_fdr,
                           val='fdr')+
  ggplot2::ggtitle(paste("FDR ≤", max_fdr)) +
  ggplot2::theme(axis.text.y = element_text(size = 7))
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

#### Reactable

Not quite sure why the reactable is not displayed normally in the .md
file, but the function should work

``` r
#hypeR.GEM::rctbls(hypeR_GEM_enrichments_unweighted, fdr_cutoff=0.05)
```

### KEGG (weighted)

``` r
max_fdr <- 0.01
hypeR_GEM_enrichments_weighted <- hypeR.GEM::enrichment(hypeR_GEM_obj,
                                         genesets = kegg$genesets,
                                         genesets_name = "KEGG",
                                         method='weighted',
                                         weights = 'one_minus_fdr',
                                         background=3068)
```

#### Plot

``` r
hypeR.GEM::enrichment_plot(hypeR_GEM_enrichments_weighted,
                           top=40,
                           abrv=50,
                           size_by="significance",
                           fdr_cutoff= max_fdr,
                           val='fdr')+
  ggplot2::ggtitle(paste("FDR ≤", max_fdr)) +
  ggplot2::theme(axis.text.y = element_text(size = 7))
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

#### Reactable

``` r
#hypeR.GEM::rctbls(hypeR_GEM_enrichments_weighted, fdr_cutoff=0.05)
```
