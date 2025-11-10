# Create datasets for machine learning

This vignette identifies and writes to disk datasets for use in
multi-dataset machine learning. For each relevant study condition, two
files are written:

1.  A `.rda` file containing a TreeSummarizedExperiment with full
    colData, rowData, a phylogenetic tree, and assay data containing
    taxonomic data.
2.  A `.csv.gz` file containing only key sample metadata, and taxonomic
    data (variables in columns)

See [source file on
GitHub](https://github.com/waldronlab/curatedMetagenomicAnalyses/blob/main/vignettes/MLdatasets.Rmd)

Do `BiocManager::install("waldronlab/curatedMetagenomicAnalyses")` to
install the `makeSEforCondition` function and to run this vignette.

Packages used here:

``` r

library(curatedMetagenomicData)
library(curatedMetagenomicAnalyses)
library(dplyr)
```

## Investigate potential response variables

These are the 10 study conditions most commonly found in
curatedMetagenomicData:

``` r

data("sampleMetadata")
availablediseases <- pull(sampleMetadata, study_condition) %>%
  table() %>%
  sort(decreasing = TRUE)
availablediseases
```

    ## .
    ##                   control                       IBD                       T2D 
    ##                     15121                      2088                       882 
    ##                       CRC                       IGT            premature_born 
    ##                       701                       563                       448 
    ##                  melanoma                       CAD                      ACVD 
    ##                       315                       269                       214 
    ##                   adenoma                       FMT                 cirrhosis 
    ##                       209                       149                       132 
    ##                       STH                    otitis             schizophrenia 
    ##                       108                       107                       106 
    ##              hypertension                        AS                        HF 
    ##                        99                        97                        95 
    ##                       T1D           acute_diarrhoea          pre-hypertension 
    ##                        89                        56                        56 
    ##                       CDI                    ME/CFS                  migraine 
    ##                        53                        50                        49 
    ##                      STEC               fatty_liver                 psoriasis 
    ##                        42                        41                        41 
    ## carcinoma_surgery_history                        AD            cephalosporins 
    ##                        40                        38                        36 
    ##                        PD                    asthma             periodontitis 
    ##                        31                        24                        24 
    ##                       SRP          peri-implantitis                        BD 
    ##                        24                        23                        20 
    ##                 mucositis                bronchitis   TKI_dependent_diarrhoea 
    ##                        20                        18                        16 
    ##            respiratoryinf        metabolic_syndrome                     donor 
    ##                        13                        10                         9 
    ##            pyelonephritis infectiousgastroenteritis                      MDRB 
    ##                         6                         5                         5 
    ##                     fever                 pneumonia               tonsillitis 
    ##                         3                         3                         3 
    ##                     cough             pyelonefritis                   skininf 
    ##                         2                         2                         2 
    ##                stomatitis                  cystitis             salmonellosis 
    ##                         2                         1                         1 
    ##                    sepsis                   suspinf 
    ##                         1                         1

And the number of studies they are found in:

``` r

studies <- lapply(names(availablediseases), function(x){
  filter(sampleMetadata, study_condition %in% x) %>%
    pull(study_name) %>%
    unique()
})
names(studies) <- names(availablediseases)
studies <- studies[-grep("control", names(studies))] #get rid of controls
studies <- studies[sapply(studies, length) > 1] #available in more than one study
studies
```

    ## $IBD
    ## [1] "HallAB_2017"     "HMP_2019_ibdmdb" "IaniroG_2022"    "IjazUZ_2017"    
    ## [5] "LiJ_2014"        "NielsenHB_2014"  "VilaAV_2018"    
    ## 
    ## $T2D
    ## [1] "HMP_2019_t2d"           "KarlssonFH_2013"        "LiJ_2014"              
    ## [4] "MetaCardis_2020_a"      "QinJ_2012"              "SankaranarayananK_2015"
    ## 
    ## $CRC
    ##  [1] "FengQ_2015"      "GuptaA_2019"     "HanniganGD_2017" "ThomasAM_2018a" 
    ##  [5] "ThomasAM_2018b"  "ThomasAM_2019_c" "VogtmannE_2016"  "WirbelJ_2018"   
    ##  [9] "YachidaS_2019"   "YuJ_2015"        "ZellerG_2014"   
    ## 
    ## $IGT
    ## [1] "HMP_2019_t2d"      "KarlssonFH_2013"   "MetaCardis_2020_a"
    ## 
    ## $premature_born
    ## [1] "BrooksB_2017" "OlmMR_2017"  
    ## 
    ## $melanoma
    ## [1] "FrankelAE_2017"       "GopalakrishnanV_2018" "LeeKA_2022"          
    ## [4] "MatsonV_2018"         "PetersBA_2019"        "WindTT_2020"         
    ## 
    ## $adenoma
    ## [1] "FengQ_2015"      "HanniganGD_2017" "ThomasAM_2018a"  "YachidaS_2019"  
    ## [5] "ZellerG_2014"   
    ## 
    ## $FMT
    ## [1] "IaniroG_2022" "LiSS_2016"   
    ## 
    ## $cirrhosis
    ## [1] "LoombaR_2017" "QinN_2014"   
    ## 
    ## $STH
    ## [1] "RosaBA_2018"  "RubelMA_2020"
    ## 
    ## $schizophrenia
    ## [1] "Castro-NallarE_2015" "ZhuF_2020"          
    ## 
    ## $T1D
    ## [1] "Heitz-BuschartA_2016" "KosticAD_2015"        "LiJ_2014"            
    ## 
    ## $acute_diarrhoea
    ## [1] "DavidLA_2015" "KieserS_2018"
    ## 
    ## $CDI
    ## [1] "IaniroG_2022"  "VincentC_2016"

Each of these datasets has six data types associated with it; for
example:

``` r

curatedMetagenomicData("JieZ_2017.+")
```

    ## 2021-03-31.JieZ_2017.gene_families
    ## 2021-03-31.JieZ_2017.marker_abundance
    ## 2021-03-31.JieZ_2017.marker_presence
    ## 2021-03-31.JieZ_2017.pathway_abundance
    ## 2021-03-31.JieZ_2017.pathway_coverage
    ## 2021-03-31.JieZ_2017.relative_abundance
    ## 2021-10-14.JieZ_2017.gene_families
    ## 2021-10-14.JieZ_2017.marker_abundance
    ## 2021-10-14.JieZ_2017.marker_presence
    ## 2021-10-14.JieZ_2017.pathway_abundance
    ## 2021-10-14.JieZ_2017.pathway_coverage
    ## 2021-10-14.JieZ_2017.relative_abundance

## Write relative abundance datasets to disk

``` r

for (i in seq_along(studies)){
  cond <- names(studies)[i]
  se <-
    curatedMetagenomicAnalyses::makeSEforCondition(cond, removestudies = "HMP_2019_ibdmdb", dataType = "relative_abundance")
  print(paste("Next study condition:", cond, " /// Body site: ", unique(colData(se)$body_site)))
  print(with(colData(se), table(study_name, study_condition)))
  cat("\n \n")
  save(se, file = paste0(cond, ".rda"))
  flattext <- select(as.data.frame(colData(se)), c("study_name", "study_condition", "subject_id"))
  rownames(flattext) <- colData(se)$sample_id
  flattext <- cbind(flattext, data.frame(t(assay(se))))
  write.csv(flattext, file = paste0(cond, ".csv"))
  system(paste0("gzip ", cond, ".csv"))
}
```

    ## [1] "Next study condition: IBD  /// Body site:  stool"
    ##                 study_condition
    ## study_name       control IBD
    ##   HallAB_2017         74 185
    ##   IaniroG_2022         3   3
    ##   IjazUZ_2017         38  56
    ##   LiJ_2014            10 140
    ##   NielsenHB_2014     248 148
    ##   VilaAV_2018          0 355
    ## 
    ## 

    ## [1] "Next study condition: T2D  /// Body site:  stool"
    ##                         study_condition
    ## study_name               control T2D
    ##   HMP_2019_t2d                46  11
    ##   KarlssonFH_2013             43  53
    ##   LiJ_2014                    10  79
    ##   MetaCardis_2020_a          642 550
    ##   QinJ_2012                  174 170
    ##   SankaranarayananK_2015      18  19
    ## 
    ## 

    ## [1] "Next study condition: CRC  /// Body site:  stool"
    ##                  study_condition
    ## study_name        control CRC
    ##   FengQ_2015           61  46
    ##   GuptaA_2019          30  30
    ##   HanniganGD_2017      28  27
    ##   ThomasAM_2018a       24  29
    ##   ThomasAM_2018b       28  32
    ##   ThomasAM_2019_c      40  40
    ##   VogtmannE_2016       52  52
    ##   WirbelJ_2018         65  60
    ##   YachidaS_2019       251 258
    ##   YuJ_2015             54  74
    ##   ZellerG_2014         61  53
    ## 
    ## 

    ## [1] "Next study condition: IGT  /// Body site:  stool"
    ##                    study_condition
    ## study_name          control IGT
    ##   HMP_2019_t2d           46 239
    ##   KarlssonFH_2013        43  49
    ##   MetaCardis_2020_a     642 275
    ## 
    ##  
    ## [1] "Next study condition: premature_born  /// Body site:  stool"     
    ## [2] "Next study condition: premature_born  /// Body site:  oralcavity"
    ## [3] "Next study condition: premature_born  /// Body site:  skin"      
    ##               study_condition
    ## study_name     control premature_born
    ##   BrooksB_2017       5            403
    ##   OlmMR_2017         0             45
    ## 
    ## 

    ## [1] "Next study condition: melanoma  /// Body site:  stool"
    ##                       study_condition
    ## study_name             melanoma
    ##   FrankelAE_2017             39
    ##   GopalakrishnanV_2018       25
    ##   LeeKA_2022                165
    ##   MatsonV_2018               39
    ##   PetersBA_2019              27
    ##   WindTT_2020                20
    ## 
    ## 

    ## [1] "Next study condition: adenoma  /// Body site:  stool"
    ##                  study_condition
    ## study_name        adenoma control
    ##   FengQ_2015           47      61
    ##   HanniganGD_2017      26      28
    ##   ThomasAM_2018a       27      24
    ##   YachidaS_2019        67     251
    ##   ZellerG_2014         42      61
    ## 
    ## 

    ## [1] "Next study condition: FMT  /// Body site:  stool"
    ##               study_condition
    ## study_name     control FMT
    ##   IaniroG_2022       3 109
    ##   LiSS_2016          5  40
    ## 
    ## 

    ## [1] "Next study condition: cirrhosis  /// Body site:  stool"
    ##               study_condition
    ## study_name     cirrhosis control
    ##   LoombaR_2017         9      36
    ##   QinN_2014          123     114
    ## 
    ## 

    ## [1] "Next study condition: STH  /// Body site:  stool"
    ##               study_condition
    ## study_name     control STH
    ##   RosaBA_2018        5  19
    ##   RubelMA_2020      86  89
    ## 
    ## 

    ## [1] "Next study condition: schizophrenia  /// Body site:  oralcavity"
    ## [2] "Next study condition: schizophrenia  /// Body site:  stool"     
    ##                      study_condition
    ## study_name            control schizophrenia
    ##   Castro-NallarE_2015      16            16
    ##   ZhuF_2020                81            90
    ## 
    ## 

    ## [1] "Next study condition: T1D  /// Body site:  stool"
    ##                       study_condition
    ## study_name             control T1D
    ##   Heitz-BuschartA_2016      26  27
    ##   KosticAD_2015             89  31
    ##   LiJ_2014                  10  31
    ## 
    ## 

    ## [1] "Next study condition: acute_diarrhoea  /// Body site:  stool"
    ##               study_condition
    ## study_name     acute_diarrhoea control
    ##   DavidLA_2015              38       9
    ##   KieserS_2018              18       9
    ## 
    ## 

    ## [1] "Next study condition: CDI  /// Body site:  stool"
    ##                study_condition
    ## study_name      CDI control
    ##   IaniroG_2022   20       3
    ##   VincentC_2016  33     196
    ## 
    ## 

## Direct link to files

Download the .csv and .rda files directly from
<https://www.dropbox.com/sh/0t0nbhj9eqm3wkq/AACZIw42WA-uHjzo97bG5tE6a?dl=0>
