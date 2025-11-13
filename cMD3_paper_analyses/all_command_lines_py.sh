#!/bin/bash

mkdir -p images
mkdir -p ML/ml_dis_rf/results

## SECTION 1: SPECIES, GENUS, PWY, and KOs META-ANALYSIS ON SEX, AGE, and BMI
python ../python_tools/metaanalyze.py clr_tables/sex-species_clr.tsv -z s__ -re -cc "0.0:1.0" --formula "gender + age + BMI + number_reads" 
python ../python_tools/metaanalyze.py clr_tables/sex-genera_clr.tsv -z g__ -re -cc "0.0:1.0" --formula "gender + age + BMI + number_reads"  
python ../python_tools/metaanalyze.py clr_tables/sex-pathways_clr.tsv -z PWY -re -cc "0.0:1.0" --formula "gender + age + BMI + number_reads"   
python ../python_tools/metaanalyze.py clr_tables/sex-kegg_clr.tsv -z K -re -cc "0.0:1.0" --formula "gender + age + BMI + number_reads" -si dataset_name

python ../python_tools/metaanalyze.py clr_tables/age-species_clr.tsv -z s__ -mc --formula "age + BMI + C(gender) + number_reads"
python ../python_tools/metaanalyze.py clr_tables/age-genera_clr.tsv -z g__ -mc --formula "age + BMI + C(gender) + number_reads"
python ../python_tools/metaanalyze.py clr_tables/age-pathways_clr.tsv -z PWY -mc --formula "age + BMI + C(gender) + number_reads"
python ../python_tools/metaanalyze.py clr_tables/age-kegg_clr_substitution.tsv -z K -mc --formula "age + BMI + C(gender) + number_reads" #-si dataset_name

python ../python_tools/metaanalyze.py clr_tables/bmi-species_clr.tsv -z s__ -mc --formula "BMI + age + C(gender) + number_reads"
python ../python_tools/metaanalyze.py clr_tables/bmi-genera_clr.tsv -z g__ -mc --formula "BMI + age + C(gender) + number_reads"
python ../python_tools/metaanalyze.py clr_tables/bmi-pathways_clr.tsv -z PWY -mc --formula "BMI + age + C(gender) + number_reads"
python ../python_tools/metaanalyze.py clr_tables/bmi-kegg_clr_substitution.tsv -z K -mc --formula "BMI + age + C(gender) + number_reads" #-si dataset_name
 
python ../python_tools/draw_figure_with_ma.py clr_tables/sex-species_clr_metaanalysis.tsv clr_tables/sex-genera_clr_metaanalysis.tsv --names SPC GN --outfile images/SPC_FIG2_meta_5precent_prevalence --how joint_top     -ps "Sex: male" -ns "Sex: female" --x_axis "Standardized mean difference" --y_axis "Microbial taxa" --title "Meta-analysis of sex-related microbial taxa" --a_single 0.2 --a_random 0.05 -re "RE_Effect"     -ci "RE_conf_int" -rq "RE_Effect_Qvalue" -es "_Effect" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.20 --neg_max_rho 0.8 --pos_max_rho 0.8  --narrowed --boxes --imp 30 --narrowed --boxes --imp 30 --check_genera -ms 2

python ../python_tools/draw_figure_with_ma.py clr_tables/age-pathways_clr_metaanalysis.tsv --outfile images/PWY_FIG3_meta_supp1 --how first -ps "Older age associated" -ns "Younger age associated" --x_axis "Partial correlation with age (yrs.)" --y_axis "Metabolic pathway" --title "Meta-analysis of age-associated metabolic pathways" --a_single 0.2 --a_random 0.05 -re "RE_Correlation" -ci "RE_conf_int" -rq "RE_Correlation_Qvalue" -es "_Correlation" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.10 --neg_max_rho 0.4 --pos_max_rho 0.4 --narrowed --boxes --imp 30 -ms 2
    
python ../python_tools/draw_figure_with_ma.py clr_tables/bmi-pathways_clr_metaanalysis.tsv --outfile images/PWY_FIG3_meta_supp2 --how first -ps "Higher BMI associated" -ns "Lower BMI associated" --x_axis "Partial correlation with BMI (kg/m^2)" --y_axis "Metabolic pathway" --title "Meta-analysis of BMI-associated metabolic pathways" --a_single 0.2 --a_random 0.05 -re "RE_Correlation" -ci "RE_conf_int" -rq "RE_Correlation_Qvalue" -es "_Correlation" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.10 --neg_max_rho 0.4 --pos_max_rho 0.4 --narrowed --boxes --imp 30 -ms 2
  
python ../python_tools/draw_figure_with_ma.py clr_tables/sex-kegg_clr_metaanalysis.tsv --outfile images/KOS_FIG2_meta_supp2 --how first -ps "Sex: male" -ns "Sex: female" --x_axis "Standardized mean difference" --y_axis "" --title "Meta-analysis of sex-related KOs" --a_single 0.2 --a_random 0.05 -re "RE_Effect" -ci "RE_conf_int" -rq "RE_Effect_Qvalue" -es "_Effect" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.20 --neg_max_rho 0.6 --pos_max_rho 0.6 --narrowed --boxes --imp 30 -ms 2
  
python ../python_tools/draw_figure_with_ma.py clr_tables/age-kegg_clr_substitution_metaanalysis.tsv --outfile images/KOS_FIG3_meta_supp3 --how first -ps "Older age associated" -ns "Younger age associated" --x_axis "Partial correlation with age (yrs.)" --y_axis "" --title "Meta-analysis of age-associated KOs" --a_single 0.2 --a_random 0.05 -re "RE_Correlation" -ci "RE_conf_int" -rq "RE_Correlation_Qvalue" -es "_Correlation" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.10 --neg_max_rho 0.4 --pos_max_rho 0.4 -a --narrowed --boxes --imp 30 -ms 2
     
python ../python_tools/draw_figure_with_ma.py clr_tables/bmi-kegg_clr_substitution_metaanalysis.tsv --outfile images/KOS_FIG3_meta_supp4 --how first -ps "Higher BMI associated" -ns "Lower BMI associated" --x_axis "Partial correlation with BMI (kg/m^2)" --y_axis "" --title "Meta-analysis of BMI-associated KOs" --a_single 0.2 --a_random 0.05 -re "RE_Correlation" -ci "RE_conf_int" -rq "RE_Correlation_Qvalue" -es "_Correlation" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.10 --neg_max_rho 0.4 --pos_max_rho 0.4 --narrowed --boxes --imp 30 -ms 2

## SECTION 1.1 - SENSITIVITY ANALYSIS CLR vs asin(sqrt(relative abundance))

cd robustness\ arcsin\ vs\ clr/
bash do_the_figures.sh


## SECTION 2: HIERARCHICAL META-ANALYSIS ON 12 DISEASES
mkdir -p clr_metaaans_nest
bzip2 -d nested_clr_tables/*.tsv.bz2

python ../python_tools/metaanalyze.py nested_clr_tables/BD_species_clr.tsv          -z s__ --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:BD -sre -H PM
python ../python_tools/metaanalyze.py nested_clr_tables/CRC_pathways_clr.tsv       -z PWY --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:positive -re -H PM

python ../python_tools/metaanalyze.py nested_clr_tables/T2D_species_clr.tsv  -z s__ --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:T2D -re -H PM
python ../python_tools/metaanalyze.py nested_clr_tables/ACVD_pathways_clr.tsv    -z PWY --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:positive -sre -H PM

python ../python_tools/metaanalyze.py nested_clr_tables/CRC_species_clr.tsv  -z s__      --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:CRC -re -H PM
python ../python_tools/metaanalyze.py nested_clr_tables/schizofrenia_pathways_clr.tsv  -z PWY --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:positive -sre -H PM

python ../python_tools/metaanalyze.py nested_clr_tables/ACVD_species_clr.tsv     -z s__ --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:ACVD -sre -H PM
python ../python_tools/metaanalyze.py nested_clr_tables/CD_pathways_clr.tsv       -z PWY  --formula "disease_subtype + C(gender) + age + BMI + number_reads" -cc control:positive -re -H PM

python ../python_tools/metaanalyze.py nested_clr_tables/schizofrenia_species_clr.tsv   -z s__ --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:schizofrenia -sre -H PM
python ../python_tools/metaanalyze.py nested_clr_tables/UC_pathways_clr.tsv -z PWY --formula "disease_subtype + C(gender) + age + BMI + number_reads" -cc control:positive -re -H PM

python ../python_tools/metaanalyze.py nested_clr_tables/CD_species_clr.tsv          -z s__ --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:IBD -re -H PM
python ../python_tools/metaanalyze.py nested_clr_tables/MECFS_pathways_clr.tsv     -z PWY --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:positive -sre -H PM

python ../python_tools/metaanalyze.py nested_clr_tables/UC_species_clr.tsv -z s__ --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:IBD -re -H PM
python ../python_tools/metaanalyze.py nested_clr_tables/asthma_pathways_clr.tsv  -z PWY --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:positive -sre -H PM

python ../python_tools/metaanalyze.py nested_clr_tables/MECFS_species_clr.tsv    -z s__   --formula "study_condition + C(gender) + age + BMI + number_reads" -cc "control:ME/CFS" -sre -H PM
python ../python_tools/metaanalyze.py nested_clr_tables/STH_pathways_clr.tsv -z PWY --formula "study_condition + C(gender) + age + BMI + number_reads" -cc neg:positive -sre -H PM

python ../python_tools/metaanalyze.py nested_clr_tables/asthma_species_clr.tsv   -z s__ --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:asthma -sre -H PM
python ../python_tools/metaanalyze.py nested_clr_tables/cirrhosis_pathways_clr.tsv  -z PWY --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:positive -sre -H PM

python ../python_tools/metaanalyze.py nested_clr_tables/STH_species_clr.tsv  -z s__ --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:STH -sre -H PM
python ../python_tools/metaanalyze.py nested_clr_tables/cirrhosis_species_clr.tsv -z s__  --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:cirrhosis -sre -H PM

python ../python_tools/metaanalyze.py nested_clr_tables/migraine_pathways_clr.tsv  -z PWY --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:positive -sre -H PM
python ../python_tools/metaanalyze.py nested_clr_tables/BD_pathways_clr.tsv      -z PWY --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:positive -sre -H PM

python ../python_tools/metaanalyze.py nested_clr_tables/migraine_species_clr.tsv   -z s__ --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:migraine -sre -H PM
python ../python_tools/metaanalyze.py nested_clr_tables/T2D_pathways_clr.tsv -z PWY --formula "study_condition + C(gender) + age + BMI + number_reads" -cc control:positive -re -H PM

python ../python_tools/hierarchical_metaanalysis.py \
    --names CRC T2D CD UC ACVD cirrhosis schizofrenia asthma STH BD "ME/CFS" migraine \
    -k re -H PM -se RE_stdErr -es RE_Effect --outfile clr_metaaans_nest/SPC_cMD3_paper_two_layers.tsv \
    -qs RE_Effect_Qvalue RE_Effect_Qvalue RE_Effect_Qvalue RE_Effect_Qvalue \
    Qvalue Qvalue Qvalue Qvalue Qvalue Qvalue Qvalue Qvalue \
    -ps RE_Pvalue RE_Pvalue RE_Pvalue RE_Pvalue \
    RE_Pvalue RE_Pvalue RE_Pvalue RE_Pvalue RE_Pvalue RE_Pvalue RE_Pvalue RE_Pvalue \
    --sheets \
    nested_clr_tables/CRC_species_clr_metaanalysis.tsv nested_clr_tables/T2D_species_clr_metaanalysis.tsv nested_clr_tables/CD_species_clr_metaanalysis.tsv \
    nested_clr_tables/UC_species_clr_metaanalysis.tsv nested_clr_tables/ACVD_species_clr_metaanalysis.tsv nested_clr_tables/cirrhosis_species_clr_metaanalysis.tsv \
    nested_clr_tables/schizofrenia_species_clr_metaanalysis.tsv nested_clr_tables/asthma_species_clr_metaanalysis.tsv nested_clr_tables/STH_species_clr_metaanalysis.tsv \
    nested_clr_tables/BD_species_clr_metaanalysis.tsv nested_clr_tables/MECFS_species_clr_metaanalysis.tsv nested_clr_tables/migraine_species_clr_metaanalysis.tsv

python ../python_tools/hierarchical_metaanalysis.py \
    --names CRC T2D CD UC ACVD cirrhosis schizofrenia asthma STH BD "ME/CFS" migraine \
    -k re -H PM -se RE_stdErr -es RE_Effect --outfile clr_metaaans_nest/PWY_cMD3_paper_two_layers.tsv \
    -qs RE_Effect_Qvalue RE_Effect_Qvalue RE_Effect_Qvalue RE_Effect_Qvalue \
    Qvalue Qvalue Qvalue Qvalue Qvalue Qvalue Qvalue Qvalue \
    -ps RE_Pvalue RE_Pvalue RE_Pvalue RE_Pvalue \
    RE_Pvalue RE_Pvalue RE_Pvalue RE_Pvalue RE_Pvalue RE_Pvalue RE_Pvalue RE_Pvalue \
    --sheets \
    nested_clr_tables/CRC_pathways_clr_metaanalysis.tsv nested_clr_tables/T2D_pathways_clr_metaanalysis.tsv nested_clr_tables/CD_pathways_clr_metaanalysis.tsv \
    nested_clr_tables/UC_pathways_clr_metaanalysis.tsv nested_clr_tables/ACVD_pathways_clr_metaanalysis.tsv nested_clr_tables/cirrhosis_pathways_clr_metaanalysis.tsv \
    nested_clr_tables/schizofrenia_pathways_clr_metaanalysis.tsv nested_clr_tables/asthma_pathways_clr_metaanalysis.tsv nested_clr_tables/STH_pathways_clr_metaanalysis.tsv \
    nested_clr_tables/BD_pathways_clr_metaanalysis.tsv nested_clr_tables/MECFS_pathways_clr_metaanalysis.tsv nested_clr_tables/migraine_pathways_clr_metaanalysis.tsv

python ../python_tools/draw_figure_with_ma.py clr_metaaans_nest/SPC_cMD3_paper_two_layers.tsv --outfile images/SPC_FIG4_meta_10precent_prevalence -a relative_abundances/usable_all_species.tsv -re RE_Effect -x "Standardized mean difference" -es "_Effect" -qs "_Qvalue" --imp 30 --neg_max_rho 1.25 --pos_max_rho 1.25 --title "Meta-analysis of disease-associated microbial species" -ps "Unhealthy microbiome" -ns "Healthy microbiome" --legloc "lower right" -ar 0.05 -as 0.2 --color_red goldenrod --color_blue cornflowerblue --color_black black -rq RE_Qvalue --confint RE_conf_int --markers -ms 2 --boxes --legloc "upper left" 
 
python ../python_tools/draw_figure_with_ma.py clr_metaaans_nest/PWY_cMD3_paper_two_layers.tsv --outfile images/PWY_FIG4_meta_supp1 -a relative_abundances/usable_all_pathways.tsv -re RE_Effect -x "Standardized mean difference" -es "_Effect" -qs "_Qvalue" --imp 30 --neg_max_rho 1.0 --pos_max_rho 1.0 --title "Meta-analysis of disease-associated metabolic pathways" -ps "Unhealthy microbiome" -ns "Healthy microbiome" --legloc "lower right" -ar 0.05 -as 0.2 --color_red goldenrod --color_blue cornflowerblue --color_black black -rq RE_Qvalue --confint RE_conf_int --markers -ms 2 --boxes --legloc "upper left" #-mp 0.01 -mna 5

## SECTION 3: MACHINE LEARNING ON DISEASES
cd ML
python ml_tests_on_diseases_rf.py
python ../figure4_complete_ml.py 
cd ..

## SECTION 4: ORAL INTROGRESSION ANALYSIS AND META-ANALYSIS
cd oral_enrichment

python 01-estimate_oral_enrichment.py
python 02-oral_enrichment.py Oral_Richness -osp 0.01

python 03-superframe_to_stats.py figures_data/panelC_data.tsv figures_data/panelC_stats.tsv
python 03-superframe_to_stats.py figures_data/panelD_data.tsv figures_data/panelD_stats.tsv
