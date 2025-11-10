#!/bin/bash

mkdir -p images


python ../../python_tools/metaanalyze.py clr_abundances/sex_species.tsv \
    -of meta_analyses/clr_sex_species.tsv \
    -z s__ -re -cc "female:male" --formula "gender + age + BMI + number_reads" -si study_identifier

python ../../python_tools/metaanalyze.py clr_abundances/age_species.tsv \
    -of meta_analyses/clr_age_species.tsv \
    -z s__ -mc --formula "age + BMI + C(gender) + number_reads" -si study_identifier

python ../../python_tools/metaanalyze.py clr_abundances/bmi_species.tsv \
    -of meta_analyses/clr_bmi_species.tsv \
    -z s__ -mc --formula "BMI + age + C(gender) + number_reads" -si study_identifier


### 

exit

###

python ../../python_tools/metaanalyze.py arcsin_abundances/sex_species.tsv -of meta_analyses/arcsin_sex_species.tsv \
    -z s__ -re -cc "female:male" --formula "gender + age + BMI + number_reads" -si study_identifier

python ../../python_tools/metaanalyze.py arcsin_abundances/age_species.tsv -of meta_analyses/arcsin_age_species.tsv \
    -z s__ -mc --formula "age + BMI + C(gender) + number_reads" -si study_identifier

python ../../python_tools/metaanalyze.py arcsin_abundances/bmi_species.tsv -of meta_analyses/arcsin_bmi_species.tsv \
    -z s__ -mc --formula "BMI + age + C(gender) + number_reads" -si study_identifier

python scatter_clr_vs_asin.py

#python ../curatedMetagenomicDataAnalyses/python_tools/metaanalyze.py ../clr_abundances/bmi_genera.tsv -of ../meta_analyses/bmi_genera.tsv \
#    -z g__ -mc --formula "BMI + age + C(gender) + number_reads" -si study_identifier

# *** ALPHA DIVERSITY
#python ../curatedMetagenomicDataAnalyses/python_tools/metaanalyze.py ../relative_abundances/sex_alpha_diversity.tsv -of ../meta_analyses/sex_alpha_diversity.tsv \
#    -z ad__ -re --formula "gender + BMI + age + number_reads" -si study_identifier -cc female:male
#python ../curatedMetagenomicDataAnalyses/python_tools/metaanalyze.py ../relative_abundances/age_alpha_diversity.tsv -of ../meta_analyses/age_alpha_diversity.tsv \
#    -z ad__ -mc --formula "age + BMI + C(gender) + number_reads" -si study_identifier
#python ../curatedMetagenomicDataAnalyses/python_tools/metaanalyze.py ../relative_abundances/bmi_alpha_diversity.tsv -of ../meta_analyses/bmi_alpha_diversity.tsv \
#    -z ad__ -mc --formula "BMI + age + C(gender) + number_reads" -si study_identifier

# *** BETA DIVERSITY
#python ../curatedMetagenomicDataAnalyses/python_tools/hierarchical_binary_test.py ../clr_abundances/sex_species.tsv -a permanova 
#python ../curatedMetagenomicDataAnalyses/python_tools/hierarchical_binary_test.py ../clr_abundances/sex_species.tsv -a anosim

## SPECIES FIGURES
#python ../curatedMetagenomicDataAnalyses/python_tools/draw_figure_with_ma.py ../meta_analyses/sex_species.tsv ../meta_analyses/sex_genera.tsv --names SPC GN \
#    --outfile images/SPC_FIG2_meta --how joint_top -ps "Sex: male" -ns "Sex: female" --x_axis "Standardized mean difference" --y_axis "Microbial taxa" \
#    --title "Meta-analysis of sex-related microbial taxa" --a_single 0.2 --a_random 0.01 -re "RE_Effect" -ci "RE_conf_int" -rq "RE_Effect_Qvalue" \
#    -es "_Effect" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.20 --neg_max_rho 0.8 --pos_max_rho 0.8 \
#    -a ../relative_abundances/sex_species.tsv ../relative_abundances/sex_genera.tsv --narrowed --imp 30 --check_genera -ms 5

#python ../curatedMetagenomicDataAnalyses/python_tools/draw_figure_with_ma.py ../meta_analyses/age_species.tsv ../meta_analyses/age_genera.tsv --names SPC GN \
#    --outfile images/SPC_FIG3_meta_age --how joint_top -ps "Older age" -ns "Younger age" --x_axis "Partial correlation with age" --y_axis "Microbial taxa" \
#    --title "Meta-analysis of age-associated microbial taxa" --a_single 0.2 --a_random 0.01 -re "RE_Correlation" -ci "RE_conf_int" -rq "RE_Correlation_Qvalue" \
#    -es "_Correlation" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.10 --neg_max_rho 0.4 --pos_max_rho 0.6 \
#    -a ../relative_abundances/age_species.tsv ../relative_abundances/age_genera.tsv --narrowed --imp 30 --check_genera -ms 5

#python ../curatedMetagenomicDataAnalyses/python_tools/draw_figure_with_ma.py ../meta_analyses/bmi_species.tsv ../meta_analyses/bmi_genera.tsv --names SPC GN \
#    --outfile images/SPC_FIG3_meta_bmi --how joint_top -ps "Higher BMI" -ns "Lower BMI" --x_axis "Partial correlation with BMI" \
#    --y_axis "Microbial taxa" --title "Meta-analysis of BMI-associated microbial taxa" --a_single 0.2 --a_random 0.01 -re "RE_Correlation" -ci "RE_conf_int" \
#    -rq "RE_Correlation_Qvalue" -es "_Correlation" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.10 --neg_max_rho 0.4 --pos_max_rho 0.4 \
#    -a ../relative_abundances/bmi_species.tsv ../relative_abundances/bmi_genera.tsv --narrowed --imp 30 --check_genera -ms 5 #--size_on_err

#exit

# *** PATHWAYS
#python ../curatedMetagenomicDataAnalyses/python_tools/metaanalyze.py ../clr_abundances/sex_pathways.tsv -of ../meta_analyses/sex_pathways.tsv \
#    -z PWY -re --formula "gender + age + BMI + number_reads" -si study_identifier -cc "female:male"
#python ../curatedMetagenomicDataAnalyses/python_tools/metaanalyze.py ../clr_abundances/age_pathways.tsv -of ../meta_analyses/age_pathways.tsv \
#    -z PWY -mc --formula "age + C(gender) + BMI + number_reads" -si study_identifier
#python ../curatedMetagenomicDataAnalyses/python_tools/metaanalyze.py ../clr_abundances/bmi_pathways.tsv -of ../meta_analyses/bmi_pathways.tsv \
#    -z PWY -mc --formula "BMI + C(gender) + age + number_reads" -si study_identifier

# FIGURES
#python ../curatedMetagenomicDataAnalyses/python_tools/draw_figure_with_ma.py ../meta_analyses/sex_pathways.tsv --names PWY \
#    --outfile images/SEX_PATHWAYS --how first -ps "Sex: male" -ns "Sex: female" --x_axis "Standardized mean difference" --y_axis "Metabolic pathway" \
#    --title "Meta-analysis of sex-related microbiome metabolic pathway" --a_single 0.2  -re "RE_Effect" -ci "RE_conf_int" -rq "RE_Effect_Qvalue" \
#    -es "_Effect" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.10 --neg_max_rho 0.4 --pos_max_rho 0.4 \
#    -a ../relative_abundances/sex_pathways.tsv --narrowed --imp 30 -ms 5 --a_random 0.01

#python ../curatedMetagenomicDataAnalyses/python_tools/draw_figure_with_ma.py ../meta_analyses/age_pathways.tsv --names PWY \
#    --outfile images/AGE_PATHWAYS --how first -ps "Older" -ns "Younger" --x_axis "Partial correlation coefficient" --y_axis "Metabolic pathway" \
#    --title "Meta-analysis of age-related microbiome metabolic pathways" --a_single 0.2 -re "RE_Correlation" -ci "RE_conf_int" -rq "RE_Correlation_Qvalue" \
#    -es "_Correlation" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.10 --neg_max_rho 0.4 --pos_max_rho 0.4 \
#    -a ../relative_abundances/age_pathways.tsv --narrowed --imp 30 -ms 5 --a_random 0.01

#python ../curatedMetagenomicDataAnalyses/python_tools/draw_figure_with_ma.py ../meta_analyses/bmi_pathways.tsv --names PWY \
#    --outfile images/BMI_PATHWAYS --how first -ps "Higher BMI" -ns "Lower BMI" --x_axis "Partial correlation coefficient" --y_axis "Metabolic pathway" \
#    --title "Meta-analysis of BMI-related microbiome metabolic pathway" --a_single 0.2 -re "RE_Correlation" -ci "RE_conf_int" -rq "RE_Correlation_Qvalue" \
#    -es "_Correlation" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.10 --neg_max_rho 0.4 --pos_max_rho 0.4 \
#    -a ../relative_abundances/bmi_pathways.tsv --narrowed --imp 30 -ms 5 --a_random 0.01

# *** KEGG
#python ../curatedMetagenomicDataAnalyses/python_tools/metaanalyze.py ../clr_abundances/sex_kegg.tsv -of ../meta_analyses/sex_kegg.tsv \
#    -z K -re --formula "gender + age + BMI + number_reads" -si study_identifier -cc "female:male"
#python ../curatedMetagenomicDataAnalyses/python_tools/metaanalyze.py ../clr_abundances/age_kegg.tsv -of ../meta_analyses/age_kegg.tsv \
#    -z K -mc --formula "age + C(gender) + BMI + number_reads" -si study_identifier
#python ../curatedMetagenomicDataAnalyses/python_tools/metaanalyze.py ../clr_abundances/bmi_kegg.tsv -of ../meta_analyses/bmi_kegg.tsv \
#    -z K -mc --formula "BMI + C(gender) + age + number_reads" -si study_identifier 

#python ../curatedMetagenomicDataAnalyses/python_tools/draw_figure_with_ma.py ../meta_analyses/sex_kegg.tsv --names KEGG \
#    --outfile images/SEX_KEGG --how first -ps "Sex: male" -ns "Sex: female" --x_axis "Standardized mean difference" --y_axis "KEGG Ortholog" \
#    --title "Meta-analysis of sex-related microbiome KEGG Orthologs" --a_single 0.2  -re "RE_Effect" -ci "RE_conf_int" -rq "RE_Effect_Qvalue" \
#    -es "_Effect" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.10 --neg_max_rho 0.4 --pos_max_rho 0.4 \
#    -a ../relative_abundances/sex_kegg.tsv --narrowed --imp 30 -ms 5 --a_random 0.01

#python ../curatedMetagenomicDataAnalyses/python_tools/draw_figure_with_ma.py ../meta_analyses/age_kegg.tsv --names KEGG \
#    --outfile images/AGE_KEGG --how first -ps "Older" -ns "Younger" --x_axis "Partial correlation coefficient" --y_axis "KEGG Ortholog" \
#    --title "Meta-analysis of age-related microbiome KEGG Orthologs" --a_single 0.2 -re "RE_Correlation" -ci "RE_conf_int" -rq "RE_Correlation_Qvalue" \
#    -es "_Correlation" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.10 --neg_max_rho 0.4 --pos_max_rho 0.4 \
#    -a ../relative_abundances/age_kegg.tsv --narrowed --imp 30 -ms 5 --a_random 0.01

#python ../curatedMetagenomicDataAnalyses/python_tools/draw_figure_with_ma.py ../meta_analyses/bmi_kegg.tsv --names KEGG \
#    --outfile images/BMI_KEGG --how first -ps "Lower BMI" -ns "Higher BMI" --x_axis "Partial correlation coefficient" --y_axis "KEGG Ortholog" \
#    --title "Meta-analysis of BMI-related microbiome KEGG Orthologs" --a_single 0.2 -re "RE_Correlation" -ci "RE_conf_int" -rq "RE_Correlation_Qvalue" \
#    -es "_Correlation" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.10 --neg_max_rho 0.4 --pos_max_rho 0.4 \
#    -a ../relative_abundances/bmi_kegg.tsv --narrowed --imp 30 -ms 5 --a_random 0.01

## figures for revision:
#python ../curatedMetagenomicDataAnalyses/python_tools/metaanalyze.py ../clr_abundances/sex_species.tsv -of ../meta_analyses/sex_species_fixed_effect.tsv \
#    -z s__ -re -cc "female:male" --formula "gender + age + BMI + number_reads" -si study_identifier -H FIX

#python ../curatedMetagenomicDataAnalyses/python_tools/metaanalyze.py ../clr_abundances/sex_genera.tsv -of ../meta_analyses/sex_genera_fixed_effect.tsv  \
    #-z g__ -re -cc "female:male" --formula "gender + age + BMI + number_reads" -si study_identifier -H FIX

#python ../curatedMetagenomicDataAnalyses/python_tools/draw_figure_with_ma.py ../meta_analyses/sex_species_fixed_effect.tsv ../meta_analyses/sex_genera_fixed_effect.tsv --names SPC GN \
#    --outfile images/sex_scaled_with_fixed_effect --how joint_top -ps "Sex: male" -ns "Sex: female" --x_axis "Standardized mean difference" --y_axis "Microbial taxa" \
#    --title "Meta-analysis of sex-related microbial taxa" --a_single 0.2 --a_random 0.01 -re "RE_Effect" -ci "RE_conf_int" -rq "RE_Effect_Qvalue" \
#    -es "_Effect" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.20 --neg_max_rho 0.3 --pos_max_rho 0.3 \
#    -a ../relative_abundances/sex_species.tsv ../relative_abundances/sex_genera.tsv --narrowed --imp 30 --check_genera -ms 5 --size_on_err
 
#python ../curatedMetagenomicDataAnalyses/python_tools/draw_figure_with_ma.py ../meta_analyses/sex_species.tsv ../meta_analyses/sex_genera.tsv --names SPC GN \
#    --outfile images/sex_scaled_with_paule_mandel --how joint_top -ps "Sex: male" -ns "Sex: female" --x_axis "Standardized mean difference" --y_axis "Microbial taxa" \
#    --title "Meta-analysis of sex-related microbial taxa" --a_single 0.2 --a_random 0.01 -re "RE_Effect" -ci "RE_conf_int" -rq "RE_Effect_Qvalue" \
#    -es "_Effect" -qs "_Qvalue" -cbl black -cr darkgoldenrod -cb dodgerblue -il 0.20 --neg_max_rho 0.8 --pos_max_rho 0.8 \
#    -a ../relative_abundances/sex_species.tsv ../relative_abundances/sex_genera.tsv --narrowed --imp 30 --check_genera -ms 5 --size_on_err

