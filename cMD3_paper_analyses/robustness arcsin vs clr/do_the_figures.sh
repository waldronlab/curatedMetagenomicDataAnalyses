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

