# curatedMetagenomicDataAnalyses main paper
 
This folder contains final command lines used to produce analyses and figures of the paper "Meta-analysis of 20,533 human metagenomes defines an index of oral to gut microbial introgression and associations with age, sex, BMI, and diseases".

## ./all_command_lines.sh

This bash wrapper contains all the analysis used to generate the paper's figures and text.
It is divided into 4 main sections, which correspond to the four main conceptual blocks that need adequate resources and time to be completed.
The command lines reported include the code generaring the figures.
The four main blocks are:

1. SECTION 1: SPECIES, GENUS, PWY, and KOs META-ANALYSIS ON SEX, AGE, and BMI:
2. SECTION 2: HIERARCHICAL META-ANALYSIS ON 12 DISEASES
3. SECTION 3: MACHINE LEARNING ON DISEASES
4. SECTION 4: ORAL INTROGRESSION ANALYSIS AND META-ANALYSIS

### Requirements
 
You need 
* [metAML](https://github.com/SegataLab/metaml/).
* [pandas](https://pandas.pydata.org/) >= 1.0.3
* [scikit-learn](https://scikit-learn.org/stable/) >= 0.24.2
* [statsmodels](https://www.statsmodels.org/stable/index.html) >= 1.11.1
* [scikit-bio](http://scikit-bio.org/) >= 0.5.6
* [numpy](https://numpy.org/) >= 1.18.1
* [scipy](https://scipy.org/) >= 1.4.1
* [matplotlib](https://matplotlib.org/) >= 3.3.4

### Getting Started

Most of the meta-analysis used can be launched with the following script:

    `python ../python_tools/metaanalyze.py -h`

* This program is a caller for several types of microbiome-dedicated meta-analysese including several that are not used in the aforementioned paper.
* Please feel free to contact paolo.manghi@unitn.it for curiosities and questions about these meta-analyses!

```
usage: metaanalyze.py [-h] [-z FEAT_ID] [-uma] [-usma] [-re] [-mc] [-sre]
                      [-sor] [-OR] [-scor] [-COR] [-cc CONTROLCASE]
                      [-fm FORMULA] [-H {FIX,PM,DL}] [-si STUDYID]
                      [-of OUTFILE]
                      metadata_data_table

This script runs a meta-analysis. It asks for: 
 A metadata + data table (with samples as columns and fields as index)
 A flag for the analysis, mandatory: can be -re , -mc or -sre (see after)
 A feature identifier to detect which are the features and which are metadata 
 A formula for the model, of the form : age + BMI + psoriasis 
	 the script will perform the meta-analysis on the first variable 
 it will include an intercept automatically. NOTE that this script work 
 only with a OLS model (no mixed, no interaction). 
 Example (sex-contrast): 
 python metaanalyze.py \
	sex_contrast_meta_analysis.tsv \
	-z s__ -re -cc female:male \
	--formula "C(gender) + age + BMI"

positional arguments:
  metadata_data_table   Samples in columns, fields in index: must include both metadata you want ot study and features. 

optional arguments:
  -h, --help            show this help message and exit
  -z FEAT_ID, --feat_id FEAT_ID
                        The feature identifier will 
                         refer to all the features. Remember that the default is s__ (species level rel.abundances)
  -uma, --uma           Activate a random effect (categorical) meta-analysis: real interpretation of the responses
  -usma, --usma         Is like the -uma, but works on a single dataset: real interpreation of the respnses
  -re, --re             Activate a random effect (categorical) meta-analysis
  -mc, --mc             Activate a meta-correlation analsysis: continouus case
  -sre, --sre           Is like the -re, but works on a single dataset.
  -sor, --sor           Activate logistic regression for single dataset (using odd ratio)
  -OR, --OR             Activate meta-analysis using logistic regression (using odd ratio)
  -scor, --scor         Activate logistic regression for a continuous predictor for single dataset
  -COR, --COR           Activate logistic regression for a continuous predictor for meta-analysis 
  -cc CONTROLCASE, --controlcase CONTROLCASE
                        Use ':' : -cc female:male will set control and cases for an analysis on sex, e.g..
  -fm FORMULA, --formula FORMULA
                        Var of study + covariates (cancer + age + BMI) 
                         **REFERRING TO COVARIATES** use patsy formulas, e.g. C(sex) for defining a catgorical variable. 
                         NOTE that: for -re and -sre is always necessary to specify which is the negative and the positive via --controlcase
  -H {FIX,PM,DL}, --heterogeneity {FIX,PM,DL}
                        Heterogeneity (not active for -sre. Default is PM (Paule-Mandel) NOTE THAT by setting '-H FIX' you use a fixed-effect model
  -si STUDYID, --studyid STUDYID
                        Name of column 'Study'. Default=[study_name]
  -of OUTFILE, --outfile OUTFILE
                        If not set, attach _metaanalysis to the name of the input.
```

