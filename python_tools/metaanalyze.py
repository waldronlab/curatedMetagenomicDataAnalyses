#!/usr/bin/env python

import argparse as ap
import sys, os
import pandas as pd
import numpy as np
sys.path.append(os.path.join(os.path.dirname(sys.path[0]), "python_modules"))
from meta_analysis_runner import singleStudyEffect
from meta_analysis_runner import meta_analysis_with_linear_model
from meta_analysis_runner import analysis_with_linear_model
from meta_analyses import paule_mandel_tau
from meta_analyses import RE_meta_binary
from meta_analyses import RE_meta


def read_params():
    p = ap.ArgumentParser(description=\
	"This script runs a meta-analysis. "
	"It asks for: "
	"\n A metadata + data table (with samples as columns and fields as index)"
	"\n A flag for the analysis, mandatory: can be -re , -mc or -sre (see after)"
	"\n A feature identifier to detect which are the features and which are metadata "
	"\n A formula for the model, of the form : age + BMI + psoriasis "
	"\n\t the script will perform the meta-analysis on the first variable "
	"\n it will include an intercept automatically. NOTE that this script work "
	"\n only with a OLS model (no mixed, no interaction). "
	"\n Example (sex-contrast): "
	"\n python metaanalyze.py \\"
	"\n\tsex_contrast_meta_analysis.tsv \\"
	"\n\t-z s__ -re -cc female:male \\"
	"\n\t--formula \"C(gender) + age + BMI\"", \
	formatter_class=ap.RawTextHelpFormatter)
    add = p.add_argument
    add("metadata_data_table", type=str, help="Samples in columns, fields in index: "
	"must include both metadata you want ot study and features. ")
    add("-z", "--feat_id", type=str, default="s__", help="The feature identifier will "
	"\n refer to all the features. Remember that the default is s__ (species level rel.abundances)")
    add("-re", "--re", action="store_true", help="Activate a random effect (categorical) meta-analysis")
    add("-mc", "--mc", action="store_true", help="Activate a meta-correlation analsysis: continouus case")
    add("-sre", "--sre", action="store_true", help="Is like the -re, but works on a single dataset.")
    add("-cc", "--controlcase", type=str, default=None, help="Use ':' : -cc female:male will set control and cases for an analysis on sex, e.g..")
    add("-fm", "--formula", type=str, default="", help="Var of study + covariates (cancer + age + BMI) "
	"\n **REFERRING TO COVARIATES** use patsy formulas, e.g. C(sex) for defining a catgorical variable. "
	"\n NOTE that: for -re and -sre is always necessary to specify which is the negative and the positive via --controlcase")
    add("-H", "--heterogeneity", type=str, choices=["FIX", "PM", "DL"], help="Heterogeneity (not active for -sre. Default is PM (Paule-Mandel) "
	"NOTE THAT by setting '-H FIX' you use a fixed-effect model", default="PM")
    add("-si", "--studyid", type=str, default="study_name", help="Name of column 'Study'. Default=[study_name]")
    add("-of", "--outfile", type=str, default="", help="If not set, attach _metaanalysis to the name of the input.")
    return p.parse_args()


def main(args):

    datatable = pd.read_csv(args.metadata_data_table, sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")

    if (not args.outfile): args.outfile = args.metadata_data_table.replace(".tsv", "_metaanalysis.tsv")
    feats = [j for j in datatable.index.tolist() if (args.feat_id in j)]
    metadata = [k for k in datatable.index.tolist() if (not k in feats)]
    
    if np.count_nonzero(np.array([args.re, args.mc, args.sre], dtype=np.int64))!=1:
        raise SyntaxError("Sorry: you have to specify ONE in -re, -mc or -sre")

    if args.re:
        if (not args.controlcase):
            raise SyntaxError("For --re you *must* include a -cc var with control and cases defined")
        neg, pos = tuple(args.controlcase.split(":"))
        ma = meta_analysis_with_linear_model(datatable, args.formula, args.studyid, feats, args.outfile, "CLS", args.heterogeneity, pos=pos, neg=neg)
        ma.random_effect_regression_model()

    elif args.mc:
        ma = meta_analysis_with_linear_model(datatable, args.formula, args.studyid, feats, args.outfile, "REG", args.heterogeneity, pos=None, neg=None)
        ma.random_effect_regression_model()
        
    elif args.sre:
        if (not args.controlcase):
            raise SyntaxError("For --sre you *must* include a -cc var with control and cases defined")
        neg, pos = tuple(args.controlcase.split(":"))
        alm = analysis_with_linear_model(datatable, args.formula, feats, args.outfile, pos, neg)
        alm.write_out()


if __name__ == "__main__":
    main(read_params())
