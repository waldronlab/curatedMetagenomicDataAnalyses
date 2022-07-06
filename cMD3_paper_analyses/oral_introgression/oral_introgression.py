#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys, os
import argparse as ap
import glob
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
import seaborn as sns
from functools import partial
import itertools as it
import subprocess as sb
from matplotlib import rc
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["font.size"] = 14
from scipy import stats as sts
from scipy.stats import rankdata
from scipy.spatial import distance
from scipy.cluster import hierarchy
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pingouin as pg
from sklearn.metrics import pairwise_distances
from sklearn.metrics import roc_auc_score
from metareg import RE_meta_binary, RE_meta, singleStudyEffect, paule_mandel_tau
from skbio.stats import composition as cps
from statsmodels.stats.diagnostic import lilliefors

## python oral_introgression.py Oral_Richness ## DISEASES
## python oral_introgression.py Oral_Richness -la age -tm AGE -dt tables/usable_for_age_species.tsv
## python oral_introgression.py Oral_Richness -la BMI -tm REG -dt tables/usable_for_BMI_species.tsv
## python oral_introgression.py Oral_Richness -la gender:male:female -tm SEX -dt tables/usable_for_sex_species.tsv

class oral_introgression(object):
    def readargs(self, args):
        p = ap.ArgumentParser()
        add = p.add_argument
        add("Oral_Richness", type=str)
        add("-dt", "--dataset", type=str, default="tables/a_usable_for_diseases_on_spp_complete.tsv")

        add("-st", "--stratify", type=str, default="study_name")
        add("-la", "--lookat", type=str, default="condition:positive:negative")
        ####add("-i", "--index", type=str, default="Oral_Ri")
        add("-os", "--oral_species", type=str, default="oral_data/Species_oral_cmd3.tsv" )
        add("-nm", "--name_of_list", type=str, default="list_of_species/Species_oral_cmd3_list.tsv")

        add("-z", "--featid", type=str, default="s__")
        add("-ma", "--min_ab", type=float, default=0.1)
        add("-mp", "--min_prev", type=float, default=5)
        add("-tm", "--type_of_meta", type=str, default="CLS", choices=["CLS", "REG", "AGE", "SEX"])
        return p.parse_args()

    def __init__(self, args):
        self.args = self.readargs(args)

        dts = [\
        "XieH_2016", "ZhuF_2020", "JieZ_2017", "QinN_2014", "RubelMA_2020", "YeZ_2018", "NagySzakalD_2017", "HMP_2019_ibdmdb", \
        "NielsenHB_2014", "QinJ_2012", "KarlssonFH_2013", "SankaranarayananK_2015", "ZellerG_2014", "YuJ_2015", "FengQ_2015", \
        "VogtmannE_2016", "ThomasAM_2019_a", "ThomasAM_2019_b", "WirbelJ_2018", "GuptaA_2019", "HanniganGD_2017", "YachidaS_2019"\
        ]
        conds = ["asthma", "migraine", "schizofrenia", "ACVD", "cirrhosis", "STH", "BD", "ME/CFS", "IBD", "T2D", "CRC", "control"]
        
        self.featid = self.args.featid

        if self.args.type_of_meta in ["CLS", "SEX"]:
            self.condition, self.positive, self.negative = tuple(self.args.lookat.split(":"))

        elif self.args.type_of_meta in ["REG", "AGE"]:
            self.condition = self.args.lookat

        ## Here U recreated the usable (3.162)
        self.input = pd.read_csv(self.args.dataset, sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")
        feats = [i for i in self.input.index if (i.startswith(self.args.featid))]

        # for s in self.input.columns:
        #    self.input.loc[ feats, s ] = 
 
        #if self.args.type_of_meta == "CLS":
        #    with open("other_files/samples_for_figure_5.tsv") as sf5:
        #        self.input = self.input[[sample.rstrip() for sample in sf5.readlines()]]
        
        if (not self.condition in self.input.index.tolist()) and (self.args.type_of_meta=="CLS"):
            self.input.loc[self.condition] = [(self.positive if c!="control" else self.negative) for c in self.input.loc["study_condition"].tolist()]

        self.oral_cMD = self.define_oral_bases(self.args.oral_species) 
        self.oral_estimate = self.quantify_oral_introgression()
        self.mapping_with_diseases = self.map_2_diseases()
       


    def define_oral_bases(self, oral_cmd):
        #if os.path.isfile(oral_cmd):
        return pd.read_csv(oral_cmd, sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")
        #else:
        #    OP = pd.read_csv("July_2021_usable_new.tsv", sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")
        #    OP = OP.loc[:, OP.loc["body_site"]=="oralcavity"]
        #    feats = [i for i in OP.index.tolist() if ((self.args.featid in i) and (np.sum(OP.loc[i].values.astype(float))>0.))]
        #    OP.index.name = "sample_id"
        #    OP.to_csv(oral_cmd, sep="\t", header=True, index=True)
        #    return OP
       

    def quantify_oral_introgression(self):
        print("Initializing oral-score estimation...")
        feats = [i for i in self.oral_cMD.index.tolist() if (self.featid in i)]

        print("Before filter you have: %i" %len(feats))

        abuns = dict([(ft, np.mean(self.oral_cMD.loc[ft].values.astype(float))) for ft in feats])
        prevs = dict([(ft, (np.count_nonzero(self.oral_cMD.loc[ft].values.astype(float))/float(self.oral_cMD.shape[1]))*100.) \
            for ft in feats])

        oral_taxa = [ft for ft in feats if ((abuns[ft]>self.args.min_ab) and (prevs[ft]>self.args.min_prev))]

        print("After filter: %i" %len(oral_taxa))

        with open(self.args.name_of_list, "w") as Os:
            for sp in oral_taxa: Os.write(sp + "\n")

        oral_sample = self.oral_cMD.loc[oral_taxa, :]        
   
        def entropy(proportions_array):
            P = np.array(proportions_array, dtype=np.float64)
            log_of_P = [ (np.log10(n) if ( np.isfinite( np.log10(n))  ) else 0.0) for n in proportions_array ]
            return -1.*(np.sum(P * np.array(log_of_P, dtype=np.float64)))
 
        def ginismp(proportions_array):
            return np.sum(proportions_array**2.)

        def oral_fraction(tot_proportions_array, or_proportions_array):
            return np.count_nonzero(or_proportions_array) / float(np.count_nonzero(tot_proportions_array))

        metadata = [i for i in self.input.index.tolist() if (not self.featid in i)]
        OSP = [osp for osp in oral_taxa if osp in (self.input.index.tolist())]
        oral_estimate = self.input.copy()

        oral_richness_corrected = dict([ (samplename, np.sum(oral_estimate.loc[OSP, samplename].values.astype(float))) \
	    for samplename in oral_estimate.columns.tolist()  ])

        oral_richness = dict([(samplename, np.sum(oral_estimate.loc[OSP, samplename].values.astype(float))) \
            for samplename in oral_estimate.columns.tolist()])

        oral_entropy = dict([(samplename, entropy(oral_estimate.loc[OSP, samplename].values.astype(float)/100.)) \
            for samplename in oral_estimate.columns.tolist()])
        oral_ginisimp = dict([(samplename, ginismp(oral_estimate.loc[OSP, samplename].values.astype(float)/100.)) \
            for samplename in oral_estimate.columns.tolist()])

        #oral_fraction_ = dict([(samplename, oral_fraction(\
	#    oral_estimate.loc[[i for i in oral_estimate.index.tolist() if (self.featid in i)], samplename].values.astype(float)/100., \
	#    oral_estimate.loc[OSP, samplename].values.astype(float)/100. \
	#    ))  for samplename in oral_estimate.columns.tolist()])

        Get = lambda dic,key : dic[key] if (not np.isnan(dic[key])) else 0.0 

        oral_estimate.loc["Oral_Richness"] = [Get(oral_richness, samplename) for samplename in oral_estimate.columns.tolist()]
        oral_estimate.loc["Oral_Richness_Corrected"] = [Get(oral_richness_corrected, samplename) for samplename in oral_estimate.columns.tolist()]
        oral_estimate.loc["Oral_Entropy"] = [Get(oral_entropy, samplename) for samplename in oral_estimate.columns.tolist()]
        oral_estimate.loc["Oral_Entropy_Corrected"] = oral_estimate.loc["Oral_Entropy"].values.astype(float)
        oral_estimate.loc["Oral_Gini"] = [Get(oral_ginisimp, samplename) for samplename in oral_estimate.columns.tolist()]
        #oral_estimate.loc["Oral_Fraction"] = [Get(oral_fraction_, samplename) for samplename in oral_estimate.columns.tolist()]
        print("Finished oral-score estimation.")

 

        #self.input
        table = oral_estimate.loc[["study_name", "study_condition", "Oral_Richness", "Oral_Entropy"]].T
        table.sort_values(["study_name", "study_condition"], inplace=True)
        table.index.name = "sample_id"

        #pd.DataFrame({"sample_id": oral_estimate.columns.tolist(), "study_name": oral_estimate.loc["study_name"].tolist(), \
        #    "study_condition": oral_estimate.loc["study_condition"].tolist(), "oral_richness": })

        #for sam in oral_estimate.columns.tolist():
        #    print(sam, )

        return oral_estimate
     


    def segregate_datasets(self):
        print("Initiating dataset segregation...")
        single_datasets = []
        for dataset in self.oral_estimate.loc["study_name"].unique():
            if (not dataset in ["XieH_2016", "HMP_2019_ibdmdb", "NielsenHB_2014"]):
                this_estimate = self.oral_estimate.loc[:, self.oral_estimate.loc["study_name"]==dataset]
                key = dataset
                single_datasets += [(key, this_estimate)]
            else:
                if dataset in ["HMP_2019_ibdmdb", "NielsenHB_2014"]:
                    for dis,anti_dis in [("UC", "CD"), ("CD", "UC")]:
                        this_estimate = self.oral_estimate.loc[:, self.oral_estimate.loc["study_name"]==dataset]
                        this_estimate = this_estimate.loc[:, this_estimate.loc["disease_subtype"]!=anti_dis]
                        key = dataset + "_" + dis
                        single_datasets += [(key, this_estimate)]
                elif dataset == "XieH_2016":
                    for dis,anti_dis in [("asthma", "migraine"), ("migraine", "asthma")]:
                        this_estimate = self.oral_estimate.loc[:, self.oral_estimate.loc["study_name"]==dataset]
                        this_estimate = this_estimate.loc[:, this_estimate.loc["study_condition"]!=anti_dis]
                        key = dataset + "_" + dis
                        single_datasets += [(key, this_estimate)]
        print("Terminated dataset aggregation")
        return dict(single_datasets)


    def regression(self, index, data, problem, NO_covars):
        datat = data.T
        datat.fillna("NA", inplace=True)
        datat = datat.loc[datat["age"] != "NA"]
        datat = datat.loc[datat["gender"] != "NA"]
        datat = datat.loc[datat["BMI"] != "NA"]
        datat[index] = datat[index].values.astype(float)
        #datat[ index ] = datat[ index ].values.astype(float)  

        datat[ index ] = np.arcsin(np.sqrt( datat[ index ].values.astype(float) /100.  ))

        Lens = len(datat)
        covariates = {\
            "gender": ["age", "BMI"], \
            "age": ["BMI", "gender"], \
            "BMI": ["age", "gender"], \
            "condition": ["BMI", "age", "gender"]}
        datat["age"] = datat["age"].values.astype(float)
        datat["BMI"] = datat["BMI"].values.astype(float)
        #datat["gender"] = [(1.0 if(p=="male") else 0.0) for p in datat["gender"].tolist()]

        #print(("condition" in datat.columns), " MA CHE CAZZO")
        #exit(1) 

        if ("condition" in datat.columns) and (not "condition" in covariates[problem]):
            covariates[problem] += ["condition"]

        datat = datat[[index, problem] + covariates[problem]]

        if NO_covars:
            formula = ("%s ~ %s" %(index, problem))
        else:
            formula = ("%s ~ " %index) + " + ".join([problem] + covariates[problem])

        md = smf.ols(formula, data=datat)
        model_fit = md.fit()
        t = model_fit.tvalues.loc[problem]
        n = float(len(datat))
        r = float(t) / np.sqrt(np.float((t**2.) + (n - 1.)))
        Zr = np.arctanh(r) #0.5 * np.log((1. + r) / float/(1. - r))
        SEr = 1/np.sqrt(n - 3)
 
        r_lw = Zr- (sts.t.ppf(0.975, Lens-1) * SEr)
        r_up = Zr+ (sts.t.ppf(0.975, Lens-1) * SEr)
        #r_lw = Zr - (1.96*SEr)
        #r_up = Zr + (1.96*SEr)
        return np.tanh(r), model_fit.pvalues.loc[problem], n, 0.0, "N"



    def simple_mean_diff_rev(self, index, data, problem, positive_class, negative_class, raw, no_covariates):
        datat = data.T
        datat.fillna("NA", inplace=True)
        datat[index] = datat[index].values.astype(float)
        datat["age"] = datat["age"].values.astype(float)
        datat["BMI"] = datat["BMI"].values.astype(float)

        covariates = {\
            "gender": ["age", "BMI"], \
            "age": ["BMI", "gender", "condition"], \
            "BMI": ["age", "gender"], \
            "condition": ["BMI", "age", "gender"]}
             
        #datat = datat[[index, problem] + covariates[problem]]
        #i = datat[ index ].values.astype(float)
        #i[ i == 0 ] = 0.00001
        #i = np.log2( i )
        #datat[ index ] = i

        #if not no_covariates:
        formula_adjusted = ("%s ~ " %problem) + " + ".join([index] + covariates[problem][:-1])
        #else:
        #formula_simple = ("%s ~ %s" %(problem, index))
        results = {}

        for fml,c in zip([formula_adjusted], ["adjusted"]): #["adjusted", "crude"]):
            md = smf.ols( fml, data=datat )
            model_fit = md.fit()

            #print(model_fit.summary())
            results[c] = {}
            #if not prob in model_fit.pvalues.index:
            #    results[c]["wald"] = 1.0
            #    results[c]["beta"] = np.nan
            #    results[c]["var"] = np.nan
            #    results[c]["ci"] = (np.nan, np.nan)
            #else:

            results[c]["wald"] = model_fit.pvalues.loc[index]
            results[c]["beta"] = model_fit.params.loc[index]
            results[c]["var"] = model_fit.bse.loc[index]**2.
            results[c]["ci"] = (results[c]["beta"]- (sts.t.ppf(0.975, len(datat)-1) * model_fit.bse.loc[index]), \
                results[c]["beta"]+ (sts.t.ppf(0.975, len(datat)-1) * model_fit.bse.loc[index]))
 
            #results[c]["wald"] = model_fit.pvalues.loc[problem]
            #results[c]["beta"] = model_fit.params.loc[problem]
            #results[c]["var"] = model_fit.bse.loc[problem]**2.
            #results[c]["ci"] = (results[c]["beta"]- (sts.t.ppf(0.975, len(datat)-1) * model_fit.bse.loc[problem]), \
            #    results[c]["beta"]+ (sts.t.ppf(0.975, len(datat)-1) * model_fit.bse.loc[problem]))
        return results



    def simple_mean_diff(self, index, data, problem, positive_class, negative_class, raw, no_covariates):
        datat = data.T
        datat.fillna("NA", inplace=True)
        datat[index] = datat[index].values.astype(float)
        datat["age"] = datat["age"].values.astype(float)
        datat["BMI"] = datat["BMI"].values.astype(float)
        covariates = {\
            "gender": ["age", "BMI"], \
            "age": ["BMI", "gender"], \
            "BMI": ["age", "gender"], \
            "condition": ["BMI", "gender", "age"]}

        datat = datat[[index, problem] + covariates[problem]]
        i = datat[ index ].values.astype(float)
 
        i[ i == 0 ] = 0.00001
        i = np.log2( i )
        datat[ index ] = i

        if not no_covariates:
            formula = ("%s ~ " %index) + " + ".join([problem] + covariates[problem])
        else:
            formula = ("%s ~ %s" %(index, problem))
 
        md = smf.ols( formula, data=datat )
        model_fit = md.fit()
        print(model_fit.summary())
        results = {}
        for c,prob in zip(["BMI", "gender", "age", "condition"], ["BMI", "%s[T.%s]" %("gender", "male"), "age", "%s[T.%s]" %("condition", self.positive)]):
            results[c] = {}
            if not prob in model_fit.pvalues.index:
                results[c]["wald"] = 1.0
                results[c]["beta"] = np.nan
                results[c]["var"] = np.nan
                results[c]["ci"] = (np.nan, np.nan)
            else:
                results[c]["wald"] = model_fit.pvalues.loc[prob]
                results[c]["beta"] = model_fit.params.loc[prob]
                results[c]["var"] = model_fit.bse.loc[prob]**2.
                results[c]["ci"] = (results[c]["beta"]- (sts.t.ppf(0.975, len(datat)-1) * model_fit.bse.loc[prob]), \
                    results[c]["beta"]+ (sts.t.ppf(0.975, len(datat)-1) * model_fit.bse.loc[prob]))
        return results




    def std_mean_diff(self, index, data, problem, positive_class, negative_class, raw, no_covariates):
        datat = data.T
        datat.fillna("NA", inplace=True)
        #datat = datat.loc[datat["age"] != "NA"]
        #datat = datat.loc[datat["gender"] != "NA"]
        #datat = datat.loc[datat["BMI"] != "NA"]
        datat[index] = datat[index].values.astype(float)
        datat["age"] = datat["age"].values.astype(float)
        datat["BMI"] = datat["BMI"].values.astype(float)

        #problem = "age"
        covariates = {\
	    "gender": ["age", "BMI"], \
	    "age": ["BMI", "gender"], \
	    "BMI": ["age", "gender"], \
	    "condition": ["BMI", "gender", "age"]}
        datat = datat[[index, problem] + covariates[problem]]

        if not raw:
            datat[ index ] = np.arcsin(np.sqrt( datat[index].values.astype(float) / 100. ))
        else:
            datat[ index ] = datat[ index ].values.astype(float)

        if not no_covariates:
            formula = ("%s ~ " %index) + " + ".join([problem] + covariates[problem])
        else:
            formula = ("%s ~ %s" %(index, problem))

        md = smf.ols( formula, data=datat )
        model_fit = md.fit()
        ## print(model_fit.summary())

        t = model_fit.tvalues.loc["%s[T.%s]" %(self.condition, self.positive)] #* (-1.)
        n1 = float(len(datat.loc[(datat[problem]==negative_class)]))
        n2 = float(len(datat.loc[(datat[problem]==positive_class)]))
        wilco = sts.ranksums(\
            datat.loc[datat[problem].isin([positive_class]), index].values.astype(float), \
            datat.loc[datat[problem].isin([negative_class]), index].values.astype(float))[1]
        wald = model_fit.pvalues.loc["%s[T.%s]" %(self.condition, self.positive)]
        d = (t*(n1+n2))/float(np.sqrt(n1*n2)*np.sqrt(n1+n2-2))
        SEd = np.sqrt(((n1+n2-1)/float(n1+n2-3)) * ((4./float(n1+n2))*(1+((d**2.)/8.))))

        d_lw = d- (sts.t.ppf(0.975, n1+n2-1) * SEd)
        d_up = d+ (sts.t.ppf(0.975, n1+n2-1) * SEd)
        return d, wald, wilco, (n1, n2), 0.0, "N", SEd



    def compute_auc(self, cohort_frame, oral_index):
        print("Starting AUC computation...", end="")
        observed = [(1.0 if c==self.positive else 0.0) for c in cohort_frame.loc[self.condition]]
        predicted = cohort_frame.loc[oral_index].values.astype(float) 
        #print(predicted, " VAFFANCULO")
        #predicted /= np.sum(predicted)
        #predicted = [np.arcsin(np.sqrt(p) if np.isfinite(p) else 0.0) for p in predicted]
        #[(np.arcsin(np.sqrt( a/np.sum(cohort_frame.loc[oral_index].values.astype(float))))) 
        #for a in cohort_frame.loc[oral_index].values.astype(float)]  #.apply(lambda a : [np.arcsin(np.sqrt( aa/np.sum(a))) for aa in a])
        print("Finished.")
        aucs = roc_auc_score(observed, predicted)
        print(aucs)
        return aucs
 
    #def fold_change(self, cohort_frame, oral_index):

    def wilcoxon(self, cohort_frame, oral_index):
        print("Starting Stat Test comput...", end="")
        positives = cohort_frame.loc[oral_index, cohort_frame.loc[self.condition]==self.positive]
        negatives = cohort_frame.loc[oral_index, cohort_frame.loc[self.condition]==self.negative]
        print("Finished.")
        return sts.ranksums(positives, negatives)[1]

    def score_arrays(self, cohort_frame, oral_index):
        positives = cohort_frame.loc[oral_index, cohort_frame.loc[self.condition]==self.positive]
        negatives = cohort_frame.loc[oral_index, cohort_frame.loc[self.condition]==self.negative]
        return positives,negatives

    def map_2_diseases(self):
        dats = [\
            "XieH_2016_asthma", "XieH_2016_migraine", \
            "ZhuF_2020", "JieZ_2017", "QinN_2014", \
            "RubelMA_2020", "YeZ_2018", \
            "NagySzakalD_2017", "HMP_2019_ibdmdb_UC", \
            "NielsenHB_2014_UC", "HMP_2019_ibdmdb_CD", \
            "NielsenHB_2014_CD", \
            "QinJ_2012", "KarlssonFH_2013", \
            "SankaranarayananK_2015", "ZellerG_2014", \
            "YuJ_2015", "FengQ_2015", \
            "VogtmannE_2016", "ThomasAM_2019_a", \
            "ThomasAM_2019_b", "WirbelJ_2018", \
            "GuptaA_2019", "HanniganGD_2017", \
            "YachidaS_2019"]
        conds = [\
            "asthma", "migraine", "schizofrenia", "ACVD", "cirrhosis", "STH", \
            "BD", "ME/CFS", "UC", "UC", "CD", "CD", "T2D", "T2D", \
            "T2D", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC"]
        mapper = dict([(dt,cn) for dt,cn in zip(dats,conds)])
        return mapper
   

    def regr_validation(self, an_oral_index):
        dats = self.input.loc[self.args.stratify].unique().tolist()

        dat_2_frame = dict([(dt, self.oral_estimate.loc[:, self.oral_estimate.loc[self.args.stratify]==dt]) for dt in dats])

        out_metaanalysis_file = an_oral_index + "_METAANALYSIS_" + self.condition

        #print("FIRST AGE WITH YES COVARIATES")
        self.meta_analysis(dat_2_frame, an_oral_index, dats, self.condition, None, \
            None, out_metaanalysis_file, self.args.type_of_meta, False, False)

        print("SECOND AGE WITH NO COVARIATES (BOTH arcsin)")
        self.meta_analysis(dat_2_frame, an_oral_index, dats, self.condition, None, \
            None, out_metaanalysis_file + "_REV_CRUDE", self.args.type_of_meta, False, True)

        self.meta_interpretable_2(dat_2_frame, an_oral_index, dats, self.condition, None, \
            None, out_metaanalysis_file + "_REV_ADJUSTED", self.args.type_of_meta, raw=True, no_covariates=False)

 

    def get_numbers_binary(self, single_datasets, index, datasets, problem, positive, negative):

        def log2(array):
            array[ array == 0 ] = 0.00001
            return np.log2(array)

        avgs = {"data": [], "lilliefors": [], positive: [], negative:[]}

        for dat in datasets:
            dat_here = single_datasets[dat]
            pos = dat_here.loc[index, dat_here.loc[problem]==positive].values.astype(float)
            neg = dat_here.loc[index, dat_here.loc[problem]==negative].values.astype(float)

            avgs["data"].append(dat)
            avgs["lilliefors"].append( (lilliefors(log2(pos)), lilliefors(log2(neg))) )
            avgs[positive].append(np.mean(pos))
            avgs[negative].append(np.mean(neg))

        for ll,a,b,c in zip(avgs["lilliefors"], avgs["data"], avgs[positive], avgs[negative]):
            print(ll, a, b, c)

        #td = 
        #tc = 

        # (sts.t.ppf(0.975, n1+n2-1) * SEd)
        #print("TOT DISEASE (mean/median)", np.mean(avgs[positive]), np.median(avgs[positive])   ) 
        #print("TOT CONTROL (mean/median)", np.mean(avgs[negative]), np.median(avgs[negative])   )







        
    def sex_validation(self, an_oral_index):
        dats = self.input.loc[self.args.stratify].unique().tolist()
        dat_2_frame = dict([(dt, self.oral_estimate.loc[:, self.oral_estimate.loc[self.args.stratify]==dt]) for dt in dats])
        out_metaanalysis_file = an_oral_index + "_METAANALYSIS_" + self.condition
        AUCs = dict([(cohort, self.compute_auc(dat_2_frame[cohort], an_oral_index)) for cohort in dats])
        for c in AUCs:
            print(c, " ==> ", AUCs[c])
        self.meta_analysis(dat_2_frame, an_oral_index, dats, self.condition, self.positive, \
            self.negative, out_metaanalysis_file, "CLS")

    def main_validation(self, an_oral_index):
        segregated = self.segregate_datasets()

        non_crc_datasets = [\
            "XieH_2016_asthma", "XieH_2016_migraine", \
            "ZhuF_2020", "JieZ_2017", "QinN_2014", \
            "RubelMA_2020", "YeZ_2018", \
            "NagySzakalD_2017", "HMP_2019_ibdmdb_UC", \
            "NielsenHB_2014_UC", "HMP_2019_ibdmdb_CD", \
            "NielsenHB_2014_CD", \
            "QinJ_2012", "KarlssonFH_2013", \
            "SankaranarayananK_2015"]

        crc_datasets = [\
            "ZellerG_2014", "YuJ_2015", "FengQ_2015", \
            "VogtmannE_2016", "ThomasAM_2019_a", \
            "ThomasAM_2019_b", "WirbelJ_2018", \
            "GuptaA_2019", "HanniganGD_2017", \
            "YachidaS_2019"]
        
        AUCs = dict([(cohort, self.compute_auc(segregated[cohort], an_oral_index )) \
	    for cohort in (non_crc_datasets + crc_datasets)])

        Pvals = dict([(cohort, self.wilcoxon(segregated[cohort], an_oral_index)) \
	    for cohort in (non_crc_datasets + crc_datasets)])

        out_metaanalysis_file = os.path.join("meta_analysis_files/", an_oral_index + "_metaanalysis_" + self.condition)
        out_metaanalysis_file_uncorrected = os.path.join("meta_analysis_files/", an_oral_index + "_metaanalysis_uncorrected_" + self.condition)
        
        out_metaanalysis_file_raw = os.path.join("meta_analysis_files/", an_oral_index + "_metaanalysis_rawabunds_" + self.condition)
        out_metaanalysis_file_fold_change = os.path.join("meta_analysis_files/", an_oral_index + "_foldchange_metaanalysis_" + self.condition)

        meandiff_metaanalysis_file_raw = os.path.join("meta_analysis_files/", an_oral_index + "_metaanalysis_meandiff_raw_" + self.condition)
        meandiff_metaanalysis_file_adjusted = os.path.join("meta_analysis_files/", an_oral_index + "_metaanalysis_meandiff_adjusted_" + self.condition)

                            #nsingle_datasets, index, datasets, problem, positive, negative
        self.get_numbers_binary(segregated, an_oral_index, crc_datasets + non_crc_datasets, self.condition, self.positive, \
            self.negative  )

        self.meta_analysis(segregated, an_oral_index, crc_datasets + non_crc_datasets, self.condition, self.positive, \
            self.negative, out_metaanalysis_file_raw, self.args.type_of_meta, raw=True, no_covariates=True)

        self.meta_analysis(segregated, an_oral_index, crc_datasets + non_crc_datasets, self.condition, self.positive, \
            self.negative, out_metaanalysis_file, self.args.type_of_meta, raw=False, no_covariates=False)

        self.meta_analysis(segregated, an_oral_index, crc_datasets + non_crc_datasets, self.condition, self.positive, \
            self.negative, out_metaanalysis_file_uncorrected, self.args.type_of_meta, raw=False, no_covariates=True)
 
        self.log_fold_change_meta_an( segregated, an_oral_index, crc_datasets + non_crc_datasets, \
            self.condition, self.positive, self.negative, out_metaanalysis_file_fold_change, self.args.type_of_meta)
         
        self.meta_interpretable(segregated, an_oral_index, crc_datasets + non_crc_datasets, self.condition, self.positive, \
            self.negative, meandiff_metaanalysis_file_raw, self.args.type_of_meta, raw=True, no_covariates=True    )

        self.meta_interpretable(segregated, an_oral_index, crc_datasets + non_crc_datasets, self.condition, self.positive, \
            self.negative, meandiff_metaanalysis_file_adjusted, self.args.type_of_meta, raw=True, no_covariates=False    )
 
        self.meta_interpretable_2(segregated, an_oral_index, crc_datasets + non_crc_datasets, self.condition, self.positive, \
             self.negative, meandiff_metaanalysis_file_adjusted, self.args.type_of_meta, raw=True, no_covariates=False)

        DATS = crc_datasets + non_crc_datasets
        score_arrays = dict([(cohort, self.score_arrays(segregated[cohort], an_oral_index)) for cohort in DATS])

        self.build_the_main_figure(\
	    segregated, DATS, crc_datasets, score_arrays, an_oral_index, \
	    AUCs, Pvals, out_metaanalysis_file, self.args.type_of_meta)




    def fold_change(self, index, name, data, problem, positive_class, negative_class, std_mean_diff):
        datat = data.T
    
        datat = data.T
        datat[index] = datat[index].values.astype(float)

        disease = datat.loc[datat[problem].isin([positive_class]), index].values.astype(float)
        control = datat.loc[datat[problem].isin([negative_class]), index].values.astype(float)

        disease[ disease == 0 ] = 0.00001
        control[ control == 0 ] = 0.00001
 
        gmean_disease = np.mean( np.log2( disease ))
        gmean_control = np.mean( np.log2( control ))

        if std_mean_diff > 0.0:
            log_fold_change = gmean_disease - gmean_control
        else:
            log_fold_change = gmean_control - gmean_disease
 
        log_fold_change_variance = np.var(np.log2( disease ), ddof=1) + np.var(np.log2( control ), ddof=1)
        se = np.sqrt( log_fold_change_variance )
 
        ci_lw = log_fold_change - (se * sts.t.ppf(0.975, len(datat)-1))
        ci_up = log_fold_change + (se * sts.t.ppf(0.975, len(datat)-1))
 
        print( name, " ==> ", log_fold_change, "95% CI [", ci_lw, ci_up, "]" )

        return log_fold_change, log_fold_change_variance




    def log_fold_change_meta_an(self, single_datasets, INDEX, datasets, \
            problem, positive_class, negative_class, out_metaanalysis_file, TYPE, smd=1.0):

        print("********************* FOLD CHANGES META ANALYSIS **************************")

        study2lens = {}
        singleStudiesClass = []

        for dat in datasets:
            fc, fc_var = self.fold_change( INDEX, dat, single_datasets[dat], problem, positive_class, negative_class, smd )
            study2lens[dat] = {"fc": fc, "fc_var": fc_var}
            sys.stdout.write("%.4f\n" %(fc))

        W = np.array([(1./study2lens[dat]["fc_var"]) for e in datasets], dtype=np.float64)
        effs = np.array([study2lens[dat]["fc"] for e in datasets], dtype=np.float64)

        PM, _ = paule_mandel_tau(\
            np.array([study2lens[dat]["fc"] for e in datasets], dtype=np.float64), \
            np.array([study2lens[dat]["fc_var"] for e in datasets], dtype=np.float64))
 
        W = [(1./float(v+PM)) for v in np.array([study2lens[dat]["fc_var"] for e in datasets], dtype=np.float64)]

        log_eff = np.sum(W*effs)/float(np.sum(W)) 
        sigma_log = 1/float(np.sum(W))

        log_se = np.sqrt(sigma_log) / np.sqrt(len(datasets))

        MU = 2**log_eff
        sigma_fc = 2**sigma_log
        se = 2**log_se

        log_ci = log_eff- (log_se*sts.t.ppf(0.975, len(datasets)-1)), log_eff+ (log_se*sts.t.ppf(0.975, len(datasets)-1))
        ci = [(2**x) for x in list( log_ci ) ]

        print( "AVERAGE\tse\tCI-lw\tCI-up"  )
        print(MU, se, ci[0], ci[1])

        print("************************************* finishing *************************************")
        return MU, se, ci



    def meta_interpretable_2(self, single_datasets, INDEX, datasets, \
        problem, positive_class, negative_class, out_metaanalysis_file, \
        TYPE, raw, no_covariates):
        
        print("************************** meta-analysis of direct responses of score (raw: %s, NO covariates: %s) ******************"%(raw, no_covariates))
        study2lens = {}
        singleStudiesClass = []
        problem = "age"

        for dat in datasets:
            study2lens[dat] = self.simple_mean_diff_rev(INDEX, single_datasets[dat], problem, positive_class, negative_class, raw, no_covariates)

        for k in ["adjusted"]: #, "crude"]:
            print("******* META ANALYSIS OF INTERPRETABLE ==> %s (raw: %s, NO covariates: %s)" %(k, raw, no_covariates))
            for dat in datasets:
                print("Estimating data %s [ %s ]  ===>" %(dat, problem), \
                    study2lens[dat][k]["beta"], study2lens[dat][k]["ci"][0], study2lens[dat][k]["ci"][1], "    (", study2lens[dat][k]["wald"], ")")
  
            W = np.array([(1./study2lens[dat][k]["var"]) for e in datasets], dtype=np.float64)
            effs = np.array([study2lens[dat][k]["beta"] for e in datasets], dtype=np.float64)
 
            PM, _ = paule_mandel_tau(\
                np.array([study2lens[dat][k]["beta"] for e in datasets], dtype=np.float64), \
                np.array([study2lens[dat][k]["var"] for e in datasets], dtype=np.float64))
 
            W = np.array([(1./float(v+PM)) for v in [study2lens[dat][k]["var"] for e in datasets]], dtype=np.float64)
            re = np.sum(W*effs)/float(np.sum(W))
            re_var = 1/float(np.sum(W))
            re_se = np.sqrt(re_var)

            Z = re/re_se
            print(re, re_se, Z)
            P = "%.5f" %(2.*(1 - sts.norm.cdf(np.abs(Z))))
 
            print("FIRST; OF THE RANDOM EFFECT\n")
            print(re, "95% CI [", re - (1.96 * re_se), re + (1.96 * re_se), "]", P, "\n")
 
            h = PM
            sum_of_squared_weights = np.sum( W**2. )
            sum_of_weights = np.sum( W )
            V0 = ( ( 1/sum_of_weights )**2. ) * ( (2 * float(len(datasets))) + ( sum_of_squared_weights / ( sum_of_weights ** 2. ) ) )
            z = 1.96 * np.sqrt(V0)
            up = h + z
            lw = np.max([ 0, h - z])
            print("** =>( HETEROGENEITY )", h, " 95% CI [", lw, up, "]")
            print("*******************************************************   finishings ******************************************************")




    def meta_interpretable( self, single_datasets, INDEX, datasets, \
            problem, positive_class, negative_class, out_metaanalysis_file, \
            TYPE, raw, no_covariates):

        print("************************************  !!! META ANALYSIS OF INTERPRETABLE (raw: %s, NO covariates: %s)" %(raw, no_covariates)) 
        study2lens = {}
        singleStudiesClass = []

        for dat in datasets:
            study2lens[dat] = self.simple_mean_diff(INDEX, single_datasets[dat], problem, positive_class, negative_class, raw, no_covariates)
 
        for k in [problem, "gender", "age", "BMI"]:
            print("******* META ANALYSIS OF INTERPRETABLE ==> %s (raw: %s, NO covariates: %s)" %(k, raw, no_covariates))
            for dat in datasets:
                print("Estimating data %s [ %s ]  ===>" %(dat, problem), \
                    study2lens[dat][k]["beta"], study2lens[dat][k]["ci"][0], study2lens[dat][k]["ci"][1], "    (", study2lens[dat][k]["wald"], ")")

            W = np.array([(1./study2lens[dat][k]["var"]) for e in datasets], dtype=np.float64)
            effs = np.array([study2lens[dat][k]["beta"] for e in datasets], dtype=np.float64)
 
            PM, _ = paule_mandel_tau(\
                np.array([study2lens[dat][k]["beta"] for e in datasets], dtype=np.float64), \
                np.array([study2lens[dat][k]["var"] for e in datasets], dtype=np.float64))
  
            W = np.array([(1./float(v+PM)) for v in [study2lens[dat][k]["var"] for e in datasets]], dtype=np.float64)
            re = np.sum(W*effs)/float(np.sum(W))
            re_var = 1/float(np.sum(W))
            re_se = np.sqrt(re_var)
  
            print("FIRST; OF THE RANDOM EFFECT\n")
            print(re, "95% CI [", re - (1.96 * re_se), re + (1.96 * re_se), "]\n")
 
            h = PM
            sum_of_squared_weights = np.sum( W**2. )
            sum_of_weights = np.sum( W )
            V0 = ( ( 1/sum_of_weights )**2. ) * ( (2 * float(len(datasets))) + ( sum_of_squared_weights / ( sum_of_weights ** 2. ) ) )
            z = 1.96 * np.sqrt(V0)
            up = h + z
            lw = np.max([ 0, h - z])
            print("** =>( HETEROGENEITY )", h, " 95% CI [", lw, up, "]")
            print("*******************************************************   finishings ******************************************************")




    def meta_analysis(self, single_datasets, INDEX, datasets, \
	    problem, positive_class, negative_class, out_metaanalysis_file, \
            TYPE, raw, no_covariates):

        study2lens = {}
        singleStudiesClass = []

        #INDEX = INDEX + "_Corrected"

        if TYPE=="CLS":

            print("************************************  META ANALYSIS OF CLASSIFICATION (raw: %s, NO covariates: %s)" %(raw, no_covariates))
            for dat in datasets:
                print("Estimating data %s [ %s ]" %(dat, problem), end="   ==> ")
                cohenD, p_value_cor, p_value_crude, Lens, fake_p, ci, SEd = \
		    self.std_mean_diff(INDEX, single_datasets[dat], problem, positive_class, negative_class, raw, no_covariates)
                
                #print(cohenD, p_value_cor, p_value_crude,  Lens)
                singleStudiesClass += [singleStudyEffect((cohenD, p_value_cor), dat, Lens, False)]
                study2lens[dat] = {}
                study2lens[dat]["cases"] = Lens[1]
                study2lens[dat]["control"] = Lens[0]
                study2lens[dat]["partialP"] = fake_p if fake_p==fake_p else 1.0
                study2lens[dat]["ci"] = ci
                study2lens[dat]["SE"] = SEd
                sys.stdout.write("%.4f         (%.3f)\n" %(cohenD, p_value_cor))

            considerable = [e for e in singleStudiesClass if (e.accepted)]

            if len(considerable) >= 4: #(len(self.datasets) / 4.):
                re = RE_meta_binary(\
                    [e.effect for e in considerable], \
                    [e.Pvalue for e in considerable], \
                    [e.Name for e in considerable], \
                    [study2lens[e.Name]["cases"] for e in considerable], \
                    [study2lens[e.Name]["control"] for e in considerable], \
                    INDEX, "D", False, False, "PM") #, \

                sys.stdout.write("%s (%i studies) Random Effect = %.3f (%.3f) \n" \
                    %(INDEX, len(considerable), re.RE, re.Pval))
 
                print("FIRST; OF THE RANDOM EFFECT")
                print(re.RE, "95% CI [", re.RE - (1.96 * re.stdErr), re.RE + (1.96 * re.stdErr), "]")

                h = re.t2_PM
                Weights = np.array([   (1./((study2lens[e.Name]["SE"]**2.)+h)) for e in considerable] , dtype=np.float64)
                sum_of_squared_weights = np.sum( Weights**2. )
                sum_of_weights = np.sum( Weights )
                V0 = ( ( 1/sum_of_weights )**2. ) * ( (2 * float(len(considerable))) + ( sum_of_squared_weights / ( sum_of_weights ** 2. ) ) )
                z = 1.96 * np.sqrt(V0)
                up = h + z
                lw = np.max([ 0, h - z])
                print("** =>( HETEROGENEITY )", h, " 95% CI [", lw, up, "]")

                for x,ci in zip(\
                    [y for y in re.result.columns.tolist() if \
			((y.endswith("_CohenD")) and (not y.startswith("RE_")))], re.CI_of_d):
                    re.result[x+"_conf_int"] = ";".join(map(str, ci))

                print("*******************************************************   finishings ******************************************************")
                re.result.to_csv(out_metaanalysis_file+".tsv", sep="\t", header=True, index=True)

        elif TYPE in ["REG", "AGE"]:
            for dat in datasets:
                print("Estimating data %s [ %s ]" %(dat, problem), end="   ==> ")
                PC, p_value, Len, fake_p, ci = self.regression(INDEX, single_datasets[dat], problem, no_covariates)
                print(PC, p_value, Len)
                singleStudiesClass += [singleStudyEffect((PC, p_value), dat, Len, True)]
                considerable = [e for e in singleStudiesClass if (e.accepted)]
            re = RE_meta(\
                [e.effect for e in considerable], \
                [e.Pvalue for e in considerable], \
                [e.Name for e in considerable], \
                [e.Len for e in considerable], INDEX, \
                het="PM", REG=True)
            sys.stdout.write("%s (%i studies) Random Effect = %.3f [%.3f, %.3f]  (%.3f)\n" \
                %(problem, len(considerable), re.RE, re.RE - (1.96 * re.stdErr), re.RE + (1.96 * re.stdErr), re.Pval))
            for x,ci in zip(\
                [y for y in re.result.columns.tolist() if \
		((y.endswith("_Correlation")) and (not y.startswith("RE_")))], re.CI_of_z):
                re.result[x+"_conf_int"] = ";".join(map(str, ci))
            re.result.to_csv(out_metaanalysis_file+".tsv", sep="\t", header=True, index=True)

     


    def build_the_main_figure(self, segregated, DATS, crc_DATS, score_arrays, \
	    index_name, AUCs, Pvals, meta_an, metaan_type):

        index_name2 = index_name #+ "_Corrected"
        ## Pvals : crude
        ## meta_an : corrected

        sns.set_style("whitegrid")
        fig = plt.figure(figsize=(24, 17)) #, constrained_layout=False)

        gs = gridspec.GridSpec(2, 17, hspace=0.5)

        ax_a = plt.subplot(gs[ 0,:10 ])
        ax_b = plt.subplot(gs[ 1,:15 ])
        ax_c = plt.subplot(gs[ 0, 13:-2 ])
        ax_d = plt.subplot(gs[ 1, 16: ])


        ### A
        super_frame = []
        for d in crc_DATS:
            sc_pos = score_arrays[d][0]
            sc_neg = score_arrays[d][1]
            frame_crc = pd.DataFrame(\
                {index_name: list(sc_pos)+list(sc_neg), \
                 "DATASET": \
                 [d for q in range(len(sc_pos) + len(sc_neg))], \
                 "CLASS": \
                 [self.positive for i in range(len(sc_pos))] + [self.negative for j in range(len(sc_neg))], \
                }, index=[d]*(len(sc_pos) + len(sc_neg)))
            if len(super_frame)==0:
                super_frame = frame_crc
            else:
                super_frame = super_frame.append(frame_crc)

        panel_a = sns.boxplot(data=super_frame, x="DATASET", y=index_name, ax=ax_a, \
            hue="CLASS", palette={self.positive: "orange", self.negative: "dodgerblue"})

        ax_a.legend_.remove()
        panel_a.set_yscale("log")

        pvals = [ (("=%.2f" %Pvals[dat]) if Pvals[dat]>0.01 else "<0.01") for dat in crc_DATS]

        ax_a.set_xticklabels([("%s\nAUC=%.2f, P%s" %(dat, AUCs[dat], pval)) \
	    for dat,pval in zip(crc_DATS, pvals)], rotation=90)

        for _,s in ax_a.spines.items():
            s.set_linewidth(1.5)
            s.set_color("black")


        ### B
        super_frame = []
        for d in DATS:
            if(not d in crc_DATS):
                sc_pos = score_arrays[d][0]
                sc_neg = score_arrays[d][1]
                frame_b = pd.DataFrame(\
                   {index_name: list(sc_pos)+list(sc_neg), \
                    "DATASET": \
                    [d for q in range(len(sc_pos) + len(sc_neg))], \
                    "CLASS": \
                    [self.positive for i in range(len(sc_pos))] + [self.negative for j in range(len(sc_neg))], \
                   }, index=[d]*(len(sc_pos) + len(sc_neg)))
                if len(super_frame)==0:
                    super_frame = frame_b
                else:
                    super_frame = super_frame.append(frame_b)

        panel_b = sns.boxplot(data=super_frame, x="DATASET", y=index_name, ax=ax_b, \
           hue="CLASS", palette={self.positive: "orange", self.negative: "dodgerblue"})

        panel_b.set_yscale("log")
        ax_b.legend_.remove()
        for _,s in ax_b.spines.items():
            s.set_linewidth(1.5)
            s.set_color("black")

        nDATS = [d for d in DATS if (not d in crc_DATS)]
        pvals = [ (("=%.2f" %Pvals[dat]) if Pvals[dat]>0.01 else "<0.01") for dat in nDATS]
        ax_b.set_xticklabels([("%s (%s)\nAUC=%.2f, P%s" %(\
	    d, self.mapping_with_diseases[d], AUCs[d], pval)) for d,pval in zip(nDATS,pvals)], rotation=90) #, ha="right")


        ### C
        #meta_an1 = "Oral_Entropy_metaanalysis_uncorrected_condition.tsv"
        #meta_an2 = "Oral_Entropy_metaanalysis_condition.tsv"

        meta_an1 = "meta_analysis_files/Oral_Richness_metaanalysis_uncorrected_condition.tsv"
        meta_an2 = "meta_analysis_files/Oral_Richness_metaanalysis_condition.tsv"

        def plot(meta_an, marker, palette, ax=ax_c):
 
            MA = pd.read_csv(meta_an, sep="\t", header=0, index_col=0)
 
            read_ma = lambda MA,index_name,d : \
                tuple(list(map(float, \
                MA.loc[index_name, [d+"_CohenD", d+"_Pvalue"]].tolist() + \
                str(MA.loc[index_name, d+("_CohenD" if d!="RE" else "")+"_conf_int"]).split(";"))))
            MA_results = [read_ma(MA, index_name2, dat) for dat in DATS + ["RE"]]
 
            MA_frame = {"Dataset": [], \
	        "Effect-size": [], "Lower-Conf": [], "Upper-Conf": [], \
                "P-value": [], "P-value-crude": [], "Enumerate": [], "Significant": []}

            for e,d in enumerate(DATS + ["RE"]):
                MA_frame["Effect-size"] += [MA_results[e][0]]
                MA_frame["Lower-Conf"] += [MA_results[e][2]]
                MA_frame["Upper-Conf"] += [MA_results[e][3]]
                MA_frame["P-value"] += [MA_results[e][1]]
                MA_frame["P-value-crude"] += [Pvals[d] if d!="RE" else 0.0]
                MA_frame["Dataset"] += [d]
                MA_frame["Enumerate"] += [e]
                MA_frame["Significant"] += ["Y" if (MA_results[e][1]<0.05) else "N"]
 
            MA_frame = pd.DataFrame(MA_frame, index=DATS + ["RE"])
            MA_frame = MA_frame.loc[MA_frame.index.tolist()[::-1]]
 
            panel_c = sns.scatterplot(data=MA_frame, x="Effect-size", y="Dataset", \
               hue="Significant", size="Significant", palette=palette, \
               sizes={"Y": 90, "B": 90, "N": 45}, ax=ax, legend=False, marker=marker) ##, order=(DATS + ["RE"])[::-1])
 
            for y_lev,dat in zip(ax.get_yticks()[::-1], MA_frame.index.tolist()[::-1]): 
                MA_frame.loc[dat, "Upper-Conf"]
                MA_frame.loc[dat, "Lower-Conf"]
 
            if float(MA_frame.loc[dat, "P-value"])<0.05:
                clr = "goldenrod"
            else:
                clr = "slategrey"

            ax.plot(
                [float(MA_frame.loc[dat, "Lower-Conf"]), float(MA_frame.loc[dat, "Upper-Conf"])], \
                [y_lev, y_lev], \
                linewidth=2, \
                color=clr, linestyle="-"\
                )

            for _,s in ax.spines.items():
                s.set_linewidth(1.5)
                s.set_color("black") 
 
            ax.set_yticklabels(["Random effect"]+DATS[::-1]) #, rotation=90) #,  ha="")
            ax.set_xlim([-1.75,1.75])
            ax.axvline(c="black", linestyle="--", linewidth=1.5)

            panel_c_tw = ax.twinx()
            panel_c_tw.set_ylim(ax.get_ylim())
            panel_c_tw.set_yticks(ax.get_yticks())
            panel_c_tw.set_yticklabels([("%.2f  [%.2f %.2f]" %\
                (MA_frame.loc[d, "Effect-size"], MA_frame.loc[d, "Lower-Conf"], \
            MA_frame.loc[d, "Upper-Conf"])) for d in (DATS + ["RE"])[::-1]])

        plot(meta_an1, "d", {"Y": "purple", "N": "slategrey"})
        plot(meta_an2, "o", {"Y": "mediumturquoise", "N": "slategrey"})

        leg_handles = [(mlines.Line2D([], [], color=color, marker=marker, linestyle='None', alpha=1.0, \
            markersize=markersize, label=label)) for label,color,marker,markersize in zip(\
            ["Disease", "Control", "Significant (crude)", "Significant (weighted)", "NS", "RE model"], \
            ["orange", "dodgerblue", "goldenrod", "slategrey", "black"], \
            ["s", "s", "o", "o", "D"], [10, 10, 12, 12, 6, 12])]
 
        legend = ax_d.legend(handles=leg_handles, fontsize=12, \
	    loc="center", frameon=True, markerfirst=True, ncol=2)
        ax_d.axis("off")

        plt.subplots_adjust(bottom=0.2, right=0.8)
        [plt.savefig("../images/NEWLY_CLR_TEST_ORAL_SCORE_abc_%s.%s" %(index_name, fmt), dpi=200) for fmt in ["png", "svg"]]



if __name__ == "__main__":
    OI = oral_introgression(sys.argv)
    if OI.args.type_of_meta == "CLS":
        OI.main_validation(OI.args.Oral_Richness)

    elif OI.args.type_of_meta == "AGE": ## python oral_introgression.py 
        OI.regr_validation(OI.args.Oral_Richness)

    elif OI.args.type_of_meta == "SEX":
        OI.sex_validation(OI.args.Oral_Richness)
