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
from metareg import RE_meta_binary, RE_meta, singleStudyEffect

## python oral_introgression.py Oral_Richness ## DISEASES
## python oral_introgression.py Oral_Richness -la age -tm REG -dt datas/age_u_species.tsv 
## python oral_introgression.py Oral_Richness -la BMI -tm REG -dt datas/bmi_u_species.tsv
## python oral_introgression.py Oral_Richness -la gender:male:female -tm SEX -dt datas/sex_u_species.tsv

class oral_introgression(object):
    def readargs(self, args):
        p = ap.ArgumentParser()
        add = p.add_argument
        add("Oral_Richness", type=str)
        add("-dt", "--dataset", type=str, default="nested/the_big_usable.tsv")
        add("-st", "--stratify", type=str, default="dataset_name")
        add("-la", "--lookat", type=str, default="your_problem:positive:negative")
        ####add("-i", "--index", type=str, default="Oral_Ri")
        add("-os", "--oral_species", type=str, default="cMD3_oral_samples.tsv")
        add("-z", "--featid", type=str, default="s__")
        add("-ma", "--min_ab", type=float, default=0.01)
        add("-mp", "--min_prev", type=float, default=5)
        add("-tm", "--type_of_meta", type=str, default="CLS", choices=["CLS", "REG", "SEX"])
        return p.parse_args()

    def __init__(self, args):
        self.args = self.readargs(args)
        self.featid = self.args.featid
        if self.args.type_of_meta in ["CLS", "SEX"]:
            self.condition, self.positive, self.negative = tuple(self.args.lookat.split(":"))
        elif self.args.type_of_meta == "REG":
            self.condition = self.args.lookat
        self.input = pd.read_csv(self.args.dataset, sep="\t", header=0, index_col=0, low_memory=False, engine="c")
        self.oral_cMD = self.define_oral_bases(self.args.oral_species) 
        self.oral_estimate = self.quantify_oral_introgression()
        self.mapping_with_diseases = self.map_2_diseases()
        
    def define_oral_bases(self, oral_cmd):
        if os.path.isfile(oral_cmd):
            return pd.read_csv(oral_cmd, sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")
        else:
            OP = pd.read_csv("March_2021_usable.tsv", sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")
            OP.loc["days_from_first_collection"] = [(0.0 if (str(c)=="NA") else float(c)) for c in OP.loc["days_from_first_collection"]]
            OP = OP.loc[:, ~(OP.loc["days_from_first_collection"]>0.0)]
            OP = OP.loc[:, OP.loc["body_site"]=="oralcavity"]
            OP = OP.T.drop_duplicates(subset="subjectID", keep="first").T
            feats = [i for i in OP.index.tolist() if ((self.args.featid in i) and (np.sum(OP.loc[i].values.astype(float))>0.))]
            OP.index.name = "sampleID"
            OP.to_csv(oral_cmd, sep="\t", header=True, index=True)
            return OP
         
    def quantify_oral_introgression(self):
        print("Initializing oral-score estimation...")
        feats = [i for i in self.oral_cMD.index.tolist() if (self.featid in i)]
        abuns = dict([(ft, np.mean(self.oral_cMD.loc[ft].values.astype(float))) for ft in feats])
        prevs = dict([(ft, (np.count_nonzero(self.oral_cMD.loc[ft].values.astype(float))/float(self.oral_cMD.shape[1]))*100.) for ft in feats])
        oral_taxa = [ft for ft in feats if ((abuns[ft]>self.args.min_ab) and (prevs[ft]>self.args.min_prev))]

        print("The amount of oral taxa you have is ", len(oral_taxa))

        with open("Oral_species_March_21.tsv", "w") as Os:
            for sp in oral_taxa: Os.write(sp + "\n")

        oral_sample = self.oral_cMD.loc[oral_taxa, :]        
   
        def entropy(proportions_array):
            P = np.array(proportions_array, dtype=np.float64)
            log_of_P = [ (np.log(n) if (float(n)!=0.) else 0.0) for n in proportions_array ]
            return -1.*(np.sum(P * np.array(log_of_P, dtype=np.float64)))
 
        def ginismp(proportions_array):
            return np.sum(proportions_array**2.)

        def oral_fraction(tot_proportions_array, or_proportions_array):
            return np.count_nonzero(or_proportions_array) / float(np.count_nonzero(tot_proportions_array))

        metadata = [i for i in self.input.index.tolist() if (not self.featid in i)]
        OSP = [osp for osp in oral_taxa if osp in (self.input.index.tolist())]
        oral_estimate = self.input.copy()
        oral_richness = dict([(samplename, np.sum(oral_estimate.loc[OSP, samplename].values.astype(float)/100.)) \
            for samplename in oral_estimate.columns.tolist()])
        oral_entropy = dict([(samplename, entropy(oral_estimate.loc[OSP, samplename].values.astype(float)/100.)) \
            for samplename in oral_estimate.columns.tolist()])
        oral_ginisimp = dict([(samplename, ginismp(oral_estimate.loc[OSP, samplename].values.astype(float)/100.)) \
            for samplename in oral_estimate.columns.tolist()])
        oral_fraction_ = dict([(samplename, oral_fraction(\
	    oral_estimate.loc[[i for i in oral_estimate.index.tolist() if (self.featid in i)], samplename].values.astype(float)/100., \
	    oral_estimate.loc[OSP, samplename].values.astype(float)/100. \
	    ))  for samplename in oral_estimate.columns.tolist()])

        Get = lambda dic,key : dic[key] if (not np.isnan(dic[key])) else 0.0 

        oral_estimate.loc["Oral_Richness"] = [Get(oral_richness, samplename) for samplename in oral_estimate.columns.tolist()]
        oral_estimate.loc["Oral_Entropy"] = [Get(oral_entropy, samplename) for samplename in oral_estimate.columns.tolist()]
        oral_estimate.loc["Oral_Gini"] = [Get(oral_ginisimp, samplename) for samplename in oral_estimate.columns.tolist()]
        oral_estimate.loc["Oral_Fraction"] = [Get(oral_fraction_, samplename) for samplename in oral_estimate.columns.tolist()]
        print("Finished oral-score estimation.")
        return oral_estimate
      
    def segregate_datasets(self):
        print("Initiating dataset segregation...")
        single_datasets = []
        for dataset in self.oral_estimate.loc["dataset_name"].unique():
            if (not dataset in ["XieH_2016", "HMP_2019_ibdmdb", "NielsenHB_2014"]):
                this_estimate = self.oral_estimate.loc[:, self.oral_estimate.loc["dataset_name"]==dataset]
                key = dataset
                single_datasets += [(key, this_estimate)]
            else:
                if dataset in ["HMP_2019_ibdmdb", "NielsenHB_2014"]:
                    for dis,anti_dis in [("UC", "CD"), ("CD", "UC")]:
                        this_estimate = self.oral_estimate.loc[:, self.oral_estimate.loc["dataset_name"]==dataset]
                        this_estimate = this_estimate.loc[:, this_estimate.loc["disease_subtype"]!=anti_dis]
                        key = dataset + "_" + dis
                        single_datasets += [(key, this_estimate)]
                elif dataset == "XieH_2016":
                    for dis,anti_dis in [("asthma", "migraine"), ("migraine", "asthma")]:
                        this_estimate = self.oral_estimate.loc[:, self.oral_estimate.loc["dataset_name"]==dataset]
                        this_estimate = this_estimate.loc[:, this_estimate.loc["study_condition"]!=anti_dis]
                        key = dataset + "_" + dis
                        single_datasets += [(key, this_estimate)]
        print("Terminated dataset aggregation")
        return dict(single_datasets)

    def regression(self, index, data, problem):
        datat = data.T

        #print(datat)
        datat.fillna("NA", inplace=True)
        datat = datat.loc[datat["age"] != "NA"]
        datat = datat.loc[datat["gender"] != "NA"]
        datat = datat.loc[datat["BMI"] != "NA"]
        datat[index] = datat[index].values.astype(float)
        Lens = len(datat)
        covariates = {\
            "gender": ["age", "BMI"], \
            "age": ["BMI", "gender"], \
            "BMI": ["age", "gender"], \
            "your_problem": ["BMI", "age", "gender"]}
        datat["age"] = datat["age"].values.astype(float)
        datat["BMI"] = datat["BMI"].values.astype(float)
        #datat["gender"] = [(1.0 if(p=="male") else 0.0) for p in datat["gender"].tolist()]
        datat = datat[[index, problem] + covariates[problem]]
        formula = ("%s ~ 1 + " %index) + " + ".join([problem] + covariates[problem])
        md = smf.ols(formula, data=datat)
        model_fit = md.fit()
        t = model_fit.tvalues.loc[problem]
        n = float(len(datat))
        r = float(t) / np.sqrt(np.float((t**2.) + (n - 1.)))
        Zr = np.arctanh(r) #0.5 * np.log((1. + r) / float/(1. - r))
        SEr = 1/np.sqrt(n - 3)
        r_lw = Zr - (1.96*SEr)
        r_up = Zr + (1.96*SEr)
        return np.tanh(r), model_fit.pvalues.loc[problem], n, 0.0, "N"

    def std_mean_diff(self, index, data, problem, positive_class, negative_class):
        datat = data.T
        datat.fillna("NA", inplace=True)
        datat = datat.loc[datat["age"] != "NA"]
        datat = datat.loc[datat["gender"] != "NA"]
        datat = datat.loc[datat["BMI"] != "NA"]
        datat[index] = datat[index].values.astype(float)
        datat["age"] = datat["age"].values.astype(float)
        datat["BMI"] = datat["BMI"].values.astype(float)
        covariates = {\
	    "gender": ["age", "BMI"], \
	    "age": ["BMI", "gender"], \
	    "BMI": ["age", "gender"], \
	    "your_problem": ["BMI", "age", "gender"]}
        datat = datat[[index, problem] + covariates[problem]]
        formula = ("%s ~ 1 + " %index) + " + ".join([problem] + covariates[problem])
        md = smf.ols(formula, data=datat)
        model_fit = md.fit()
        #print(model_fit.summary())
        t = model_fit.tvalues.loc["%s[T.%s]" %(self.condition, self.positive)] * (-1.)
        n1 = float(len(datat.loc[(datat[problem]==negative_class)]))
        n2 = float(len(datat.loc[(datat[problem]==positive_class)]))
        wilco = sts.ranksums(\
            datat.loc[datat[problem].isin([positive_class]), index].values.astype(float), \
            datat.loc[datat[problem].isin([negative_class]), index].values.astype(float))[1]
        wald = model_fit.pvalues.loc["%s[T.%s]" %(self.condition, self.positive)]
        d = (t*(n1+n2))/float(np.sqrt(n1*n2)*np.sqrt(n1+n2-2))
        SEd = np.sqrt(((n1+n2-1)/float(n1+n2-3)) * ((4./float(n1+n2))*(1+((d**2.)/8.))))
        d_lw = d-(1.96*SEd)
        d_up = d+(1.96*SEd)
        return d, wald, wilco, (n1, n2), 0.0, "N"

    def compute_auc(self, cohort_frame, oral_index):
        print("Starting AUC computation...", end="")
        observed = [(1.0 if c==self.positive else 0.0) for c in cohort_frame.loc[self.condition]]
        predicted = cohort_frame.loc[oral_index]
        print("Finished.")
        return roc_auc_score(observed, predicted)
 
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
        self.meta_analysis(dat_2_frame, an_oral_index, dats, self.condition, None, \
            None, out_metaanalysis_file, self.args.type_of_meta)
        
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
        
        AUCs = dict([(cohort, self.compute_auc(segregated[cohort], an_oral_index)) \
	    for cohort in (non_crc_datasets + crc_datasets)])
        Pvals = dict([(cohort, self.wilcoxon(segregated[cohort], an_oral_index)) \
	    for cohort in (non_crc_datasets + crc_datasets)])

        out_metaanalysis_file = an_oral_index + "_METAANALYSIS_" + self.condition
        self.meta_analysis(segregated, an_oral_index, crc_datasets + non_crc_datasets, self.condition, self.positive, \
            self.negative, out_metaanalysis_file, self.args.type_of_meta)

        DATS = crc_datasets + non_crc_datasets
        score_arrays = dict([(cohort, self.score_arrays(segregated[cohort], an_oral_index)) for cohort in DATS])
        self.build_the_main_figure(\
	    segregated, DATS, crc_datasets, score_arrays, an_oral_index, \
	    AUCs, Pvals, out_metaanalysis_file, self.args.type_of_meta)

    def meta_analysis(self, single_datasets, INDEX, datasets, \
	    problem, positive_class, negative_class, out_metaanalysis_file, TYPE):
        study2lens = {}
        singleStudiesClass = []
        if TYPE=="CLS":
            for dat in datasets:
                print("Estimating data %s [ %s ]" %(dat, problem), end="   ==> ")
                cohenD, p_value_cor, p_value_crude, Lens, fake_p, ci = \
		    self.std_mean_diff(INDEX, single_datasets[dat], problem, positive_class, negative_class)
                print(cohenD, p_value_cor, p_value_crude,  Lens)
                singleStudiesClass += [singleStudyEffect((cohenD, p_value_cor), dat, Lens, False)]
                study2lens[dat] = {}
                study2lens[dat]["cases"] = Lens[1]
                study2lens[dat]["control"] = Lens[0]
                study2lens[dat]["partialP"] = fake_p if fake_p==fake_p else 1.0
                study2lens[dat]["ci"] = ci
            considerable = [e for e in singleStudiesClass if (e.accepted)]
            if len(considerable) >= 4: #(len(self.datasets) / 4.):
                re = RE_meta_binary(\
                    [e.effect for e in considerable], \
                    [e.Pvalue for e in considerable], \
                    [e.Name for e in considerable], \
                    [study2lens[e.Name]["cases"] for e in considerable], \
                    [study2lens[e.Name]["control"] for e in considerable], \
                    INDEX, "D", False, False, "PM") #, \
                sys.stdout.write("%s (%i studies) Random Effect = %.3f (%.3f) [H: %.3f %.3f  %.3f]\n" \
                    %(INDEX, len(considerable), re.RE, re.Pval, re.t2_DL, re.t2_PM, re.I2))
                for x,ci in zip(\
                    [y for y in re.result.columns.tolist() if \
			((y.endswith("_CohenD")) and (not y.startswith("RE_")))], re.CI_of_d):
                    re.result[x+"_conf_int"] = ";".join(map(str, ci))
                re.result.to_csv(out_metaanalysis_file+".tsv", sep="\t", header=True, index=True)
        elif TYPE=="REG":
            for dat in datasets:
                print("Estimating data %s [ %s ]" %(dat, problem), end="   ==> ")
                PC, p_value, Len, fake_p, ci = self.regression(INDEX, single_datasets[dat], problem)
                print(PC, p_value, Len)
                singleStudiesClass += [singleStudyEffect((PC, p_value), dat, Len, True)]
                considerable = [e for e in singleStudiesClass if (e.accepted)]
            re = RE_meta(\
                [e.effect for e in considerable], \
                [e.Pvalue for e in considerable], \
                [e.Name for e in considerable], \
                [e.Len for e in considerable], INDEX, \
                het="PM", REG=True)
            sys.stdout.write("%s (%i studies) Random Effect = %.3f (%.3f) [H: %.3f %.3f  %.3f]\n" \
                %(problem, len(considerable), re.RE, re.Pval, re.t2_DL, re.t2_PM, re.I2))
            for x,ci in zip(\
                [y for y in re.result.columns.tolist() if \
		((y.endswith("_Correlation")) and (not y.startswith("RE_")))], re.CI_of_z):
                re.result[x+"_conf_int"] = ";".join(map(str, ci))
            re.result.to_csv(out_metaanalysis_file+".tsv", sep="\t", header=True, index=True)

     
    def build_the_main_figure(self, segregated, DATS, crc_DATS, score_arrays, \
	    index_name, AUCs, Pvals, meta_an, metaan_type):

        ## Pvals : crude
        ## meta_an : corrected

        sns.set_style("whitegrid")
        fig = plt.figure(figsize=(24, 17)) #, constrained_layout=False)

        gs = gridspec.GridSpec(2, 17, hspace=0.5)
        #gs = fig.add_gridspec(nrows=2, ncols=17, left=0.1, right=0.1, wspace=0.6)
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
            hue="CLASS", palette={self.positive: "darkorange", self.negative: "skyblue"})
        panel_a_swarm = sns.swarmplot(data=super_frame, x="DATASET", y=index_name, ax=ax_a, \
           hue="CLASS", palette={self.positive: "darkorange", self.negative: "skyblue"}, \
           edgecolor="gray", dodge=True, size=6)

        ax_a.legend_.remove()
        #panel_b_swarm_plot()
        #crc_swm = sns.swarmplot(data=data_crc, x="dataset_name", y=index, orient="v", dodge=True, \
        #       hue=self.condition, palette={self.positive: "goldenrod", self.negative: "lightskyblue"}, \
        #       ax=crc_ax, edgecolor="black", size=2, order=crc_datasets)
        panel_a.set_yscale("log")
        panel_a_swarm.set_yscale("log")
        #ax_a.set_yticklabels([])
        #ax_a.set_xticklabels(crc_DATS, rotation=90) #, ha="right")
        #panel_a_tw = ax_a.twiny()
        #panel_a_tw.set_xlim(ax_a.get_xlim())
        pvals = [ (("=%.2f" %Pvals[dat]) if Pvals[dat]>0.01 else "<0.01") for dat in crc_DATS]
        #panel_a_tw.set_xticks(ax_a.get_xticks())
        #panel_a_tw.set_xticklabels([("AUC=%.2f, P%s" %(AUCs[dat], pval)) \
	#    for dat,pval in zip(crc_DATS, pvals)], rotation=90) #, ha="left")
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
           hue="CLASS", palette={self.positive: "darkorange", self.negative: "skyblue"})
        panel_b_swarm = sns.swarmplot(data=super_frame, x="DATASET", y=index_name, ax=ax_b, \
	   hue="CLASS", palette={self.positive: "darkorange", self.negative: "skyblue"}, \
	   edgecolor="gray", dodge=True, size=6)

        panel_b.set_yscale("log")
        panel_b_swarm.set_yscale("log")
        ax_b.legend_.remove()
        for _,s in ax_b.spines.items():
            s.set_linewidth(1.5)
            s.set_color("black")

        nDATS = [d for d in DATS if (not d in crc_DATS)]
        pvals = [ (("=%.2f" %Pvals[dat]) if Pvals[dat]>0.01 else "<0.01") for dat in nDATS]
        ax_b.set_xticklabels([("%s (%s)\nAUC=%.2f, P%s" %(\
	    d, self.mapping_with_diseases[d], AUCs[d], pval)) for d,pval in zip(nDATS,pvals)], rotation=90) #, ha="right")
        #panel_b_tw = ax_b.twiny()
        #pvals = [ (("%.3f" %Pvals[dat]) if Pvals[dat]>0.0001 else "<0.0001") for dat in crc_DATS]

        ### C
        MA = pd.read_csv(meta_an+".tsv", sep="\t", header=0, index_col=0)
        read_ma = lambda MA,index_name,d : \
            tuple(list(map(float, \
            MA.loc[index_name, [d+"_CohenD", d+"_Pvalue"]].tolist() + \
            str(MA.loc[index_name, d+("_CohenD" if d!="RE" else "")+"_conf_int"]).split(";"))))
        MA_results = [read_ma(MA, index_name, dat) for dat in DATS + ["RE"]]
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
            MA_frame["Significant"] += ["B" if (MA_results[e][1]<0.05) else ("Y" if (Pvals[d]<0.05) else "N")]
        MA_frame = pd.DataFrame(MA_frame, index=DATS + ["RE"])
        MA_frame = MA_frame.loc[MA_frame.index.tolist()[::-1]]

        panel_c = sns.scatterplot(data=MA_frame, x="Effect-size", y="Dataset", \
           hue="Significant", size="Significant", palette={"B": "goldenrod", "Y": "red", "N": "slategrey"}, \
           sizes={"Y": 90, "B": 90, "N": 45}, ax=ax_c, legend=False) ##, order=(DATS + ["RE"])[::-1])

        for y_lev,dat in zip(ax_c.get_yticks()[::-1], MA_frame.index.tolist()[::-1]): 
            #nnprint(str(y_lev), dat)
            MA_frame.loc[dat, "Upper-Conf"]
            MA_frame.loc[dat, "Lower-Conf"]
 
            #print("FORSE CHE MI PRENDI PER IL CULO???? ", dat, float(MA_frame.loc[dat, "P-value"]))
            if float(MA_frame.loc[dat, "P-value-crude"])<0.05:
                if float(MA_frame.loc[dat, "P-value"])<0.05:
                    clr = "goldenrod"
                else:
                    clr = "red"
            else:
                clr = "slategrey"

            ax_c.plot(
                [float(MA_frame.loc[dat, "Lower-Conf"]), float(MA_frame.loc[dat, "Upper-Conf"])], \
                [y_lev, y_lev], \
                linewidth=2, \
                color=clr, linestyle="-"\
                )

        for _,s in ax_c.spines.items():
            s.set_linewidth(1.5)
            s.set_color("black") 

        ax_c.set_yticklabels(["Random effect"]+DATS[::-1]) #, rotation=90) #,  ha="")
        ax_c.set_xlim([-1.75,1.75])
        
        #ax_c.axhline(c="")
        ax_c.axvline(c="black", linestyle="--", linewidth=1.5)
        #ax_c.axhline(c="")

        panel_c_tw = ax_c.twinx()
        panel_c_tw.set_ylim(ax_c.get_ylim())
        panel_c_tw.set_yticks(ax_c.get_yticks())
        panel_c_tw.set_yticklabels([("%.2f  [%.2f %.2f]" %\
	  (MA_frame.loc[d, "Effect-size"], MA_frame.loc[d, "Lower-Conf"], \
          MA_frame.loc[d, "Upper-Conf"])) for d in (DATS + ["RE"])[::-1]]) #, rotation=90)  #, ha="center")
        #ax_c.xaxis.set_label_position("")
        #ax_c.xaxis.set_ticks_position("right")
         
        leg_handles = [(mlines.Line2D([], [], color=color, marker=marker, linestyle='None', alpha=1.0, \
            markersize=markersize, label=label)) for label,color,marker,markersize in zip(\
            ["Disease", "Control", "Significant (crude)", "Significant (weighted)", "NS", "RE model"], \
            ["darkorange", "skyblue", "red", "goldenrod", "slategrey", "black"], \
            ["s", "s", "o", "o", "o", "D"], [10, 10, 12, 12, 6, 12])]
 
        legend = ax_d.legend(handles=leg_handles, fontsize=12, \
	    loc="center", frameon=True, markerfirst=True, ncol=2)
        ax_d.axis("off")

        plt.subplots_adjust(bottom=0.2, right=0.8)
        [plt.savefig("images/TEST_ORAL_SCORE_abc_%s.%s" %(index_name, fmt), dpi=200) for fmt in ["png", "svg"]]

if __name__ == "__main__":
    OI = oral_introgression(sys.argv)
    if OI.args.type_of_meta == "CLS":
        OI.main_validation(OI.args.Oral_Richness)
    elif OI.args.type_of_meta == "REG": ## python oral_introgression.py 
        OI.regr_validation(OI.args.Oral_Richness)
