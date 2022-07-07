#!/usr/bin/env python

import time
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
from scipy import stats as sts
import tempfile
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from skbio.stats.composition import ilr_inv,ilr
from skbio.stats.composition import multiplicative_replacement
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import r2_score
from statsmodels.regression.mixed_linear_model import MixedLMParams
import pingouin as pg

from RFbox import RFcls

class meta_an_of_meta_an(object):
    def read_command_line(self, command_line):
        p = ap.ArgumentParser()
        add = p.add_argument
        add("-new", "--new", action="store_true", help="Restart meta-analysis from scratch.")
        add("-pn", "--problem_name", type=str, default="condition")
        add("-pa", "--problem_area", type=str, default="disease_population")

        add("-af", "--analysis_folder", type=str, default="ml_dis_rf/")
        add("-mu", "--make_unique", action="store_true")
        add("-lit", "--literals", action="store_true")
        return p.parse_args()



    def __init__(self, args):
        self.args = self.read_command_line(args)
        self.features_id = "s__"
        self.usable = pd.read_csv("../relative_abundances/usable_all_species_with_condition.tsv", sep="\t", header=0, index_col=0, low_memory=False, engine="c")
        self.usable.fillna("NA", inplace=True)
        if "s__[Collinsella]_massiliensis" in self.usable.index:
            self.usable.loc["s__Collinsella_massiliensis"] = self.usable.loc[\
                "s__[Collinsella]_massiliensis"].values.astype(float)
            self.usable.drop("s__[Collinsella]_massiliensis", inplace=True)

        self.feats = [i for i in self.usable.index.tolist() if (self.features_id in i )]
        self.clinical_data = [i for i in self.usable.index.tolist() if (not i in self.feats)]

        self.meta_analysis = {}
        self.usableS = {}
        self.sample_lists = {}
        self.catS = {}
        self.monos = {}
        self.mono_set = set()
        self.meta_set = set()
        self.signs = {}
        self.datasetS = {}

        feats = [f for f in self.usable.index.tolist() if(self.features_id in f)]
        self.set_up_lodo_datasets()
 


    def set_up_lodo_datasets(self):

        #self.datasetS["hypertension"] = []
        #self.datasetS["fatty_liver"] = []
        ## THOMAS c is predicted by:
        # - all CRC (usable_for_lodo_of_ => CRC)
        # - all e basta : self.usable
        # - all non CRC:  
        ## - cv
  
        #self.usable.loc["Condition"] = [("positive" if c!="control" else "negative") for c in self.usable.loc["study_condition"].tolist()]
        
        self.meta_set = {"CRC", "T2D", "CD", "UC"}
        self.mono_set = {"schizofrenia", "asthma", "migraine", "BD", "ACVD", "STH", "cirrhosis", "ME/CFS"}

        self.datasetS = {
		"CRC": [ "FengQ_2015", "GuptaA_2019", "HanniganGD_2017", "ThomasAM_2019_a", \
			 "ThomasAM_2019_b", "VogtmannE_2016", "WirbelJ_2018", "YachidaS_2019", \
			 "ZellerG_2014", "YuJ_2015" ],
		"CD": [ 	"HMP_2019_ibdmdb", "NielsenHB_2014" 			  ],
		"UC": [ "HMP_2019_ibdmdb", "NielsenHB_2014"				  ],
		"T2D": [ "KarlssonFH_2013", "QinJ_2012", "SankaranarayananK_2015"], 
		"ACVD": ["JieZ_2017"],
		"ME/CFS": ["NagySzakalD_2017"],
		"cirrhosis": ["QinN_2014"],
		"STH": ["RubelMA_2020"],
		"migraine": ["XieH_2016"],
		"asthma": ["XieH_2016"],
		"BD": ["YeZ_2018"],
		"schizofrenia": ["ZhuF_2020"]
	}


        for name in self.meta_set:
            usable_of_name = self.usable.loc[:, self.usable.loc["study_name"].isin(self.datasetS[name])]
            feats = [f for f in usable_of_name.index.tolist() if(self.features_id in f)]
            usable_of_name.drop([f for f in feats if(np.sum(usable_of_name.loc[f].values.astype(float))==0.)], inplace=True)
            feats = [f for f in usable_of_name.index.tolist() if(self.features_id in f)]

            if name == "CD":
                usable_of_name = usable_of_name[[c for c,st in zip(usable_of_name.columns.tolist(), usable_of_name.loc["disease_subtype"]) if (st!="UC")]]
                usable_of_name.drop([f for f in feats if(np.sum(usable_of_name.loc[f].values.astype(float))==0.)], inplace=True)

            elif name == "UC":
                usable_of_name = usable_of_name[[c for c,st in zip(usable_of_name.columns.tolist(), usable_of_name.loc["disease_subtype"]) if (st!="CD")]]
                usable_of_name.drop([f for f in feats if(np.sum(usable_of_name.loc[f].values.astype(float))==0.)], inplace=True)

            #elif name == "migraine":
            #    usable_of_name = usable_of_name[[c for c,st in zip(usable_of_name.columns.tolist(), usable_of_name.loc["disease_subtype"]) if (st!="CD")]]
            #elif name == "asthma":
            #    usable_of_name = usable_of_name[[c for c,st in zip(usable_of_name.columns.tolist(), usable_of_name.loc["disease_subtype"]) if (st!="CD")]]
 
            #if name in ["asthma", "migraine", "UC", "CD"]:
            #    usable_of_name.to_csv(os.path.join(self.args.analysis_folder, "usable_for_lodo_of_%s.tsv" %name.replace("/","_")), sep="\t", header=True, index=True)
            #else:

            usable_of_name.to_csv(os.path.join(self.args.analysis_folder, "usable_for_lodo_of_%s.tsv" %name.replace("/","_")), sep="\t", header=True, index=True)

            for dataset in self.datasetS[name]:
                the_other_datasets = [d for d in self.datasetS[name] if (d!=dataset)]
 
                #this_dataset_in_this_diease = d

                usable_only_for_this_dataset = self.usable.loc[:, self.usable.loc["study_name"].isin([\
                    dat for dat in self.usable.loc["study_name"].tolist() if(not dat in the_other_datasets)])]

                feats = [f for f in usable_only_for_this_dataset.index.tolist() if(self.features_id in f)]

                usable_only_for_this_dataset.drop([f for f in feats \
                     if (np.sum(usable_only_for_this_dataset.loc[f].values.astype(float))==0.)], inplace=True)

                usable_only_for_this_dataset.to_csv(os.path.join(self.args.analysis_folder, "usable_for_lodo_of_%s_only_of_%s.tsv" %(name, \
		    dataset if (not name in ["asthma", "migraine", "UC", "CD"]) else dataset + "_" + name)), sep="\t", header=True, index=True)
 
        usable_only_corrected = self.usable.copy() #.loc[:, self.usable.loc["dataset_name"].isin([dat for dat in self.usable.loc["dataset_name"].tolist()])]
        usable_only_corrected.fillna("NA", inplace=True)
 
        ## if self.features_id=="s__":
        ##  usable_only_corrected.loc["gender"] = [(1.0 if(g=="male") else (0.0 if(g=="female") else "NA")) for g in usable_only_corrected.loc["gender"].tolist()]

        usable_only_corrected = usable_only_corrected.loc[:, \
	    ~(usable_only_corrected.loc["BMI"].isin(["NA"]) | usable_only_corrected.loc["gender"].isin(["NA"]) | usable_only_corrected.loc["age"].isin(["NA"]))]
        
        feats = [f for f in usable_only_corrected.index.tolist() if(self.features_id in f)]
        usable_only_corrected.drop([f for f in feats if(np.sum(usable_only_corrected.loc[f].values.astype(float))==0.)], inplace=True)

        usable_only_corrected.loc[self.features_id + "age"] = usable_only_corrected.loc["age"].tolist()
        usable_only_corrected.loc[self.features_id + "BMI"] = usable_only_corrected.loc["BMI"].tolist()
        usable_only_corrected.loc[self.features_id + "gender"] = usable_only_corrected.loc["gender"].tolist()

        usable_only_corrected.to_csv("../relative_abundances/usable_all_species_with_condition.tsv".replace(\
            ".tsv", "only_corrected.tsv"), sep="\t", header=True, index=True)


        for name in self.datasetS:
            for dataset in self.datasetS[name]: 

                usable_of_name = self.usable.loc[:, self.usable.loc["study_name"].isin([ dataset ])]

                if name == "migraine":
                    usable_of_name = usable_of_name[[c for c,st in zip(usable_of_name.columns.tolist(), usable_of_name.loc["study_condition"]) if (st!="asthma")]]
                elif name == "asthma":
                    usable_of_name = usable_of_name[[c for c,st in zip(usable_of_name.columns.tolist(), usable_of_name.loc["study_condition"]) if (st!="migraine")]]
                elif name == "CD":
                    usable_of_name = usable_of_name[[c for c,st in zip(usable_of_name.columns.tolist(), usable_of_name.loc["disease_subtype"]) if (st!="UC")]]
                elif name == "UC":
                    usable_of_name = usable_of_name[[c for c,st in zip(usable_of_name.columns.tolist(), usable_of_name.loc["disease_subtype"]) if (st!="CD")]]

                feats = [f for f in usable_of_name.index.tolist() if(self.features_id in f)]
                usable_of_name.drop([f for f in feats if(np.sum(usable_of_name.loc[f].values.astype(float))==0.)], inplace=True)
                feats = [f for f in usable_of_name.index.tolist() if(self.features_id in f)]

                filename = os.path.join(self.args.analysis_folder, "cv_dataset_%s_and_%s.tsv" %(name.replace("/", "_"), \
			(dataset if (not dataset in ["HMP_2019_ibdmdb", "NielsenHB_2014", "XieH_2016"]) else dataset + "_" + name)))
 
                usable_of_name.to_csv(filename, sep="\t", header=True, index=True)
 






 



    def machine_learning_main(self, printout=True):

        ofolder = "ml_dis_rf/results/"
        ifolder = "ml_dis_rf/"
        usable_name = "usable_all_species_with_condition.tsv"
        usable_name_metadata = "usable_all_species_with_condition_only_corrected.tsv"
        nruns = 10 
        featcuts = "1"
        ncores = 50

        self.meta_set = {"CRC", "T2D", "CD", "UC"}
        self.mono_set = {"schizofrenia", "asthma", "migraine", "BD", "ACVD", "STH", "cirrhosis", "ME/CFS"}

        self.datasetS = {
                "CRC": [ "FengQ_2015", "GuptaA_2019", "HanniganGD_2017", \
			"ThomasAM_2019_a", "ThomasAM_2019_b", "VogtmannE_2016", \
			"WirbelJ_2018", "YachidaS_2019", "ZellerG_2014", "YuJ_2015" ],
                "CD": [         "HMP_2019_ibdmdb", "NielsenHB_2014"                       ],
                "UC": [ "HMP_2019_ibdmdb", "NielsenHB_2014"                               ],
                "T2D": [ "KarlssonFH_2013", "QinJ_2012", "SankaranarayananK_2015"],
                "ACVD": ["JieZ_2017"],
                "ME/CFS": ["NagySzakalD_2017"],
                "cirrhosis": ["QinN_2014"],
                "STH": ["RubelMA_2020"],
                "migraine": ["XieH_2016"],
                "asthma": ["XieH_2016"],
                "BD": ["YeZ_2018"],
                "schizofrenia": ["ZhuF_2020"]
        }
 
        ## LODO for META ONLY IN CATGEORY
        for name in self.meta_set:
            for blockable in self.datasetS[name]:
              if blockable != "ThomasAM_2019_c":
                ##if name in ["asthma", "migraine", "UC", "CD"]:

                input_file = os.path.join( ifolder, "usable_for_lodo_of_%s.tsv" %name)
                problem = ":".join(["1", "condition", "positive"])
                
                if name in ["asthma", "migraine", "UC", "CD"]:
                    out_file = os.path.join( ofolder , "usable_for_lodo_of_%s_on_%s_rf" %(name, blockable + "_" + name))
                else:
                    out_file = os.path.join( ofolder , "usable_for_lodo_of_%s_on_%s_rf" %(name, blockable))

                rf = RFcls(problem, input_file, "O", None, 1000, 5, "0.1", featcuts, ncores, nruns, self.features_id)
                rf.disable_feats_ranking()

                target = "study_name:" + blockable

                if ((not os.path.exists(out_file + ".txt")) or (os.path.getsize(out_file + ".txt") == 0)):
                    rf.run(input_file, out_file, target)      


        ## LODO for META ONLY OUTSIDE CATEGORY
        for name in self.meta_set:  ## QUESTO SEMBRA OK
            for dataset in self.datasetS[name]:
              if dataset != "ThomasAM_2019_c":
                if name in ["asthma", "migraine", "UC", "CD"]:
                    input_file = os.path.join(  ifolder, "usable_for_lodo_of_%s_only_of_%s.tsv" %(name, dataset + "_" + name))
                else:
                    input_file = os.path.join(  ifolder, "usable_for_lodo_of_%s_only_of_%s.tsv" %(name, dataset))

                if name in ["asthma", "migraine", "UC", "CD"]:
                    out_file = os.path.join( ofolder, "usable_for_lodo_of_%s_only_of_%s_rf" %(name, dataset + "_" + name))
                else:
                    out_file = os.path.join(  ofolder, "usable_for_lodo_of_%s_only_of_%s_rf" %(name, dataset))

                problem = ":".join(["1", "condition", "positive"])
                rf = RFcls(problem, input_file, "O", None, 1000, 5, "0.1", featcuts, ncores, nruns, self.features_id)
                rf.disable_feats_ranking()
                target = "study_name:" + dataset
 
                if ((not os.path.exists(out_file + ".txt")) or (os.path.getsize(out_file + ".txt") == 0)):
                    rf.run(input_file, out_file, target)

        ## LODO FOR MONOS and for META
        for dataset in self.usable.loc["study_name"].unique(): # + ["", "", "", ""]:
            input_file = os.path.join( ifolder, usable_name )

            if (not dataset in ["HMP_2019_ibdmdb", "NielsenHB_2014", "XieH_2016"]):
                out_file = os.path.join( ifolder, usable_name ).replace("%s/" %ifolder, "%s/" %ofolder).replace(".tsv", "_on_%s_rf" %dataset)

                problem = ":".join(["1", "condition", "positive"])

                rf = RFcls(problem, input_file, "O", None, 1000, 5, "0.1", featcuts, ncores, nruns, self.features_id)
                rf.disable_feats_ranking()
                target = "study_name:" + dataset

                if ((not os.path.exists(out_file + ".txt")) or (os.path.getsize(out_file + ".txt") == 0)):
                    rf.run(input_file, out_file, target)

            else:
                if dataset in ["HMP_2019_ibdmdb", "NielsenHB_2014"]:
                    for dis,anti_dis in zip(["CD", "UC"], ["UC", "CD"]):
                        out_file = os.path.join(ifolder, usable_name).replace("%s/" %ifolder, "%s/" %ofolder).replace(".tsv", "_on_%s_rf" %(dataset + "_" + dis))

                        iii_input_file = pd.read_csv(input_file, sep="\t", header=0, index_col=0, engine="c", low_memory=False)
                        feats = [f for f in iii_input_file.index.tolist() if(self.features_id in f)]
                        iii_input_file = iii_input_file.loc[:, iii_input_file.loc["disease_subtype"]!=anti_dis]
                        iii_input_file.drop([f for f in feats if(np.sum(iii_input_file.loc[f].values.astype(float))==0.)], inplace=True)
                        iii_input_file.to_csv(  input_file.replace(".tsv", "_%s.tsv") %dis, sep="\t", header=True, index=True  )

                        problem = ":".join(["1", "condition", "positive"])
                        rf = RFcls(problem, input_file.replace(".tsv", "_%s.tsv") %dis, "O", None, 1000, 5, "0.1", featcuts, ncores, nruns, self.features_id)
                        rf.disable_feats_ranking()
                        target = "study_name:" + dataset
 
                        if ((not os.path.exists(out_file + ".txt")) or (os.path.getsize(out_file + ".txt") == 0)):
                            rf.run(input_file.replace(".tsv", "_%s.tsv") %dis, out_file, target)

                elif dataset in ["XieH_2016"]:
                    for dis,anti_dis in [("asthma", "migraine"), ("migraine", "asthma")]:
                        out_file = os.path.join(ifolder, usable_name).replace("%s/" %ifolder, "%s/" %ofolder).replace(".tsv", "_on_%s_rf" %(dataset + "_" + dis))

                        iii_input_file = pd.read_csv(input_file, sep="\t", header=0, index_col=0, engine="c", low_memory=False)
                        feats = [f for f in iii_input_file.index.tolist() if(self.features_id in f)]
                        iii_input_file = iii_input_file.loc[:, iii_input_file.loc["study_condition"]!=anti_dis]
                        iii_input_file.drop([f for f in feats if(np.sum(iii_input_file.loc[f].values.astype(float))==0.)], inplace=True)
                        iii_input_file.to_csv(  input_file.replace(".tsv", "_%s.tsv") %dis, sep="\t", header=True, index=True  )

                        problem = ":".join(["1", "condition", "positive"])
                        rf = RFcls(problem, input_file.replace(".tsv", "_%s.tsv") %dis , "O", None, 1000, 5, "0.1", featcuts, ncores, nruns, self.features_id)
                        rf.disable_feats_ranking()
                        target = "study_name:" + dataset
 
                        if ((not os.path.exists(out_file + ".txt")) or (os.path.getsize(out_file + ".txt") == 0)):
                            rf.run(input_file.replace(".tsv", "_%s.tsv") %dis, out_file, target)


        ### questo e corrected
        for dataset in self.usable.loc["study_name"].unique():   ## QUESTO OK
            input_file = os.path.join( ifolder, usable_name_metadata ) 

            if (not dataset in ["HMP_2019_ibdmdb", "NielsenHB_2014", "XieH_2016"]):
                out_file = input_file.replace("%s/" %ifolder, "%s/" %ofolder).replace(".tsv", "_on_%s_rf" %dataset)
                problem = ":".join(["1", "condition", "positive"])

                rf = RFcls(problem, input_file, "O", None, 1000, 5, "0.1", featcuts, ncores, nruns, self.features_id)
                rf.disable_feats_ranking()
                target = "study_name:" + dataset
 
                if ((not os.path.exists(out_file + ".txt")) or (os.path.getsize(out_file + ".txt") == 0)):
                    rf.run(input_file, out_file, target)

            else:
                if dataset in ["HMP_2019_ibdmdb", "NielsenHB_2014"]:
                    for dis,anti_dis in zip(["CD", "UC"], ["UC", "CD"]):
                        out_file = input_file.replace("%s/" %ifolder, "%s/" %ofolder).replace(".tsv", "_on_%s_rf" %(dataset + "_" + dis))
   
                        iii_input_file = pd.read_csv(input_file, sep="\t", header=0, index_col=0, engine="c", low_memory=False)
                        feats = [f for f in iii_input_file.index.tolist() if(self.features_id in f)]
                        iii_input_file = iii_input_file.loc[:, iii_input_file.loc["disease_subtype"]!=anti_dis]
                        iii_input_file.drop([f for f in feats if(np.sum(iii_input_file.loc[f].values.astype(float))==0.)], inplace=True)
                        iii_input_file.to_csv(  input_file.replace(".tsv", "_%s.tsv") %dis, sep="\t", header=True, index=True  )

                        problem = ":".join(["1", "condition", "positive"])
                        rf = RFcls(problem, input_file.replace(".tsv", "_%s.tsv") %dis, "O", None, 1000, 5, "0.1", featcuts, ncores, nruns, self.features_id)
                        rf.disable_feats_ranking()
                        target = "study_name:" + dataset
 
                        if ((not os.path.exists(out_file + ".txt")) or (os.path.getsize(out_file + ".txt") == 0)):
                            rf.run(input_file.replace(".tsv", "_%s.tsv") %dis, out_file, target)

                elif dataset in ["XieH_2016"]: ## QUESTO OK
                    for dis,anti_dis in [("asthma", "migraine"), ("migraine", "asthma")]:
                        out_file = input_file.replace("%s/" %ifolder, "%s/" %ofolder).replace(".tsv", "_on_%s_rf" %(dataset + "_" + dis))

                        iii_input_file = pd.read_csv(input_file, sep="\t", header=0, index_col=0, engine="c", low_memory=False)
                        feats = [f for f in iii_input_file.index.tolist() if(self.features_id in f)]
                        iii_input_file = iii_input_file.loc[:, iii_input_file.loc["study_condition"]!=anti_dis]
                        iii_input_file.drop([f for f in feats if(np.sum(iii_input_file.loc[f].values.astype(float))==0.)], inplace=True)
                        iii_input_file.to_csv(  input_file.replace(".tsv", "_%s.tsv") %dis, sep="\t", header=True, index=True  )

                        problem = ":".join(["1", "condition", "positive"])
                        rf = RFcls(problem, input_file.replace(".tsv", "_%s.tsv") %dis, "O", None, 1000, 5, "0.1", featcuts, ncores, nruns, self.features_id)
                        rf.disable_feats_ranking()
                        target = "study_name:" + dataset
 
                        if ((not os.path.exists(out_file + ".txt")) or (os.path.getsize(out_file + ".txt") == 0)):
                            rf.run( input_file.replace(".tsv", "_%s.tsv" ) %dis, out_file, target)


        ### SINFGLE CORSSS VALIDATIONS:
        for name in self.datasetS:
            for dataset in self.datasetS[name]:
                input_file = os.path.join( ifolder, "cv_dataset_%s_and_%s.tsv" %(name if (not name.startswith("ME")) else "ME_CFS", \
			(dataset if (not dataset in ["HMP_2019_ibdmdb", "NielsenHB_2014", "XieH_2016"]) else dataset + "_" + name)))

                out_file = (input_file.replace("%s/" %ifolder, "%s/" %ofolder)[:-4]) + "_rf"

                #print( input_file )
                #print( out_file )

                problem = ":".join(["1", "condition", "positive"])
                rf = RFcls(problem, input_file, out_file, None, 1000, 5, "0.1", featcuts, ncores, nruns, self.features_id)
                rf.disable_feats_ranking()

                if ((not os.path.exists(out_file + ".txt")) or (os.path.getsize(out_file + ".txt") == 0)):
                    rf.runCV()


if __name__ == "__main__":
    maoma = meta_an_of_meta_an(sys.argv)
    maoma.machine_learning_main()
