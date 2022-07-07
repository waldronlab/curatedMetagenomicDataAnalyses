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
from scipy import stats as sts
from scipy.stats import rankdata
from scipy.spatial import distance
from scipy.cluster import hierarchy

FEATID = "s__"
#FEATID="K"

class healthy_microbiome_plot(object):
    def get_aucs(self, file):
            nfeat = {1: "all"} #, 12: "10"} #, 23:"80", 34:"120"}
            aucs = []
            f = open(file)
            e = 0
            lines = f.readlines()
            for line in lines:
                if e in nfeat:
                    aucs += [lines[e+7].rstrip().split()[1]] #[nfeat[e]]
                e += 1
            f.close()
            return list(map(float, aucs)) if len(aucs) else [.5]

    def get_species(self):
        with open(self.args.species_list) as sp_ls:
            species_list = [species.rstrip() for species in sp_ls]
        return species_list

    def get_feats(self, filename, dataset_name, line="line", skip=0):
        f = open(filename) #self.ml_cross_results[self.datasets.index(ds)])
        while (not line.startswith("Feature")):
             line, skip = f.readline(), skip+1
        f.close()
        if filename:
            d = pd.read_csv(filename, sep='\t', header=None, index_col=0, skiprows=skip)[[1,2]]
        d[2] = d[2].apply(lambda n : float(n))
        d.columns = ["feat", dataset_name]
        d.set_index("feat", inplace = True)
        d[dataset_name] = [i+1 for i,f in enumerate(d[dataset_name].tolist())]
        return d

    def __init__(self, modality="complete"):

        self.corrected_datasets = ["XieH_2016_asthma", "XieH_2016_migraine", "ZhuF_2020", "JieZ_2017", "QinN_2014", \
                "RubelMA_2020", "YeZ_2018", "NagySzakalD_2017", "HMP_2019_ibdmdb_UC", "NielsenHB_2014_UC", "HMP_2019_ibdmdb_CD", "NielsenHB_2014_CD", \
                "QinJ_2012", "KarlssonFH_2013", "SankaranarayananK_2015", "ZellerG_2014", \
                "YuJ_2015", "FengQ_2015", "VogtmannE_2016", "ThomasAM_2019_a", "ThomasAM_2019_b", \
                "WirbelJ_2018", "GuptaA_2019", "HanniganGD_2017", "YachidaS_2019"]
        self.corrected_names = ["asthma", "migraine", "schizofrenia", "ACVD", "cirrhosis", \
                "STH", "BD", "ME_CFS", "UC", "UC", "CD", "CD", "T2D", "T2D", "T2D", "CRC", \
                "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC"]
  
        self.multiple_datasets = ["HMP_2019_ibdmdb_UC", "NielsenHB_2014_UC", "HMP_2019_ibdmdb_CD", "NielsenHB_2014_CD", \
                "QinJ_2012", "KarlssonFH_2013", "SankaranarayananK_2015", "ZellerG_2014", \
                "YuJ_2015", "FengQ_2015", "VogtmannE_2016", "ThomasAM_2019_a", "ThomasAM_2019_b", \
                "WirbelJ_2018", "GuptaA_2019", "HanniganGD_2017", "YachidaS_2019"]
        self.multiple_names = ["UC", "UC", "CD", "CD", "T2D", "T2D", "T2D", "CRC", \
                "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC"]
 
        self.datasets = ["XieH_2016_asthma", "XieH_2016_migraine", "ZhuF_2020", "JieZ_2017", "QinN_2014", \
                "RubelMA_2020", "YeZ_2018", "NagySzakalD_2017", "HMP_2019_ibdmdb_UC", "NielsenHB_2014_UC", "HMP_2019_ibdmdb_CD", "NielsenHB_2014_CD", \
                "QinJ_2012", "KarlssonFH_2013", "SankaranarayananK_2015", "ZellerG_2014", \
                "YuJ_2015", "FengQ_2015", "VogtmannE_2016", "ThomasAM_2019_a", "ThomasAM_2019_b", \
                "WirbelJ_2018", "GuptaA_2019", "HanniganGD_2017", "YachidaS_2019"]
        self.names = ["asthma", "migraine", "schizofrenia", "ACVD", "cirrhosis", \
                "STH", "BD", "ME_CFS", "UC", "UC", "CD", "CD", "T2D", "T2D", "T2D", "CRC", \
                "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC"]

        self.cross_datasets = ["XieH_2016_asthma", "XieH_2016_migraine", "ZhuF_2020", "JieZ_2017", "QinN_2014", \
                "RubelMA_2020", "YeZ_2018", "NagySzakalD_2017", "HMP_2019_ibdmdb_UC", "NielsenHB_2014_UC", "HMP_2019_ibdmdb_CD", "NielsenHB_2014_CD", \
                "QinJ_2012", "KarlssonFH_2013", "SankaranarayananK_2015", "ZellerG_2014", \
                "YuJ_2015", "FengQ_2015", "VogtmannE_2016", "ThomasAM_2019_a", "ThomasAM_2019_b", \
                "WirbelJ_2018", "GuptaA_2019", "HanniganGD_2017", "YachidaS_2019"]
        self.cross_names = ["asthma", "migraine", "schizofrenia", "ACVD", "cirrhosis", \
                "STH", "BD", "ME_CFS", "UC", "UC", "CD", "CD", "T2D", "T2D", "T2D", "CRC", \
                "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC", "CRC"]

        self.colormap = "cividis"
        self.ab_colormap = "plasma"

        data2name = dict([(x,y) for x,y in zip(self.datasets, self.names)])
        dat2_multiple_names = dict([(x,y) for x,y in zip(self.multiple_datasets, self.multiple_names)])
        dat2_cross_names = dict([(x,y) for x,y in zip(self.cross_datasets, self.cross_names)])
 
        #data2name["HMP_2019_ibdmdb_CD"] = "CD"
        #data2name["NielsenHB_2014_CD"] = "CD"
        #dat2_multiple_names["HMP_2019_ibdmdb_CD"] = "CD"
        #dat2_multiple_names["NielsenHB_2014_CD"] = "CD"
        #dat2_cross_names["HMP_2019_ibdmdb_CD"] = "CD"
        #dat2_cross_names["NielsenHB_2014_CD"] = "CD"

        #data2name["XieH_2016_ashtma"] = "ashtma"
        #dat2_cross_names["XieH_2016_ashtma"] = "ashtma"
  
        def lodo_corrected(dat):
            return os.path.join("ml_dis_rf", "usable_all_species_with_condition_only_corrected_on_%s_rf.txt" %dat)

        def lodo_multiple(dat): 
            return os.path.join("ml_dis_rf/results", "usable_for_lodo_of_%s_on_%s_rf.txt" %(dat2_multiple_names[dat], dat))

        def lodo_direct(dat):
            return os.path.join("ml_dis_rf", "usable_all_species_with_condition_on_%s_rf.txt" %dat)

        def lodo_one_at_the_time(dat): #ML_dis_RF/results/usable_for_lodo_of_CRC_only_of_YuJ_2015_rf.txt
            return os.path.join("ml_dis_rf/results", "usable_for_lodo_of_%s_only_of_%s_rf.txt" %(dat2_multiple_names[dat], dat))

        def cross_valids(dat):
            return os.path.join("ml_dis_rf", "cv_dataset_%s_and_%s_rf.txt" %(dat2_cross_names[dat], dat))

        data_aucs_direct = dict([(name, [self.get_aucs(lodo_direct(name))[0]]) for name in self.datasets])
        data_aucs_multiple_only = dict([(name, [self.get_aucs(lodo_multiple(name))[0]]) for name in self.multiple_datasets])
        data_aucs_one_at_the_time = dict([(name, [self.get_aucs(lodo_one_at_the_time(name))[0]]) for name in self.multiple_datasets])

        data_aucs_cross = dict([(name, [self.get_aucs(cross_valids(name))[0]]) for name in self.cross_datasets])
        data_corrected = dict([(name, [self.get_aucs(lodo_corrected(name))[0]]) for name in self.corrected_datasets])
        
        #data_aucs_multiple_only.update({\
        #    "HMP_2019_ibdmdb_CD": self.get_aucs(lodo_multiple("HMP_2019_ibdmdb_CD"))[0], \
        #    "NielsenHB_2014_CD": self.get_aucs(lodo_multiple("NielsenHB_2014_CD"))[0]})
        
        #data_aucs_one_at_the_time.update({\
        #    "HMP_2019_ibdmdb_CD": self.get_aucs(lodo_one_at_the_time("HMP_2019_ibdmdb_CD"))[0], \
        #    "NielsenHB_2014_CD": self.get_aucs(lodo_one_at_the_time("NielsenHB_2014_CD"))[0]})

        #data_aucs_cross.update({\
        #    "HMP_2019_ibdmdb_CD": self.get_aucs(cross_valids("HMP_2019_ibdmdb_CD"))[0], \
        #    "NielsenHB_2014_CD":  self.get_aucs(cross_valids("NielsenHB_2014_CD"))[0], \
        #    "XieH_2016_ashtma": self.get_aucs(cross_valids("XieH_2016_ashtma"))[0]})
 
        self.feats_frame = self.get_feats(lodo_direct(self.datasets[0]), self.datasets[0])

        def median_of_medians( datasets, dat2_multiple_names ):

            fake_feats_frame = self.get_feats(lodo_direct(self.datasets[0]), self.datasets[0])

            for dt in self.datasets[1:]:
                features = self.get_feats(lodo_direct(dt), dt)
                fake_feats_frame = fake_feats_frame.join( features, how="outer" )
  
            median_multiples = []
            median_multiples = np.array([[np.median(fake_feats_frame.loc[feat, [dt for dt in datasets \
                    if (dat2_multiple_names.get(dt, "NOT_PRESENT")==CLASS)]].values.astype(float)) \
                    for feat in fake_feats_frame.index.tolist()] for CLASS in list(set(list(dat2_multiple_names.values())))\
                    ], dtype=np.float64).T

            #print( list(set(list(dat2_multiple_names.values()))) )
            #print( list(set(list(dat2_multiple_names.values()))) )
            #print( list(set(list(dat2_multiple_names.values()))) )

            for dt in datasets:
                if (not dt in dat2_multiple_names):
                    to_stack = np.array([[fake_feats_frame.loc[ft, dt]] for ft in fake_feats_frame.index.tolist()], dtype=np.float64)
                    median_multiples = np.hstack((median_multiples, to_stack))

            medians = [np.median(row) for row in median_multiples]  #median_multiples[x, :]) for x in median_multiples]
            return medians

        for dt in self.datasets[1:]:
            features = self.get_feats(lodo_direct(dt), dt)
            self.feats_frame = self.feats_frame.join( features, how="outer" )

        Median = median_of_medians( self.datasets, dat2_multiple_names )
        self.feats_frame["median"] = Median
        self.feats_frame.sort_values("median", inplace=True, ascending=True)

        ## Implementa la direzione e l effetto e dividile in base ad essa
        #if modality == "no_crc": self.usable = pd.read_csv("the_big_usableno_crc.tsv", sep="\t", header=0, index_col=0)
        #elif modality == "only_corrected": self.usable = pd.read_csv("the_big_usableonly_corrected.tsv", sep="\t", header=0, index_col=0) 
        #elif no_crc == "only_multiple": = self.
        #else:
        self.usable = pd.read_csv("ml_dis_rf/usable_all_species_with_condition.tsv", sep="\t", header=0, index_col=0, low_memory=False, engine="c")
  
        #if modality == "only_multiple":
        #    self.usable = self.usable.loc[:, self.usable.loc["dataset_name"].isin(self.datasets)]
            
        self.feats = self.feats_frame.index.tolist()[:20]
 
        mean_pos = lambda feat : np.mean(self.usable.loc[feat, self.usable.loc["condition"]=="positive"].values.astype(float))
        mean_neg = lambda feat : np.mean(self.usable.loc[feat, self.usable.loc["condition"]=="negative"].values.astype(float))

        self.sign = dict([(feat, (-1 if (mean_neg(feat) > mean_pos(feat)) else 1)) for feat in self.feats])

	#        if modality == "only_multiple":
        #    self.P, self.N, self.dirP, self.dirN, self.A, self.A2 = self.prepare_data(data2aucs, data2aucs_one_at_the_time, True)
        #else:
        self.P, self.N, self.dirP, self.dirN, self.wilc_pos, self.wilc_neg, \
            self.Auc_direct, self.Auc_multiple, self.Auc_oneatthe, self.Auc_corrected, self.Auc_cross_val = \
	    self.prepare_data(data_aucs_direct, data_aucs_multiple_only, data_aucs_one_at_the_time, data_aucs_cross, data_corrected)
 
        #if modality == "no_crc": self.title = "healthy_microbiome_random_forest_modelization_no_crc"
        #elif modality == "only_corrected": self.title = "healthy_microbiome_random_forest_modelization_with_metadata"
        #elif modality == "only_multiple": self.title = "healthy_microbiome_random_forest_modelization_only_multiple"
        #else: self.title = "healthy_microbiome_random_forest_modelization"
        self.title = "combined_machine_learning_analysis_for_the_unhealthy_microbiome_RF"


    def prepare_data(self, data_aucs_direct, data_aucs_multiple_only, data_aucs_one_at_the_time, data_aucs_cross, data_corrected):
 
        #if multiple:
        #self.multiple_datasets += ["HMP_2019_ibdmdb_CD", "NielsenHB_2014_CD"]
        #self.multiple_names += ["CD", "CD"]

        #self.usable.loc["dataset_name"] = [( dat if (subt!="CD") else dat+"_CD"  ) for dat,subt in zip( \
        #        self.usable.loc["dataset_name"].tolist(), self.usable.loc["disease_subtype"].tolist() )]

        def get_abundance(feat, dataset):
            ab = np.mean(( self.usable.loc[feat, self.usable.loc["study_name"]==dataset].values.astype(float))) ## **2.)*100.)
            return ab if ab<=1.0 else np.log2(ab)

        def get_prevalence(feat, dataset):
            u = self.usable.loc[feat, self.usable.loc["study_name"]==dataset]
            prev = ( np.count_nonzero(u.values.astype(float))/float(len(u)) )  * 100.
  
            direction_n = np.mean(self.usable.loc[feat, ((self.usable.loc["study_name"]==dataset) & \
                (self.usable.loc["condition"]=="negative"))].values.astype(float))
            direction_p = np.mean(self.usable.loc[feat, ((self.usable.loc["study_name"]==dataset) & \
                (self.usable.loc["condition"]=="positive"))].values.astype(float))

            color = "blue" if (direction_p>direction_n) else "red"
            return prev, color


        ## this is the modification nicola suggested
        def get_prevalence_OR(feat, dataset):
            #u = self.usable.loc[feat, self.usable.loc["dataset_name"]==dataset].values.astype(float)

            print(dataset, feat, "  => This is the dataset & feat")
            usable_this_feat = self.usable.copy() #loc[:, :]
            #print(usable_this_feat) 

            if (dataset.startswith("HMP_2019_ibdmdb") or dataset.startswith("NielsenHB_2014")):

                if dataset.endswith("_CD"):
                    usable_this_feat = usable_this_feat.loc[:, usable_this_feat.loc["disease_subtype"]!="UC"]
                elif dataset.endswith("_UC"):
                    usable_this_feat = usable_this_feat.loc[:, usable_this_feat.loc["disease_subtype"]!="CD"]

            elif dataset.startswith("XieH_2016"):

                if dataset.endswith("_migraine"):
                    usable_this_feat = usable_this_feat.loc[:, usable_this_feat.loc["study_condition"]!="asthma"]
                elif dataset.endswith("_asthma"):
                    usable_this_feat = usable_this_feat.loc[:, usable_this_feat.loc["study_condition"]!="migraine"]

            else:
                usable_this_feat = usable_this_feat.loc[:, usable_this_feat.loc["study_name"]==dataset]
            
            array_positive = usable_this_feat.loc[feat, usable_this_feat.loc["condition"]=="positive"].values.astype(float)
            array_negative = usable_this_feat.loc[feat, usable_this_feat.loc["condition"]=="negative"].values.astype(float)
             
            ##array_positive = array_positive)**2.)*100.
            #array_negative = array_negative)**2.)*100.

            print(array_positive, " QUello dei positivi")
            print(array_negative, " Quello dei negativi")

            prevalence_positive = float(np.count_nonzero(array_positive)/float(len(array_positive)))
            prevalence_negative = float(np.count_nonzero(array_negative)/float(len(array_negative)))
 
            #print("p", prevalence_positive)
            #print("n", prevalence_negative)
 
            odd_p = (prevalence_positive * (1. - prevalence_positive)) #if (not float(prevalence_positive)==1.) else 1.0
            odd_n = (prevalence_negative * (1. - prevalence_negative)) #if (not float(prevalence_negative)==1.) else 1.0
 
            if float(odd_n) == 0.:
                OR = odd_p
            elif float() == 0.:
                OR = odd_p
            else:
                OR = odd_p / odd_n

            LOR = np.log10(OR) if (OR>0.0) else 0.0
            #print(LOR, " sono il LOR di %s " %feat)

            smd = (np.mean( array_positive ) - np.mean(array_negative))/\
                float(np.sqrt((np.std(array_positive, ddof=1)**2 + np.std( array_negative, ddof=1)**2) / 2.0))

            color = "blue" if (float(smd)>0.) else "red"
            wilcox = sts.ranksums(array_positive, array_negative)[1]

            #print("POR %s = " %feat, LOR, color)
            return LOR,color,("" if (wilcox>0.05) else ("+" if ( smd>0.0 ) else "-"))

        pos_feats = [feat for feat in self.feats if(self.sign[feat]>0.)]
        neg_feats = [feat for feat in self.feats if(self.sign[feat]<0.)]

        suppl_tab = {"Dataset": self.datasets}
        sign_tab = {"Dataset": self.datasets}

        for feat in pos_feats+neg_feats:
            prev_or = [ get_prevalence_OR(feat, dt) for dt in self.datasets ]

            lor = [ ch[0] for ch in prev_or ]
            #wilcs = [("" if (ch[0]>0.05) else ("+" if ((ch[0]<0.) else "-"))) for ch in prev_or]

            suppl_tab.update({feat: lor})
            
        if FEATID=="s__":
            suppl_tab = pd.DataFrame(suppl_tab).to_csv("Supplementary_Table_Machine_Features.tsv", sep="\t", header=True, index=False)
        else:
            suppl_tab = pd.DataFrame(suppl_tab).to_csv("Supplementary_Table_Machine_Features_KO.tsv", sep="\t", header=True, index=False)

        datasp = dict([(dat, [get_prevalence_OR(feat, dat) for feat in pos_feats]) for dat in self.datasets ])
        datasn = dict([(dat, [get_prevalence_OR(feat, dat) for feat in neg_feats]) for dat in self.datasets ])

        for k in datasp:
            print(k, datasp[k])

        positive = pd.DataFrame(dict([(k,[vv[0] for vv in v]) for k,v in datasp.items()]), index=pos_feats).T
        negative = pd.DataFrame(dict([(k,[vv[0] for vv in v]) for k,v in datasn.items()]), index=neg_feats).T

        positive_annot = pd.DataFrame(dict([(k, [vv[2] for vv in v]) for k,v in datasp.items()]), index=pos_feats).T
        negative_annot = pd.DataFrame(dict([(k, [vv[2] for vv in v]) for k,v in datasn.items()]), index=neg_feats).T

        # data_aucs_direct, data_aucs_multiple_only, data_aucs_one_at_the_time, data_aucs_cross, data_corrected

        Auc_direct = pd.DataFrame(dict([(dat, data_aucs_direct[dat]) for dat in self.datasets]), index=["Auc"]).T
        Auc_multiple = pd.DataFrame(dict([(dat, data_aucs_multiple_only[dat]) for dat in self.multiple_datasets]), index=["Auc"]).T
        Auc_oneatthe = pd.DataFrame(dict([(dat, data_aucs_one_at_the_time[dat]) for dat in self.multiple_datasets]), index=["Auc"]).T
        Auc_corrected = pd.DataFrame(dict([(dat, data_corrected[dat]) for dat in self.datasets if \
            (not dat in ["LiJ_2017", "LoombaR_2017"])]), index=["Auc"]).T
        Auc_cross_val = pd.DataFrame(dict([(dat, data_aucs_cross[dat]) for dat in self.datasets]), index=["Auc"]).T #+ \
            #["HMP_2019_ibdmdb_CD", "NielsenHB_2014_CD", "XieH_2016_ashtma"]]), index=["Auc"]).T

        dirp = pd.DataFrame(dict([(k,[vv[1] for vv in v]) for k,v in datasp.items()]), index=pos_feats).T
        dirn = pd.DataFrame(dict([(k,[vv[1] for vv in v]) for k,v in datasn.items()]), index=neg_feats).T

        print(Auc_direct.shape, " one shape")
        print(Auc_multiple.shape, " two shape")
        print(Auc_oneatthe.shape, " three shape")
        print(Auc_corrected.shape, " four shape")
        print(Auc_cross_val.shape, " five shape")

        return positive, negative, dirp, dirn, positive_annot, negative_annot, Auc_direct, Auc_multiple, Auc_oneatthe, Auc_corrected, Auc_cross_val



    def figure(self):
        fig = plt.figure(figsize=(12,24))
 
        self.datasets = ["XieH_2016_asthma", "XieH_2016_migraine", "ZhuF_2020", "JieZ_2017", "QinN_2014", \
            "RubelMA_2020", "YeZ_2018", "NagySzakalD_2017", "HMP_2019_ibdmdb_UC", "NielsenHB_2014_UC", "HMP_2019_ibdmdb_CD", "NielsenHB_2014_CD", \
            "QinJ_2012", "KarlssonFH_2013", "SankaranarayananK_2015", "ZellerG_2014", \
            "YuJ_2015", "FengQ_2015", "VogtmannE_2016", "ThomasAM_2019_a", "ThomasAM_2019_b", \
            "WirbelJ_2018", "GuptaA_2019", "HanniganGD_2017", "YachidaS_2019"] #
            ### "NielsenHB_2014_CD", "HMP_2019_ibdmdb_CD", "XieH_2016_ashtma"]
       
        #self.Auc_direct = self.Auc_direct.loc[self.datasets]

        self.Auc_multiple = self.Auc_multiple.loc[self.datasets[8:]]
        self.Auc_oneatthe = self.Auc_oneatthe.loc[self.datasets[8:]]
        #self.Auc_corrected = #self.Auc_corrected.loc[self.datasets[:-3]]
        #self.Auc_cross_val = self.Auc_cross_val.loc[self.datasets] #[:-3]]

        gs = gridspec.GridSpec(len(self.datasets) + 2, 6) #, width_ratios=[3, 3, 3, 3, 3, 3, 12, 6])

        ax_heatmap_auc_direct = plt.subplot(gs[ :-2, 1 ])
        ax_heatmap_auc_correct = plt.subplot(gs[ :-2, 2 ])

        ax_heatmap_auc_multiple = plt.subplot(gs[ 8:-2, 3 ])
        ax_heatmap_auc_oneatthe = plt.subplot(gs[ 8:-2, 4 ]) 
        ax_heatmap_auc_cross = plt.subplot(gs[ :-2, 5 ])
        #ax_heatmap_pos = plt.subplot(gs[ :-2, 6 ])
        #ax_heatmap_neg = plt.subplot(gs[ :-2, 7 ])
        #ax_colorbar_abundance_pos = plt.subplot(gs[ -1, 3:7 ])
        #ax_colorbar_abundance_neg = plt.subplot(gs[ -2, 3:7 ])
        ax_colorbar_auc = plt.subplot(gs[ 9:-6, 0 ])

        ### DEFINITIONS
        # self.P, self.N, self.dirP, self.dirN, self.Auc_direct, self.Auc_multiple, self.Auc_oneatthe, self.Auc_corrected, self.Auc_cross_val

        max_odd_value = 2.0

        cbar = matplotlib.colorbar.ColorbarBase(ax_colorbar_auc, cmap=matplotlib.cm.get_cmap(self.colormap), \
            norm=matplotlib.colors.Normalize(0.5, 0.95), orientation="vertical", extend="both")
 
        #cbar_ab = matplotlib.colorbar.ColorbarBase(ax_colorbar_abundance, cmap=matplotlib.cm.get_cmap(self.ab_colormap), \
        #    norm=matplotlib.colors.Normalize(0.0, np.log2(100.)), orientation="vertical", \
        #    ticks=[0.0, 0.5, 1.0, np.log2(3.0), np.log2(5.), np.log2(10.), np.log2(25), np.log2(50), np.log2(100.)])
 
        #cbar_ab_p = matplotlib.colorbar.ColorbarBase(ax_colorbar_abundance_pos, cmap=matplotlib.cm.get_cmap("PuBu_r"), \
        #    norm=matplotlib.colors.Normalize(-max_odd_value, 0.), orientation="horizontal", extend="min") #

        #cbar_ab_n = matplotlib.colorbar.ColorbarBase(ax_colorbar_abundance_neg, cmap=matplotlib.cm.get_cmap("Oranges"), \
        #    norm=matplotlib.colors.Normalize(0.0, max_odd_value), orientation="horizontal", extend="max")

            #ticks=[0.0, 0.5, 1.0, np.log2(3.0), np.log2(5.), np.log2(10.), np.log2(25), np.log2(50), np.log2(100.)])

        #cbar_ab.ax.set_yticklabels([( str(int(2.**float(x + 0.0001))) if float(x)>1.0 else str(float(x))[1:]) \
        #    for x in cbar_ab.ax.get_yticks()])

        ### self.Auc_direct, self.Auc_multiple, self.Auc_oneatthe, self.Auc_corrected, self.Auc_cross_val

        hm_direct = sns.heatmap(data=self.Auc_direct, annot=True, vmin=0.5, vmax=.95, ax=ax_heatmap_auc_direct, \
		cbar=False, square=False, cmap=self.colormap)
        hm_mul = sns.heatmap(data=self.Auc_multiple, annot=True, vmin=0.5, vmax=.95, ax=ax_heatmap_auc_multiple, \
		cbar=False, square=False, cmap=self.colormap, yticklabels=[])
                #["" for e in range(len(self.datasets)-2)] + self.datasets[-2:])
        hm_oneatthe = sns.heatmap(data=self.Auc_oneatthe, annot=True, vmin=0.5, vmax=.95, ax=ax_heatmap_auc_oneatthe, \
		cbar=False, square=False, cmap=self.colormap, yticklabels=[])
        hm_corrected = sns.heatmap(data=self.Auc_corrected, annot=True, vmin=0.5, vmax=.95, ax=ax_heatmap_auc_correct, \
		cbar=False, square=False, cmap=self.colormap)#, yticklabels=[])
        hm_cross = sns.heatmap(data=self.Auc_cross_val, annot=True, vmin=0.5, vmax=.95, ax=ax_heatmap_auc_cross, \
		cbar=False, square=False, cmap=self.colormap)#, yticklabels=[])

        #heatmap_auc = sns.heatmap(data=self.A, annot=True, vmin=0.5, vmax=1.0, ax=ax_heatmap_auc, cbar=False, square=False, cmap=self.colormap)
        #heatmap_auc2 = sns.heatmap(data=self.A2, annot=True, vmin=0.5, vmax=1.0, ax=ax_heatmap_auc2, cbar=False, square=False, \
        #    yticklabels=[], cmap=self.colormap)
 
        #mask_ = np.empty(self.N.shape)
        #for i,d in enumerate(self.N.index.tolist()):
        #    for j,feat in enumerate(self.N.columns.tolist()):
        #        if self.dirN.loc[d,feat] == "red": mask_[i,j] = True
        #        else: mask_[i,j] = False

        #print(mask_)

        ## self.wilc_pos, self.wilc_neg

        #heatmap_health = sns.heatmap(data=self.N*(-1.), vmin=0.0, vmax=(max_odd_value), ax=ax_heatmap_neg, cbar=False, square=False, cmap="Oranges", \
        #     yticklabels=[], mask=mask_, linewidths=1.5, linecolor="goldenrod", annot=self.wilc_neg, fmt="s")
   
        #heatmap_health2 = sns.heatmap(data=self.N*(-1.), vmin=(0), vmax=(max_odd_value), ax=ax_heatmap_neg, cbar=False, square=False, cmap="PuBu", \
        #    yticklabels=[], mask=np.logical_not(mask_), linewidths=1.5, linecolor="goldenrod", annot=self.wilc_neg, fmt="s")

        #mask_2 = np.empty(self.P.shape)
        #for i,d in enumerate(self.P.index.tolist()):
        #    for j,feat in enumerate(self.P.columns.tolist()):
        #        if self.dirP.loc[d,feat] == "blue": mask_2[i,j] = True
        #        else: mask_2[i,j] = False

        #heatmap_dysb = sns.heatmap(data=self.P*(-1.), vmin=0.0, vmax=(max_odd_value), ax=ax_heatmap_pos, cbar=False, square=False, cmap="PuBu", \
        #    yticklabels=[], mask=mask_2, linewidths=1.5, linecolor="goldenrod", annot=self.wilc_pos, fmt="s")

        #heatmap_dysb2 = sns.heatmap(data=self.P*(-1.), vmin=(0.), vmax=(max_odd_value), ax=ax_heatmap_pos, cbar=False, square=False, cmap="Oranges", \
        #    yticklabels=[], mask=np.logical_not(mask_2), linewidths=1.5, linecolor="goldenrod", annot=self.wilc_pos, fmt="s")

        #ax_heatmap_pos.xaxis.set_ticks_position("top") #, rotation=38)
        #ax_heatmap_neg.xaxis.set_ticks_position("top") #, rotation=38)
        #plt.setp(ax_heatmap_pos.get_xticklabels(), rotation=38, ha="left")
        #plt.setp(ax_heatmap_neg.get_xticklabels(), rotation=38, ha="left")
 
        [plt.savefig("../images/CLR_%s.%s"%(self.title,fmt), dpi=300) for fmt in ["svg", "png"]]


if __name__ == "__main__":
    hmp = healthy_microbiome_plot()
    hmp.figure()
    
