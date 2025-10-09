#!/usr/bin/env python

import pandas as pd
import numpy as np
from scipy import stats as sts
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
import seaborn as sns
import sys, os
import statsmodels.api as sm
import statsmodels.formula.api as smf
sys.path.append("../../python_modules/")
from meta_analyses import RE_meta_binary, RE_meta, singleStudyEffect, paule_mandel_tau
from statsmodels.stats.multitest import  fdrcorrection
import itertools as it
matplotlib.rcParams["svg.fonttype"] = "none"
from sklearn.metrics import roc_auc_score
from pathlib import Path

RANGE_ = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5] # 0.64, 0.128, 0.256, 0.512]
#RANGE = list(np.arange(0.0, .9, 0.05))


def meta_analysis(data, study_identifier, segregated, datasets, orals_, prevalence):
    study2lens = {}
    singleStudiesClass = []
    meta_an = {}
    print("************************************  META-ANALYSIS OF CLASSIFICATION AT PREVALENCE " + \
        (("%.3f" if isinstance(prevalence, float) else "%s") %prevalence))

    for dat in segregated: ## datasets:
        segregated[dat].loc["Oral_Richness"] = [ np.sum(segregated[dat].loc[[i for i in orals_ if i in segregated[dat].index], sam].values.astype(float)) \
	    for sam in segregated[dat].columns ]
        print("Estimating data %s" %dat, end="   ==> ")
        cohenD, p_value_cor, Lens, fake_p, ci, SEd, P, N, T = std_mean_diff(segregated[dat])
        singleStudiesClass += [singleStudyEffect((cohenD, p_value_cor), dat, Lens, False)]
        study2lens[dat] = {}
        study2lens[dat]["cases"] = Lens[1]
        study2lens[dat]["control"] = Lens[0]
        study2lens[dat]["partialP"] = fake_p if fake_p==fake_p else 1.0
        study2lens[dat]["ci"] = ci
        study2lens[dat]["SE"] = SEd
        study2lens[dat]["positive_percent"] = P
        study2lens[dat]["negative_percent"] = N
        sys.stdout.write("%.4f   (%.3f)\n" %(cohenD, p_value_cor))
        meta_an[dat + "_effect"] = cohenD
        meta_an[dat + "_pvalue"] = p_value_cor
        meta_an[dat + "_positive_perc"] = P
        meta_an[dat + "_negative_perc"] = N
        meta_an[dat + "_total_perc"] = T
    considerable = [e for e in singleStudiesClass if (e.accepted)]
    if len(considerable) >= 2:
        re = RE_meta_binary(\
          [e.effect for e in considerable], \
          [e.Pvalue for e in considerable], \
          [e.Name for e in considerable], \
          [study2lens[e.Name]["cases"] for e in considerable], \
          [study2lens[e.Name]["control"] for e in considerable], "Oral_Richness", "D", False, False, "PM") 
        sys.stdout.write("%s (%i studies) Random Effect = %.3f (%.3f) \n" %("Oral_Richness", len(considerable), re.RE, re.Pval))
        print("FIRST; OF THE RANDOM EFFECT")
        print(re.RE, "95% CI [", re.RE - (1.96 * re.stdErr), re.RE + (1.96 * re.stdErr), "]")
        meta_an["RE_effect"] = re.RE
        meta_an["RE_pvalue"] = re.Pval
        meta_an["RE_CI"] = "%.4f;%.4f" %(re.RE - (1.96 * re.stdErr), re.RE + (1.96 * re.stdErr))
        return pd.DataFrame(meta_an, index=[prevalence])


def std_mean_diff(data, problem="presence_of_disease", positive="positive", negative="negative", index="Oral_Richness", no_covariates=False):
    datat = data.T
    datat.fillna("NA", inplace=True)
    datat[index] = datat[index].values.astype(float)
    datat["age"] = datat["age"].values.astype(float)
    datat["BMI"] = datat["BMI"].values.astype(float)
    datat["number_reads"] = datat["number_reads"].values.astype(float)
    covariates = ["gender", "age", "BMI", "number_reads"]
    datat = datat[[index, problem] + covariates]

    #datat[ index ] = np.log( datat[index].values.astype(float) + 0.0000001 ) #
    P = np.count_nonzero(datat.loc[ datat[problem]=="positive", index ].values.astype(float))  
    N = np.count_nonzero(datat.loc[ datat[problem]=="negative", index ].values.astype(float)) 
   
    ## datat[ index ] = datat[index].values.astype(float)
    i = datat[ index ].values.astype(float).copy()
    datat[ index ] = np.log( i + np.min( i[i != 0] ) )

    #datat[ index ][ datat[ index ] == 0.0 ] = 0.0000001
    ## datat[ index ] = sts.zscore(datat[ index ].values.astype(float))
 
    #i = datat[ index ].values.astype(float)
    #i[ i == 0 ] = 0.000000001
    #i[ i != 0 ] = np.log( i[ i != 0 ] )
    #datat[ index ] = i

    #i = datat[ index ].values.astype(float).copy()
    #datat[ index ] = np.where(  i != 0, np.log(i), 0  )

    if not no_covariates: formula = ("%s ~ " %index) + " + ".join([problem] + covariates)
    else: formula = ("%s ~ %s" %(index, problem))

    md = smf.ols( formula, data=datat )
    model_fit = md.fit()
    t = model_fit.tvalues.loc["%s[T.%s]" %(problem, positive)]
    n1 = float(len(datat.loc[(datat[problem]==negative)]))
    n2 = float(len(datat.loc[(datat[problem]==positive)]))
    wald = model_fit.pvalues.loc["%s[T.%s]" %(problem, positive)]
    d = (t*(n1+n2))/float(np.sqrt(n1*n2)*np.sqrt(n1+n2-2))
    SEd = np.sqrt(((n1+n2-1)/float(n1+n2-3)) * ((4./float(n1+n2))*(1+((d**2.)/8.))))
    d_lw = d-(sts.t.ppf(0.975, n1+n2-1) * SEd)
    d_up = d+(sts.t.ppf(0.975, n1+n2-1) * SEd)
    #return d, wald, (n1, n2), 0.0, "N", SEd, (P/n2)*100, (N/n1)*100, ((P+N)/(n1+n2))*100.
    return model_fit.params.loc["%s[T.%s]" %(problem, positive)], wald, (n1, n2), \
 	0.0, "N", model_fit.bse.loc["%s[T.%s]" %(problem, positive)], (P/n2)*100, (N/n1)*100, ((P+N)/(n1+n2))*100.



def segregate_datasets(data): #, study_id = "study_identifier"):
    print("Initiating dataset segregation...")
    single_datasets = []
    dts = [ \
            "XieH_2016", "ZhuF_2020", "JieZ_2017", "QinN_2014", "RubelMA_2020", "YeZ_2018", "NagySzakalD_2017", "MetaCardis_2020_a", \
            "BedarfJR_2017", "HMP_2019_ibdmdb", "NielsenHB_2014", "HeQ_2017", "QinJ_2012", "KarlssonFH_2013", "SankaranarayananK_2015", \
            "ZellerG_2014", "YuJ_2015", "FengQ_2015", "VogtmannE_2016", "ThomasAM_2019_a", "ThomasAM_2019_b", "WirbelJ_2018", "GuptaA_2019", \
            "HanniganGD_2017", "YachidaS_2019" \
          ]

    conds = [ \
	    "asthma", "migraine", "schizofrenia", "ACVD", "cirrhosis", "STH", "BD", "ME/CFS", "CAD", "HF", "PD", \
            "IBD", "T2D", "CRC", "control" \
	  ]

    for dataset in data.loc["study_name"].unique():
        if (not dataset in ["XieH_2016", "HMP_2019_ibdmdb", "NielsenHB_2014", "MetaCardis_2020_a"]):
            this_estimate = data.loc[:, data.loc["study_name"]==dataset]
            key = dataset
            single_datasets += [(key, this_estimate)]
        else:
            if dataset in ["HMP_2019_ibdmdb", "NielsenHB_2014"]:
                for dis,anti_dis in [("UC", "CD"), ("CD", "UC")]:
                    this_estimate = data.loc[:, data.loc["study_name"]==dataset]
                    this_estimate = this_estimate.loc[:, this_estimate.loc["disease_subtype"]!=anti_dis]
                    key = dataset + "_" + dis
                    single_datasets += [(key, this_estimate)]
            elif dataset == "XieH_2016":
                for dis,anti_dis in [("asthma", "migraine"), ("migraine", "asthma")]:
                    this_estimate = data.loc[:, data.loc["study_name"]==dataset]
                    this_estimate = this_estimate.loc[:, this_estimate.loc["study_condition"]!=anti_dis]
                    key = dataset + "_" + dis
                    single_datasets += [(key, this_estimate)]
            elif dataset == "MetaCardis_2020_a":
                for dis,anti_dis in [("T2D", ("HF", "CAD")),  ("HF", ("T2D", "CAD")) , ("CAD", ("T2D", "HF"))]:
                    this_estimate = data.loc[:, data.loc["study_name"]==dataset]
                    this_estimate = this_estimate.loc[:, ~this_estimate.loc["study_condition"].isin(anti_dis)]
                    key = dataset + "_" + dis
                    single_datasets += [(key, this_estimate)]
    print("Terminated dataset aggregation")
    return dict(single_datasets)


def oral_species():
    cmd = pd.read_csv("../relative_abundances/oral_metagenomics.tsv", sep="\t", header=0, index_col=0, low_memory=False)    

    #cmd = cmd.loc[:, cmd.loc["age_category"]!="newborn"]
    #cmd = cmd.loc[:, cmd.loc["study_name"].isin(["FerrettiP_2018", "HMP_2012", "BritoIL_2016"])]
    #cmd = cmd.loc[:, ~((cmd.loc["study_name"]=="HMP_2012") & (cmd.loc["body_subsite"]!="tongue_dorsum"))]

    print("Available taxa and samples for oral samples in cMD3: ", cmd.shape)
    feats = [i for i in cmd.index if ("k__" in i)]
    abuns = dict([(ft, np.mean(cmd.loc[ft].values.astype(float))) for ft in feats])
    prevs = dict([(ft, np.mean([ (np.count_nonzero(cmd.loc[ft, cmd.loc["study_name"]==dt].values.astype(float))\
        /float(cmd.loc[:, cmd.loc["study_name"]==dt].shape[1])) for dt in cmd.loc["study_name"].unique().tolist() ]))  for ft in feats])
    return [(spp, prevs[spp]) for spp in feats], [(spp, abuns[spp]) for spp in feats]

#def coupled_samples(spps, abuns):
#    coupled = pd.read_csv("coupled_samples_stool_oral.tsv", sep="\t", header=0, index_col=0, low_memory=False)
#    coupled = coupled.loc[:, coupled.loc["study_name"].isin(["FerrettiP_2018", "HMP_2012", "BritoIL_2016"])]
#    coupled = coupled.loc[:, coupled.loc["age_category"]!="newborn"]
#    coupled = coupled.loc[:, coupled.loc["body_site"].isin(["stool", "oralcavity"])]
#    coupled = coupled.loc[:, ~((coupled.loc["study_name"]=="HMP_2012") & (coupled.loc["body_subsite"]!="tongue_dorsum"))]
#    print("Initializing oral-score estimation with a coupled dataset of: ", coupled.shape)
#    species = [i for i in coupled.index.tolist() if ("s__" in i)]
#    print("Before filter you have: %i" %len(species)) 
#    sam_2_body = dict([(s,b) for s,b in zip(coupled.columns.tolist(), coupled.loc["body_site"].tolist())])
#    subs = dict([(sub, coupled.loc[:, coupled.loc["subject_id"]==sub]) for sub in coupled.loc["subject_id"].unique().tolist()])
#    N = len(subs)
#    print("Tu hai un totale di ", N, " soggetti")
#    spps = dict(spps)
#    abuns = dict(abuns)
#    if not os.path.isfile("Oral_enrichment_supplementary_table.tsv"):
#        o = open("Oral_enrichment_supplementary_table.tsv", "w")
#        header = "\t".join(["Spp", "Prevalence", "N.stool", "Prev.stool", "N.oral", "Prev.oral\n"])
#        o.write(header)
#        for spp in spps:
#            n_stool = 0.
#            n_oral = 0.
#            for sub in subs:
#                sams_stool = np.count_nonzero([(float(coupled.loc[spp, sam]) \
#		    if spp in coupled.index else 0.0) for sam in subs[sub] if sam_2_body[sam]=="stool"])
#                sams_oral = np.count_nonzero([(float(coupled.loc[spp, sam]) \
#		    if spp in coupled.index else 0.0) for sam in subs[sub] if sam_2_body[sam]=="oralcavity"])
#                if sams_stool>0.:
#                    n_stool += 1.
#                if sams_oral>0:
#                    n_oral += 1.
#            line = "\t".join([spp, str(spps[spp]), str(n_stool), str(n_stool/N), str(n_oral), str(n_oral/N)])
#            o.write(line)
#        o.close()
#    table = pd.read_csv("Oral_enrichment_supplementary_table.tsv", sep="\t", header=0, index_col=0, low_memory=False)
#    return table


def write_all_lists(dataset_healthy, species_list_dir):
    Spps = [i for i in dataset_healthy.index.tolist() if i.startswith("k__")]
    preva_health = dict([(spp, np.mean([  (np.count_nonzero(dataset_healthy.loc[spp, dataset_healthy.loc["study_identifier"]==d].values.astype(float))/\
	float(len(dataset_healthy.loc[spp, dataset_healthy.loc["study_identifier"]==d]))) for d in dataset_healthy.loc["study_identifier"].unique().tolist()  ])) \
	for spp in Spps])

    prevs, abuns = oral_species()
    prevs = dict(prevs)
    
    outdir = Path(species_list_dir)
    outdir.mkdir(exist_ok=True)

    #table = coupled_samples(prevs, abuns)  ## this writes
    #for p in np.arange(0.0, 1.0, 0.05):
    #for p in RANGE_: #[0., .2, .4, .8, .16, .32, .64]:
    #    #table_here = table.loc[(table["Prev.oral"]>=p)]  #  & (table["Prev.stool"]<0.1)]
    #    species_list_at_p = [ k for k in prevs if ((prevs[k]>=p) and (preva_health.get(k, 0.0) < 0.2)) ] #table_here.index.tolist()
    #    print(p, "  oral species that meet conditions at this threshold ==> ", len(species_list_at_p))
    #    with open(os.path.join(species_list_dir, "species_list_at_%.2f.tsv" %p), "w") as osl:
    #        for spp in species_list_at_p: osl.write(spp + "\n")
    
    for p in RANGE_:   # Make sure RANGE_ is defined and an iterable
        species_list_at_p = [k for k in prevs if ((prevs[k] >= p) and (preva_health.get(k, 0.0) < 0.2))]
        print(p, "oral species that meet conditions at this threshold ==> ", len(species_list_at_p))
        output_file = outdir / ("species_list_at_%.2f.tsv" % p)
        with output_file.open("w") as osl:
            for spp in species_list_at_p:
                osl.write(spp + "\n")

def get_list_(p, species_list_dir):
    input_file = Path(species_list_dir) / ("species_list_at_%.2f.tsv" % p)
    with input_file.open() as osl:
        return [o.rstrip() for o in osl]
    
#def get_list_(p, species_list_dir):
#    with open(os.path.join(species_list_dir, "species_list_at_%.2f.tsv" %p)) as osl:
#        return [o.rstrip() for o in osl.readlines()]


def do_the_figure(dataset_healthy, species_list_dir, all_combined, dataset_sick, datasets_diseases, segregated, study_id="study_identifier"):
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(32, 6))
    gs = gridspec.GridSpec(1,4, width_ratios=[2,3,3,2])
    ax_nspp = plt.subplot(gs[0,0])
    ax_heal = plt.subplot(gs[0,1:3])
    ax_mean = plt.subplot(gs[0,3])
    #ax_sick = plt.subplot(gs[0,2])
    
    all_lists = dict([(p, [l.split("|")[-1] for l in get_list_(p, species_list_dir)]) for p in RANGE_])

    RANGE = list(RANGE_)

    for p in RANGE:
        print(p, " ====> ", len(all_lists[p]))

        #exit(1)

        frame_scatter_species = pd.DataFrame({\
            "Prevalence_threshold": [p for p in RANGE], \
        "N.species": [len(all_lists[p]) for p in RANGE]})

        frame_healthy = []
        study_unique = sorted(dataset_healthy.loc[study_id].unique().tolist())
        samples = dataset_healthy.columns.tolist()
        studies = dataset_healthy.loc[study_id].tolist()
        sam_2_std = dict([(sam,std) for sam,std in zip(samples, studies)])

        frame_sick = []
        study_unique_sick = sorted(dataset_sick.loc[study_id].unique().tolist())
        samples_sick = dataset_sick.columns.tolist()
        studies_sick = dataset_sick.loc[study_id].tolist()
        sam_2_std_sick = dict([(sam,std) for sam,std in zip(samples_sick, studies_sick)])

        for p in RANGE:
            L = [l for l in all_lists[p] if l in dataset_healthy.index]
            print(p, " ====> ", len(L))
            if len(L):
                frame_this_one = pd.DataFrame({ \
                    "Prevalence_threshold": [p for t in range(len(study_unique))], \
                    "Cohort": study_unique, \
                    "Cumulative": [(np.sum(dataset_healthy.loc[L, [sam for sam in samples if sam_2_std[sam]==std]].values.astype(float))/ \
                len([sam for sam in samples if sam_2_std[sam]==std])) for std in study_unique]}, \
                        index=[p for t in range(len(study_unique))])
                if not len(frame_healthy):
                    frame_healthy = frame_this_one
                else:
                    frame_healthy = frame_healthy.append(frame_this_one)
                print(p, "  ==> ", frame_healthy.shape)

            LL = [l for l in all_lists[p] if l in dataset_sick.index]
            if len(LL):
                frame_this_one_sick = pd.DataFrame({ \
                    "Prevalence_threshold": [p for t in range(len(study_unique_sick))], \
                    "Cohort": study_unique_sick, \
                    "Cumulative": [(np.sum(dataset_sick.loc[LL, [sam for sam in samples_sick if sam_2_std_sick[sam]==std]].values.astype(float))/ \
                        len([sam for sam in samples_sick if sam_2_std_sick[sam]==std])) for std in study_unique_sick]}, \
                        index=[p for t in range(len(study_unique_sick))])

                if not len(frame_sick):
                    frame_sick = frame_this_one_sick
                else:
                    frame_sick = frame_sick.append(frame_this_one_sick)
                print(p, "  ==> ", frame_sick.shape)
            
        #frame_healthy.to_csv("this_is_the_healthy_distribution.tsv", sep="\t", header=True, index=True)

        frame_scatter_species.loc[7, "Prevalence_threshold"] = 1.0 #.loc[frame_scatter_species.loc]
        #all_lists[]
        scatter_uno = sns.scatterplot(data=frame_scatter_species.loc[frame_scatter_species.index.tolist()], x="N.species", \
        y="Prevalence_threshold", ax=ax_nspp) #, palette=[("red" if (("%.2f" %p)[-1]=="0") else "gold") for p in RANGE])
      
        #ax_nspp.set_yticks([(x*100.) for x in RANGE][::-1])
        ax_nspp.set_xlim([0, 300])
        
        outdir=Path("meta_analyses_output")
        outdir.mkdir(exist_ok = True)

        if not os.path.isfile("meta_analyses_output/meta_analyses_at_all_thresholds.tsv"):
        #if not os.path.isfile("meta_analyses_at_all_thresholds.tsv"):
            ma = []
            for p in RANGE: #[0., .2, .4, .8, .16, .32, .64]:
                L = [l for l in all_lists[p] if l in all_combined.index]
                if len(L):
                    mma = meta_analysis(all_combined, study_id, segregated, datasets_diseases, L, p)
                    if not len(ma): ma = mma
                    else: ma = ma.append(mma)
            #print(ma)
            ma.fillna("NA").to_csv("meta_analyses_output/meta_analyses_at_all_thresholds.tsv", sep="\t", header=True, index=True)
            #ma.fillna("NA").to_csv("meta_analyses_at_all_thresholds.tsv", sep="\t", header=True, index=True)
        else:
            ma = pd.read_csv("meta_analyses_output/meta_analyses_at_all_thresholds.tsv", sep="\t", header=0, index_col=0).fillna("NA")
            #ma = pd.read_csv("meta_analyses_at_all_thresholds.tsv", sep="\t", header=0, index_col=0).fillna("NA")

        ma.index = [(("%.2f" if isinstance(i, float) else "%s") %i) for i in ma.index]
        ma["Prevalence_threshold"] = ma.index.tolist()
        frame_healthy["Cumulative"] = frame_healthy["Cumulative"].values.astype(float)
        frame_healthy["Prevalence_threshold"] = frame_healthy["Prevalence_threshold"].values.astype(str)

        frame_sick["Cumulative"] = frame_sick["Cumulative"].values.astype(float)
        frame_sick["Prevalence_threshold"] = frame_sick["Prevalence_threshold"].values.astype(str)

        frame_healthy["Diagnosis"] = "control"
        frame_sick["Diagnosis"] = "disease"
        frame_healthy = frame_healthy.append(frame_sick)
        
        scatter_due = sns.boxplot(data=frame_healthy.loc[frame_healthy.index.tolist()], \
        y="Cumulative", x="Prevalence_threshold", color="white", ax=ax_heal, hue="Diagnosis")
        #scatter_str = sns.stripplot(data=frame_healthy.loc[frame_healthy.index.tolist()], dodge=True, \
        #    y="Cumulative", x="Prevalence_threshold", palette={"control": "steelblue", "disease": "brown"}, ax=ax_heal, hue="Diagnosis")

        import pingouin as pg
        # initialize list with output statistics
        records = []
        for th in frame_healthy["Prevalence_threshold"].unique():
            fr = frame_healthy.loc[frame_healthy["Prevalence_threshold"]==th]
            print(fr)
            print(fr.columns)
            #exit()
            eff_disease_vs_control = pg.compute_effsize(
                fr.loc[fr["Diagnosis"] == "disease", "Cumulative"],
                fr.loc[fr["Diagnosis"] == "control", "Cumulative"]
            )
            
            eff_control_vs_disease = pg.compute_effsize(
                fr.loc[fr["Diagnosis"] == "control", "Cumulative"],
                fr.loc[fr["Diagnosis"] == "disease", "Cumulative"]
            )
        
            mw_test = sts.mannwhitneyu(
                fr.loc[fr["Diagnosis"] == "disease", "Cumulative"],
                fr.loc[fr["Diagnosis"] == "control", "Cumulative"]
            )

            auc = roc_auc_score(
                [1.0] * len(fr.loc[fr["Diagnosis"] == "disease", "Cumulative"]) + 
                [0.0] * len(fr.loc[fr["Diagnosis"] == "control", "Cumulative"]),
                fr.loc[fr["Diagnosis"] == "disease", "Cumulative"].tolist() + 
                fr.loc[fr["Diagnosis"] == "control", "Cumulative"].tolist()
            )

            records.append({
                "Prevalence_threshold": th,
                "Effect_size_disease_vs_control": eff_disease_vs_control,
                "Effect_size_control_vs_disease": eff_control_vs_disease,
                "Mann_Whitney_statistic": mw_test.statistic,
                "Mann_Whitney_pvalue": mw_test.pvalue,
                "AUC": auc
            })

            print("THRESHOLD: ", th, " EFFECT: ", pg.compute_effsize( fr.loc[fr["Diagnosis"]=="disease", "Cumulative"], fr.loc[fr["Diagnosis"]=="control", "Cumulative"]), \
                  "MEAN D.: ", pg.compute_effsize(fr.loc[fr["Diagnosis"]=="control", "Cumulative"], fr.loc[fr["Diagnosis"]=="disease", "Cumulative"]), \
                  "PVALUE: ", sts.mannwhitneyu(fr.loc[fr["Diagnosis"]=="disease", "Cumulative"], fr.loc[fr["Diagnosis"]=="control", "Cumulative"]), \
                  "AUC: ", \
                     roc_auc_score( \
                     [1. for i in range(len(fr.loc[fr["Diagnosis"]=="disease", "Cumulative"]))] + \
                     [0.0 for i in range(len(fr.loc[fr["Diagnosis"]=="control", "Cumulative"].tolist()))], \
                     fr.loc[fr["Diagnosis"]=="disease", "Cumulative"].tolist() + fr.loc[fr["Diagnosis"]=="control", "Cumulative"].tolist() \
                     ))
        
        results_df = pd.DataFrame(records)
        
    
    ax_heal.set_yscale("log")
    #ax_sick.set_xscale("log")
    ax_nspp.set_xticks(np.arange(0.0, 550, 50))
    pops = [c.replace("_effect", "") for c in ma.columns if (c.endswith("_effect") and (not c.startswith("RE_")))]

    print(ma)    
    meta_an_data = pd.DataFrame({ \
	"Population": list(it.chain.from_iterable([ pops for i in RANGE ])),
	"Prevalence_threshold": list(it.chain.from_iterable([ [  (("%.2f" if isinstance(i, float) else "%s") %i) for j in range(len(pops))] for i in RANGE  ])), 
	"Effects": list(it.chain.from_iterable([ [float(ma.loc[  (("%.2f" if isinstance(i, float) else "%s") %i), pop+"_effect"]) for pop in pops] for i in RANGE ])), 
	"Pvalues": list(it.chain.from_iterable([ [float(ma.loc[  (("%.2f" if isinstance(i, float) else "%s") %i), pop+"_pvalue"]) for pop in pops] for i in RANGE ])), 
        "Positive": list(it.chain.from_iterable([ [float(ma.loc[ (("%.2f" if isinstance(i, float) else "%s") %i), pop+"_positive_perc"]) for pop in pops] for i in RANGE ])),
        "Negative": list(it.chain.from_iterable([ [float(ma.loc[ (("%.2f"if isinstance(i, float) else "%s")  %i), pop+"_negative_perc"]) for pop in pops] for i in RANGE ])), 
        "Total": list(it.chain.from_iterable([ [float(ma.loc[    (("%.2f" if isinstance(i, float) else "%s") %i), pop+"_total_perc"]) for pop in pops] for i in RANGE ])) #_total_perc
	}, index=list(it.chain.from_iterable([ pops for i in RANGE ])))
    
    #meta_an_data["Percentage"]

    mapper_disease = dict([
        ("XieH_2016_asthma", "asthma"), \
	("XieH_2016_migraine", "migraine"), \
	("ZhuF_2020", "schizofrenia"), \
	("MetaCardis_2020_a_CAD", "CAD"), \
        ("JieZ_2017", "ACVD"), \
	("QinN_2014", "cirrhosis"), \
	("RubelMA_2020", "STH"), \
	("YeZ_2018", "BD"), \
	("NagySzakalD_2017", "ME/CFS"), \
	("MetaCardis_2020_a_T2D", "T2D"), \
	("MetaCardis_2020_a_HF", "HF"), \
        ("BedarfJR_2017", "PD"), \
	("HMP_2019_ibdmdb_CD", "CD"), \
	("NielsenHB_2014_CD", "CD"), \
	("HMP_2019_ibdmdb_UC", "UC"), \
        ("NielsenHB_2014_UC", "UC"), \
	("HeQ_2017", "CD"), \
	("QinJ_2012", "T2D"), \
	("KarlssonFH_2013", "T2D"), \
	("SankaranarayananK_2015", "T2D"), \
        ("ZellerG_2014", "CRC"), \
	("YuJ_2015", "CRC"), \
	("FengQ_2015", "CRC"), \
	("VogtmannE_2016", "CRC"), \
	("ThomasAM_2019_a", "CRC"), \
	("ThomasAM_2019_b", "CRC"), \
	("WirbelJ_2018", "CRC"), \
	("GuptaA_2019", "CRC"), \
        ("HanniganGD_2017", "CRC"), \
	("YachidaS_2019", "CRC") \
      ])

    meta_an_data["Disease"] = [mapper_disease[ch] for ch in meta_an_data["Population"]]
    meta_an_data["Significance"] = [("yes" if p<0.05 else "no") for p in meta_an_data["Pvalues"].values]
    #scatter_man = sns.stripplot(data=meta_an_data, x="Total", marker="s", y="Prevalence_threshold", ax=ax_mean, color="red", alpha=0.7)
    scatter_man2 = sns.boxplot(data=meta_an_data, x="Total", y="Prevalence_threshold", ax=ax_mean, color="steelblue", boxprops=dict(alpha=.3))
    plt.tight_layout()
    [plt.savefig("../images/oral_enrichment_validation.%s" %fmt, dpi=300) for fmt in ["svg", "png"]]
    


def main():
    species_list_dir = "oral_species_by_prevalence"
    dataset_healthy = pd.read_csv("../relative_abundances/ALL_adult_controls.tsv", sep="\t", header=0, index_col=0, low_memory=False)
    dataset_sick = pd.read_csv("../relative_abundances/ALL_adult_sicks.tsv", sep="\t", header=0, index_col=0, low_memory=False)

    dataset_disease = pd.read_csv("../relative_abundances/usable_all_species_with_condition.tsv", sep="\t", header=0, index_col=0, low_memory=False) 

    write_all_lists(dataset_healthy, species_list_dir)
    
    dataset_healthy.index = [ (i if (not "s__" in i) else i.split("|")[-1]) for i in dataset_healthy.index ]
    dataset_sick.index    = [ (i if (not "s__" in i) else i.split("|")[-1]) for i in dataset_sick.index ]
    dataset_disease.index = [ (i if (not "s__" in i) else i.split("|")[-1]) for i in dataset_disease.index ]

    segregated = segregate_datasets(dataset_disease)
    dts = [ \
        "XieH_2016", "ZhuF_2020", "JieZ_2017", "QinN_2014", "RubelMA_2020", "YeZ_2018", "NagySzakalD_2017", "MetaCardis_2020_a", \
        "BedarfJR_2017", "HMP_2019_ibdmdb", "NielsenHB_2014", "HeQ_2017", "QinJ_2012", "KarlssonFH_2013", "SankaranarayananK_2015", \
        "ZellerG_2014", "YuJ_2015", "FengQ_2015", "VogtmannE_2016", "ThomasAM_2019_a", "ThomasAM_2019_b", "WirbelJ_2018", "GuptaA_2019", \
        "HanniganGD_2017", "YachidaS_2019" ]
    conds = [ \
        "asthma", "migraine", "schizofrenia", "ACVD", "cirrhosis", "STH", "BD", "ME/CFS", "CAD", "HF", "PD", \
        "IBD", "T2D", "CRC", "control" ]
    do_the_figure(dataset_healthy, species_list_dir, dataset_disease, dataset_sick, dts, segregated)

if __name__ == "__main__":
    main()
