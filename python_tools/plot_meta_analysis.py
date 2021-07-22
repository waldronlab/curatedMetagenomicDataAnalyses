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

def readargs(args):
    arp = ap.ArgumentParser()
    add = arp.add_argument
    add("meta_analysis", type=str)
    add("-u", "--usabledata", type=str, default="MODS/usable_multisub_14age.tsv")
    add("-re", "--random_effect", type=str, default="RE_SpearmanR")
    add("--min_rho", type=float, default=0.05)
    add("--neg_max_rho", type=float, default=0.30)
    add("--pos_max_rho", type=float, default=0.50)
    add("--top", type=int, default=20)
    add("--random_effect_FDR", type=str, default="RE_Qvalue") ## FDR_Qvalue  ## FDR_Qvalue
    add("--single_data_suffs", type=str, default="_SpearmanR") # _FDR-P-value
    add("--single_data_FDR_suffs", type=str, default="_Qvalue")
 
    #add("--conf_") # "--single_data_FDR_suffs ".FDR_Qvalue"

    add("--qvalue", type=float, default=0.05)
    add("--dotsize", type=int, default=9)
    add("--diamsize", type=int, default=17)
    add("--outfile", type=str, default="meta_analysis_image")
    add("-x", "--effectname", type=str, default="Spearman Correlation Coefficient")
    add("-lp", "--linepositive", type=str, default="elderly")
    add("-ln", "--linenegative", type=str, default="young")
    add("-ll", "--legloc", type=str, default="best")
    add("-sp", "--suptitle", type=str, default="RE meta-analysis")
    add("-z", "--feature_id", type=str, default="s__")
    add("-cb", "--color_blue", type=str, default="darkcyan")
    add("-cr", "--color_red", type=str, default="darkorange")
    add("-b", "--color_black", type=str, default="black")

    add("-sm", "--stack_metaanalysis", type=str, default="")
    add("-su", "--stack_usable", type=str, default="")

    add("-R", "--REC", action="store_true")
    add("-A", "--A", action="store_true")
    add("-KO", action="store_true")
    add("-PWY", action="store_true")
    add("-se", "--suppl_effect", type=str, default="", help="aname:afilename")
    add("--min_prev", type=float, default=0.01)
    return arp.parse_args() 
 
###self.args.feat

def dataset_cardinality(usable_filename):
    ub = pd.read_csv(usable_filename, sep="\t", header=0, index_col=0, low_memory=False)
    cardinalities = {}
    for d in ub.loc["study_name"].unique():
        cardinalities[d] = float(ub.loc[:, ub.loc["study_name"].isin([d])].shape[1])
    return cardinalities
    


def readinput(args):
    table = pd.read_csv(args.meta_analysis_res, sep="\t", header=0, index_col=0, low_memory=False)
    table.fillna(np.nan, inplace=True)
    table = table.astype(float, errors="ignore")


    if args.REC:
        disease_columns = [c for c in table.columns if c.endswith(args.single_data_suffs)]
        acceptable = lambda spc : [effect for effect in table.loc[spc, disease_columns].tolist() if (float(effect) == float(effect))]
        table = table.loc[[spc for spc in table.index.tolist() if len(acceptable(spc))>=4]]
 
    if (not args.REC) and args.KO:
        disease_columns = [c for c in table.columns if c.endswith(args.single_data_suffs)]
        acceptable = lambda spc : [effect for effect in table.loc[spc, disease_columns].tolist() if (float(effect) == float(effect))]
        table = table.loc[[spc for spc in table.index.tolist() if len(acceptable(spc))>=4]]
 
    if (not args.REC) and args.PWY:
        disease_columns = [c for c in table.columns if c.endswith(args.single_data_suffs)]
        acceptable = lambda spc : [effect for effect in table.loc[spc, disease_columns].tolist() if (float(effect) == float(effect))]
        table = table.loc[[spc for spc in table.index.tolist() if len(acceptable(spc))>=4]]

   
    usables = pd.read_csv(args.usabledata, sep="\t", header=0, index_col=0, low_memory=False).fillna("NA")
    if args.stack_usable:
        usables = usables.append(pd.read_csv(args.stack_usable, sep="\t", header=0, index_col=0, low_memory=False).fillna("NA"))
 
    features_to_exclude = [j for j in [i for i in usables.index.tolist() if \
	    (("s__" in i) or ("g__" in i) or ("K" in i) or ("PWY" in i))] if \
            (np.count_nonzero( usables.loc[j].values.astype(float))/float(usables.shape[1]))< args.min_prev ]
    
    table.drop([ft for ft in features_to_exclude if ft in table.index.tolist()], inplace=True)

    if args.stack_metaanalysis:
        table2 = pd.read_csv(args.stack_metaanalysis, sep="\t", header=0, index_col=0, low_memory=False)
        table2.fillna(np.nan, inplace=True)
        table2 = table2.astype(float, errors="ignore")
        table = table.append(table2)
        #if args.REC:
        #table.drop(features_to_exclude, inplace=True)
        table2.drop([ft for ft in features_to_exclude if ft in table2.index.tolist()], inplace=True, errors="ignore")
        genera = [x.replace("g__", "s__") for x in table2.index.tolist()]
        species = [y for y in table.index.tolist() if y.startswith("s__")]
        to_drop = []
        for gn in genera:
            #print(gn, gn.replace("s__", "g__") in table2.index, gn.replace("s__", "g__") in table.index)
            for spc in species:
                if ((spc.startswith(gn)) and (\
                    ("%.2f" %table.loc[spc, args.random_effect]) == \
                    ("%.2f" %table.loc[gn.replace("s__", "g__"), args.random_effect]))):
                    to_drop += [gn.replace("s__", "g__")]
        table.drop(to_drop, inplace=True)
        
        disease_columns = [c for c in table.columns if c.endswith(args.single_data_suffs)]
        acceptable = lambda spc : [effect for effect in table.loc[spc, disease_columns].tolist() if (float(effect) == float(effect))]
        table = table.loc[[spc for spc in table.index.tolist() if len(acceptable(spc))>=4]]
      
        features_to_exclude = [j for j in [i for i in usables.index.tolist() if \
        (("s__" in i) or ("g__" in i) or ("K" in i) or ("PWY" in i))] if \
        (np.count_nonzero( usables.loc[j].values.astype(float))/float(usables.shape[1]))< args.min_prev ]
        table.drop([ft for ft in features_to_exclude if ft in table.index.tolist()], inplace=True)
       
    #for x,y in zip(table[args.random_effect], table[args.random_effect_FDR]): print(x,y)

    if args.top<0:
        result_features = [i for i in table.index.tolist() if np.abs(table.loc[i, args.random_effect])>=args.min_rho \
            and (table.loc[i, args.random_effect_FDR]<args.qvalue)]  ##  and (np.abs(table.loc[i, args.random_effect]>=args.min_rho)))]
    else:
        table = table.loc[(table[args.random_effect_FDR]<args.qvalue), :]
        table["ABS"] = np.abs(table[args.random_effect].values)
        table.sort_values("ABS", ascending=False, inplace=True)

        result_features = table.index.tolist()[:args.top]

        #[table.loc[i, args.random_effect] ]
        #result_features = sorted([i for i in table.index.tolist() if (table.loc[i, args.random_effect_FDR]<args.qvalue)])


    positive_correlated = [f for f in result_features if (table.loc[f, args.random_effect]>0.)]
    negative_correlated = [f for f in result_features if (table.loc[f, args.random_effect]<0.)]
    sorted_table = table.loc[negative_correlated].sort_values(args.random_effect, ascending=True).append(\
        table.loc[positive_correlated].sort_values(args.random_effect, ascending=True))


    result_features = sorted_table.index.tolist()
    Dataset_Names = []

    for e,x in enumerate(sorted_table.columns.tolist()):
        if x.endswith(args.single_data_FDR_suffs) and (not x.startswith("RE_")):
            Dataset_Names += [x.replace(args.single_data_FDR_suffs, "")]
            sorted_table.insert(e+1, x.replace(args.single_data_FDR_suffs, "") + "_color", \
                [(args.color_red if (float(sorted_table.loc[i, x] if \
                str(sorted_table.loc[i, x])[0].isdigit() else 1000)<0.2) else args.color_blue) \
                for i in sorted_table.index.tolist()])

    return sorted_table, result_features, positive_correlated, negative_correlated, Dataset_Names



def plot_meta_analysis(args):
    sns.set_style("whitegrid")
    sorted_table, result_features, positive_correlated, negative_correlated, Dataset_Names = readinput(args)
    sorted_table.index = result_features

    cardinalities = dataset_cardinality(args.usabledata)

    REframe = pd.DataFrame({"Microbial Species": [(x[3:].replace("_"," ") if (args.feature_id=="s__") else x) for x in result_features], \
        args.effectname: [sorted_table.loc[feat, args.random_effect] for feat in result_features], \
        "RESG": [("yes" if (sorted_table.loc[x, args.random_effect_FDR]<args.qvalue) else "no") for x in result_features], \
        "Q-Value": [("%.4f" %sorted_table.loc[q, args.random_effect_FDR]) for q in result_features], \
        "lower": [float( sorted_table.loc[ x, "CI_RE" ].split(";")[0]) for x in result_features], \
        "upper": [float( sorted_table.loc[ x, "CI_RE" ].split(";")[1]) for x in result_features] })

    fig = plt.figure(figsize=(16 if (args.feature_id=="s__") else 22, 16))
    gs = gridspec.GridSpec(2,1,height_ratios=[1]+[len(result_features)])
    ax = plt.subplot(gs[1,0])
    ax_legend = plt.subplot(gs[0,0])

    FRAME = pd.DataFrame({\
            "Microbial Species": [ (x[3:].replace("_"," ") if (args.feature_id=="s__") else x) for x in \
            list(it.chain.from_iterable([[feat for t in range(len(Dataset_Names))] for feat in result_features]))],\
            args.effectname: \
            list(it.chain.from_iterable([[sorted_table.loc[feat, d+args.single_data_suffs] for d in Dataset_Names] \
            for feat in result_features])), \
            "Significance": \
            list(it.chain.from_iterable([[sorted_table.loc[feat, d+"_color"] for d in Dataset_Names] \
            for feat in result_features])), \
            "Dataset_Name": \
            list(it.chain.from_iterable([Dataset_Names for feat in result_features]))})

    FRAME["Num-Species"] = np.arange(0.0, len(FRAME), 1.0)

    FRAME.dropna(inplace=True)
 
    FRAME.to_csv("%s_dataframe.tsv" %args.outfile, sep="\t", header=True, index=False)

    strip_one = sns.scatterplot(\
            x=args.effectname, y="Microbial Species", data=FRAME, size="Significance", sizes={args.color_blue: args.dotsize*8, args.color_red: args.dotsize*16}, \
        palette={args.color_red: args.color_red, args.color_blue: args.color_blue}, edgecolor="black", ax=ax, hue="Significance") ## size="Cardinality"
 
    box_one = sns.boxplot(\
        x=args.effectname, y="Microbial Species", data=FRAME, saturation=0.25, \
        color="palegoldenrod", ax=ax, width=0.45)

    strip_RE = sns.stripplot(\
        x=args.effectname, y="Microbial Species", data=REframe, \
        marker="D", ax=ax, color=args.color_black, size=args.diamsize) #palette={"no": "goldenrod", "yes": "black"}, size=args.diamsize)

    REframe.index = REframe["Microbial Species"]

    for e,feat in enumerate(REframe["Microbial Species"].tolist()):
        lw, upp = float(REframe.loc[feat, "lower"]), float(REframe.loc[feat, "upper"])
        ax.plot([ lw, lw, lw, upp, upp, upp ], [ e-0.01, e+0.01, e, e, e+0.01, e-0.01 ], c=args.color_red, linewidth=3.6)

    ax.set_xlim([-args.neg_max_rho, args.pos_max_rho])
    #ax.set_yticklabels([item.get_text()[4:].replace("_"," ") for item in ax.get_xticklabels()])
    #ax.set_yticklabels([x[4:].replace("_"," ") for x in ax.get_yticks()])


    red_line = ax.axvline(ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], linewidth=3.5, color=args.color_blue, linestyle="--", alpha=1.0)
    cyan_linem2 = ax.axvline(x=-0.1, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], linewidth=2.0, color=args.color_red, linestyle="--", \
        alpha=1.0)
    cyan_linep2 = ax.axvline(x=0.1, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], linewidth=2.0, color=args.color_red, linestyle="--", \
        alpha=1.0)

    ax_legend.set_xlim([-args.neg_max_rho,args.pos_max_rho])
    ax_legend.set_ylim([-1.0,1.0])
    leg_pos = ax_legend.arrow(x=0.02, y=0.5, dx=args.pos_max_rho*0.80, dy=0.0, color=args.color_blue, \
        shape="full", width=.075, length_includes_head=True, linestyle="-", \
        head_width=0.5, head_length=0.07)
    leg_neg = ax_legend.arrow(x=-0.02, y=0.5, dx=-args.neg_max_rho*0.80, dy=0.0, color=args.color_blue, \
        shape="full", width=.075, length_includes_head=True, linestyle="-", \
        head_width=0.5, head_length=0.07)

    #zero = ax_legend.get_position()
    #one = [zero.x0, zero.y0-2.0, zero.height, zero.width]
    #ax_legend.set_position(one)
    ax_legend.annotate(args.linepositive, xy=(args.pos_max_rho*0.80 - 0.01*len(args.linepositive) - 0.015, -1.0), fontsize=16)
    ax_legend.annotate(args.linenegative, xy=(-args.neg_max_rho*0.80, -1.0), fontsize=16)

    ax_legend.axis("off")

    leg_handles = [(mlines.Line2D([], [], color=color, marker=marker, linestyle='None', alpha=1.0, \
        markersize=markersize, label=label)) for label,color,marker,markersize in zip(\
        ["Non-Significant", "Significant", "Random-Effect-Model"], \
        [args.color_blue, args.color_red, args.color_black], ["o", "o", "D"], [args.dotsize, args.dotsize, args.diamsize])]

    legend = ax.legend(handles=leg_handles, fontsize=11, loc=args.legloc, frameon=True, markerfirst=True, ncol=1)
    #ax_legend.axis("off")

    #print(list(ax.get_xticklabels()))
    #ax.set_xticklabels(ax.get_xticklabels(), fontsize=16)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=16)
    ax.set_ylabel(ax.get_ylabel(), fontsize=16)
    ax.set_xlabel(ax.get_xlabel(), fontsize=16)

    #for d in dir(box_one): print(d)
    #axQ = ax.twinyi()

    #print(axQ.get_yticks())
    #axQ.set_yticks(np.arange(0.5, 10., 1.225))#len(ax.get_yticks()+1)))#ax.get_yticks()[::-1])

    ax.tick_params(bottom=True, top=False, left=True, right=True)
    plt.subplots_adjust(left=0.32 if args.feature_id=="s__" else 0.48)
    for lab in ax.get_yticklabels():
        lab.set_style('italic')

    plt.suptitle(args.suptitle, fontsize=18)
    [plt.savefig("%s.%s" %(args.outfile, fmt), dpi=400) for fmt in ["png", "svg"]]
    return REframe




def plot_meta_analysis_of_meta_analysis(args, mma="diseases"):
    sns.set_style("whitegrid")
    sorted_table, result_features, positive_correlated, negative_correlated, Dataset_Names = readinput(args)
 
    #REframe = pd.DataFrame({"Microbial Species": \
    #    [(x[3:].replace("_"," ") if (not args.KO) else x) for x in result_features], \
    #    args.effectname: [sorted_table.loc[feat, args.random_effect] for feat in result_features], \
    #    "RESG": [("yes" if (sorted_table.loc[x, args.random_effect_FDR]<args.qvalue) else "no") \
    #	 for x in result_features], \
    #    "Q-Value": [("%.4f" %sorted_table.loc[q, args.random_effect_FDR]) for q in result_features]})

    REframe = pd.DataFrame({"Microbial Species": [(x[3:].replace("_"," ") if (args.feature_id=="s__") else x) for x in result_features], \
        args.effectname: [sorted_table.loc[feat, args.random_effect] for feat in result_features], \
        "RESG": [("yes" if (sorted_table.loc[x, args.random_effect_FDR]<args.qvalue) else "no") for x in result_features], \
        "Q-Value": [("%.4f" %sorted_table.loc[q, args.random_effect_FDR]) for q in result_features], \
        "lower": [float( sorted_table.loc[ x, "CI_RE" ].split(";")[0]) for x in result_features], \
        "upper": [float( sorted_table.loc[ x, "CI_RE" ].split(";")[1]) for x in result_features] })


    fig = plt.figure(figsize=(16 if (not args.KO) else 24,16))
    gs = gridspec.GridSpec(2,1,height_ratios=[1]+[len(result_features)])
    ax = plt.subplot(gs[1,0])
    ax_legend = plt.subplot(gs[0,0])

    if mma=="diseases":
        markers = ["X", "h", "p", "P", "*", "v", "8", "<", ">", "^", "s", "o"]
        colors = ["darkgrey"] * 12

    d2markers = dict([(n,m) for n,m in zip(Dataset_Names, markers)])
    d2colors = dict([(n,c) for n,c in zip(Dataset_Names, colors)])
 

    FRAME = pd.DataFrame({\
        "Microbial Species": [(x[3:].replace("_"," ") if (not args.KO) else x) for x in \
        list(it.chain.from_iterable([[feat for t in range(len(Dataset_Names))] for feat in result_features]))],\
        args.effectname: \
        list(it.chain.from_iterable([[sorted_table.loc[feat, d+args.single_data_suffs] \
	    for d in Dataset_Names] for feat in result_features])), \
        "FDR":  \
        list(it.chain.from_iterable([[sorted_table.loc[feat, d + args.single_data_FDR_suffs ] \
	    for d in Dataset_Names] for feat in result_features])), \
	"Significance": \
        list(it.chain.from_iterable([[sorted_table.loc[feat, d+"_color"] for d in Dataset_Names] \
	    for feat in result_features])), \
	"Markers":\
	list(it.chain.from_iterable([[d for d in Dataset_Names] \
	    for i in range(len(result_features))]))})

    FRAME["Num-Species"] = np.arange(0.0, len(FRAME), 1.0)

    FRAME.dropna(inplace=True)

    strip_one = sns.scatterplot(\
        x=args.effectname, y="Microbial Species", data=FRAME, s=args.diamsize*15, \
        palette={args.color_red: args.color_red, args.color_blue: args.color_blue}, markers=d2markers, \
	edgecolor="black", ax=ax, hue="Significance", style="Markers")

    box_one = sns.boxplot(\
        x=args.effectname, y="Microbial Species", data=FRAME, saturation=0.25, \
        color="palegoldenrod", ax=ax, width=0.45)

    strip_RE = sns.stripplot(\
        x=args.effectname, y="Microbial Species", data=REframe, \
        marker="D", ax=ax, color=args.color_black, size=args.diamsize) 

    REframe.index = REframe["Microbial Species"]

    for e,feat in enumerate(REframe["Microbial Species"].tolist()):
        lw, upp = float(REframe.loc[feat, "lower"]), float(REframe.loc[feat, "upper"])
        ax.plot([ lw, lw, lw, upp, upp, upp ], [ e-0.01, e+0.01, e, e, e+0.01, e-0.01 ], c=args.color_red, linewidth=3.6)

    ax.set_xlim([-args.neg_max_rho, args.pos_max_rho])

    red_line = ax.axvline(ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], \
	linewidth=3.5, color=args.color_blue, linestyle="--", alpha=1.0)
    cyan_linem2 = ax.axvline(x=-0.1, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], \
	linewidth=2.0, color=args.color_red, linestyle="--", alpha=1.0)
    cyan_linep2 = ax.axvline(x=0.1, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], \
	linewidth=2.0, color=args.color_red, linestyle="--", alpha=1.0)

    ax_legend.set_xlim([-args.neg_max_rho,args.pos_max_rho])
    ax_legend.set_ylim([-1.0,1.0])
    leg_pos = ax_legend.arrow(x=0.02, y=0.5, dx=args.pos_max_rho*0.80, dy=0.0, color=args.color_blue, \
        shape="full", width=.075, length_includes_head=True, linestyle="-", \
        head_width=0.5, head_length=0.07)
    leg_neg = ax_legend.arrow(x=-0.02, y=0.5, dx=-args.neg_max_rho*0.80, dy=0.0, color=args.color_blue, \
        shape="full", width=.075, length_includes_head=True, linestyle="-", \
        head_width=0.5, head_length=0.07)

    ax_legend.annotate(args.linepositive, \
	xy=(args.pos_max_rho*0.80 - 0.01*len(args.linepositive) - 0.015, -1.0), fontsize=16)
    ax_legend.annotate(args.linenegative, \
	xy=(-args.neg_max_rho*0.80, -1.0), fontsize=16)

    ax_legend.axis("off")
    leg_handles = [(mlines.Line2D([], [], color=color, marker=marker, linestyle='None', alpha=1.0, \
        markersize=args.dotsize, label=label)) for label,color,marker in zip(\
        Dataset_Names + ["Significant", "Non-significant", "Random effect"], \
	[d2colors[d] for d in Dataset_Names] + [args.color_red, args.color_blue, args.color_black], \
	[d2markers[d] for d in Dataset_Names] + ["o", "o", "D"])]

    legend = ax.legend(handles=leg_handles, fontsize=11, loc=args.legloc, \
	frameon=True, markerfirst=True, ncol=1)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=16)
    ax.set_ylabel(ax.get_ylabel(), fontsize=16)
    ax.set_xlabel(ax.get_xlabel(), fontsize=16)
    ax.tick_params(bottom=True, top=False, left=True, right=True)

    plt.subplots_adjust(left=(0.3 if (not args.KO) else 0.45))

    if (not args.KO):
        for lab in ax.get_yticklabels(): lab.set_style('italic')
    plt.suptitle(args.suptitle, fontsize=18)
    [plt.savefig("%s.%s" %(args.outfile, fmt), dpi=400) for fmt in ["png", "svg"]]
    return REframe





def abundances(frame, usab, usab2=False, Class="gender", class0="female", class1="male", title="sex", featid="s__", fancy_name="Biological sex", PREVALENCE=False):

    if(not isinstance(usab2, bool)):
        usab = usab.append(usab2.loc[[i for i in usab2.index if i.startswith("g__")]]) #, left_index=True, right_index=True, how="outer")
        usab.fillna("NA", inplace=True)

    males = usab.loc[:, usab.loc[Class].isin(["1.0"])]
    females = usab.loc[:, usab.loc[Class].isin(["0.0"])]

    classes = []
    values = []
    species = []
    pvalues = []
    species_level = []
   

    for en,microbe in enumerate(frame["Microbial Species"].tolist()):
        if featid!="PWY":
            Microbe = [x for x in usab.index.tolist() if x.endswith(microbe.replace(" ","_"))][0]
        else: ##if featid=="PWY":
            Microbe = microbe #[x for x in usab.index.tolist() if x.endswith(microbe)]


        if(not PREVALENCE):
            mean_in_males = [np.mean([y for y in list(((np.sin(males.loc[Microbe, \
                males.loc["study_name"].isin([d])].values.astype(float)))**2.0)*100.)]) for d in males.loc["study_name"].unique()]
            mean_in_females = [np.mean([x for x in list(((np.sin(females.loc[Microbe, \
                females.loc["study_name"].isin([d])].values.astype(float)))**2.0)*100)]) for d in females.loc["study_name"].unique()]


        else:

            mean_in_males = [(np.count_nonzero( males.loc[Microbe, males.loc["study_name"].isin([d])].values.astype(float))\
                /float( males.loc[Microbe, males.loc["study_name"].isin([d])].shape[0] ))*100. \
                for d in males.loc["study_name"].unique()]
            mean_in_females = [(np.count_nonzero( females.loc[Microbe, females.loc["study_name"].isin([d])].values.astype(float))\
                /float( females.loc[Microbe, females.loc["study_name"].isin([d])].shape[0] ))*100. \
                for d in females.loc["study_name"].unique()]


        classes += [class0 for i in range(len(mean_in_females))]
        classes += [class1 for i in range(len(mean_in_males))]
        values += mean_in_females #[mean_in_females ]
        values += mean_in_males #[mean_in_males]
        species += [microbe for i in range(len(mean_in_females))]
        species += [microbe for i in range(len(mean_in_males))]
        species_level += [en for e in range(len(mean_in_females))]
        species_level += [en for e in range(len(mean_in_males))] 

    frame = pd.DataFrame({fancy_name: classes, "Relative Abundance" if (not PREVALENCE) else "Prevalence": values, \
        "Microbial Species": species, "species_level": species_level})

    fig,ax = plt.subplots(figsize=[8 if featid=="s__" else 12 ,13])

    if(not PREVALENCE):

        boxes = sns.boxplot(x="Relative Abundance" if (not PREVALENCE) else "Prevalence", y="Microbial Species", \
            data=frame, hue=fancy_name, ax=ax, dodge=True, \
            palette={class0: "goldenrod", class1: "steelblue"}, boxprops=dict(alpha=.25))

        bars = sns.swarmplot(x="Relative Abundance" if (not PREVALENCE) else "Prevalence", \
            y="Microbial Species", data=frame, hue=fancy_name, ax=ax, dodge=True, \
            palette={class0: "goldenrod", class1: "steelblue"}) #, whis=3)A
 
        ax.set(xscale="log")

    else: 
        bars = sns.violinplot(x="Relative Abundance" if (not PREVALENCE) else "Prevalence", \
            y="Microbial Species", data=frame, hue=fancy_name, ax=ax, dodge=True, inner="stick", \
            palette={class0: "goldenrod", class1: "steelblue"}, scale="area", scale_hue=True, split=True)

    ax.set_title(("%s Bugs" %fancy_name) + ("Relative Abundances" if (not PREVALENCE) else "Prevalence"), fontsize=16)
    plt.subplots_adjust(left=0.32 if featid=="s__" else 0.45)

    [plt.savefig("%s_relative_abundances.%s" %(title, fmt), dpi=300) for fmt in ["svg", "png"]]    


if __name__ == "__main__":
    Args = readargs(sys.argv)
    if (not Args.REC):
        REFrame = plot_meta_analysis(Args)
        if Args.A:
            abundances(REFrame, pd.read_csv(Args.usabledata, sep="\t", header=0, index_col=0), False \
                if(not Args.stack_usable) else pd.read_csv(Args.stack_usable, sep="\t", header=0, index_col=0), \
                Class="gender", featid=Args.feature_id)
    else:
        REFrame = plot_meta_analysis_of_meta_analysis(Args)
