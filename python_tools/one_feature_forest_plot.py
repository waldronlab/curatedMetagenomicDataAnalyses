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
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
matplotlib.rcParams["svg.fonttype"] = "none"
from scipy import stats as sts
sns.set_style("whitegrid")
import argparse as ap

def read_params():
    p = ap.ArgumentParser() 
    add = p.add_argument
    add("analysis", type=str, help="The table with the meta-analys[i/e]s generated")
    add("outcome", type=str, help="The outcome on which the meta-analysis is done (must be an row-name in 'analysis')")
    add("-ld", "--lenght_data", type=str, default=None, \
	help="The data on which the analysis has been comnputed can be used to estimate CI." \
	"If not speficid, normal distribution will be used")
    add("-sid", "--study_id", type=str, default="study_name", help="Used when data are specific for confidence intervals")    
    add("-of", "--outpref", type=str, help="Out prefix: facultative")
    add("-ess", "--effect_suff", type=str, default="_Effect")
    add("-qs", "--q_suff", type=str, default="_Qvalue")
    add("-a", "--significance", type=float, default=0.05)
    add("-tl", "--title", type=str, default=None)
    add("-pd", "--pos_dir", type=str, default="Positive")
    add("-nd", "--neg_dir", type=str, default="Negative")
    add("-mx", "--max_rho", type=float, default=None)
    add("-mn", "--min_rho", type=float, default=None)
    return p.parse_args()

def confIntMeanT(mn, se, n, conf=0.95):
    m = sts.t.ppf((1+conf)/2., n-1)
    return mn - m*se, mn + m*se

def confIntMeanN(mn, se, conf=0.95):
    return sts.norm.interval(0.95, loc=mn, scale=se)

def forest_plot(data, outfile, prevotella, significance, positive, negative, title, min_rho, max_rho):
    fs = 32
    fig = plt.figure(figsize=(10, 10)) 
    gs = gridspec.GridSpec(2,2,height_ratios=[len(data), 0.2])
    ax = plt.subplot(gs[0,:])

    ax_arrows_p = plt.subplot(gs[1,1])
    ax_arrows_n = plt.subplot(gs[1,0])
 
    #ax_arrows_p.set_xlim([0,ax.get_xlim()[1]])
   # ax_arrows_n.set_xlim([ax.get_xlim()[0], 0])
 
    #if min_rho and max_rho:
    #    ax.set_xlim([-min_rho, max_rho])
    #    ax_arrows_p.set_xlim([-min_rho, 0.])
    #    ax_arrows_n.set_xlim([0., max_rho])

    ax2 = ax.twinx()
    styles = {"pop": "s", "ma": "D"}
    data["Style"] = (["pop"] * (len(data) -1)) + ["ma"]

    a = sns.scatterplot(\
        data=data, x="Effect", y="Study", ax=ax, hue="Significance", \
        palette={"significant": "cornflowerblue", "non-significant": "black"}, \
        style="Style", markers=styles, size="Weight", sizes=(40, 800))

    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(ax.get_yticks()[::-1])

    ax_arrows_p.set_xlim([0,ax.get_xlim()[1]])
    ax_arrows_n.set_xlim([ax.get_xlim()[0], 0])

    for std,i,lw,up,Q in zip( \
        data["Study"].tolist(), \
        data["Idx"].tolist(), \
        data["CI_l"].tolist(), \
        data["CI_u"].tolist(), \
        data["Q"].tolist()):
        ax.plot([lw, up], [i, i], c="black" if Q>significance else "cornflowerblue", linewidth=2.2)

    ax.axvline(0.0, linestyle="--", color="black")
    fmt_tck = lambda l,u,w : "%.2f [%.3f, %.3f]" %(w, l, u) if isinstance(w, float) else "[%.3f, %.3f]" %(l, u)
    yticks = []
    for std,i,lw,up,w in zip(\
        data["Study"].tolist(), data["Idx"].tolist(), data["CI_l"].tolist(), data["CI_u"].tolist(), (data["Weight"].tolist()[:-1] + ["MA"])):
        yticks += [fmt_tck(lw,up,w)]

    a.get_legend().remove()
    ax2.set_yticklabels(yticks[::-1])
    ax_arrows_p.set_xlim([ax.get_xlim()[0], 0.])
    ax_arrows_n.set_xlim([0., ax.get_xlim()[1]])

    leg_pos = ax_arrows_p.arrow(x=0.02, y=0.5, dx=ax.get_xlim()[1]*0.60, dy=0.0, color="cornflowerblue", \
        shape="full", width=.075, length_includes_head=True, linestyle="-", \
        head_width=0.5, head_length=0.07)
    leg_neg = ax_arrows_n.arrow(x=-0.02, y=0.5, dx=ax.get_xlim()[0]*0.60, dy=0.0, color="cornflowerblue", \
        shape="full", width=.075, length_includes_head=True, linestyle="-", \
        head_width=0.5, head_length=0.07)

    ax_arrows_p.set_title(positive)
    ax_arrows_n.set_title(negative) 
    ax_arrows_p.axis("off")
    ax_arrows_n.axis("off")

    if title:
        ax.set_title(title)

    plt.tight_layout()
    [plt.savefig("%s.%s"%(outfile, fmt), dpi=200) for fmt in ["svg", "png"]]
    plt.close()


def get_lengths(data, study_id):
    a = pd.read_csv(data, sep="\t", header=0, index_col=0, low_memory=False, engine="c")
    d = dict([(s, float(a.loc[:, a.loc[study_id]==s].shape[1])) for s in a.loc[study_id].unique()])
    return d

def main( args ):
    data = pd.read_csv(args.analysis, sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna(np.nan)

    if not args.outcome in data.index:
        raise SyntaxError("%s not in the row-names of the analysis! exiting..." %args.outcome)
        exit(1)
 
    data = data.loc[args.outcome]
    studies = [ c.replace( args.effect_suff, "") for c in data.index if (c.endswith( args.effect_suff ) and not c.startswith("RE")) ] + [ "Random-effects model" ]
    effects = [ (s + args.effect_suff) for s in studies[:-1] ] + [ "RE" + args.effect_suff ]
    ses = [ (s+"_SE") for s in studies[:-1] ] + [ "RE_stdErr" ]
    qs = [ (s + args.q_suff) for s in studies[:-1] ] + [ "RE" + args.effect_suff + args.q_suff ]

    if args.lenght_data:
        lenghts = get_lengths(args.lenght_data, args.study_id)
        CI = [confIntMeanT( mn,se,l ) for mn,se,l in zip(\
	    data.loc[effects[:-1]].values.astype(float), data.loc[ses[:-1]].values.astype(float), [lenghts[s] for s in studies[:-1]] )]
    else:
        CI = [confIntMeanN( mn,se ) for mn,se in zip(\
            data.loc[effects[:-1]].values.astype(float), data.loc[ses[:-1]].values.astype(float))]

    print(data)
    print(studies)
    print(effects)
    print(ses)
    print(qs)
 
    rCI = data.loc["RE_conf_int"].split(";")

    ci_lows = [x[0] for x in CI] + [float(rCI[0])]
    ci_upps = [x[1] for x in CI] + [float(rCI[1])]

    plottable_data = pd.DataFrame({ \
        "Effect": data.loc[effects].values.astype(float), "CI_l": ci_lows, "CI_u": ci_upps, \
        "Weight": [(((1/x)/np.sum([(1/float(y)) for y in (data.loc[ses].values[:-1].astype(float)**2.)]))*100.) for x in (data.loc[ses].values[:-1].astype(float)**2.)] + [50.], \
	"Study": studies, "Q": data.loc[qs].values.astype(float), "Idx": [i for i in range(len(studies))] }, index=studies)

    plottable_data[ "Significance" ] = [ ("significant" if q<args.significance else "non-significant") \
	for q in data.loc[qs].values.astype(float) ]

    outfile = args.outpref + "_" + args.outcome
    print(plottable_data)
 
    forest_plot( plottable_data, outfile, args.outcome, args.significance, args.pos_dir, args.neg_dir, args.title, args.min_rho, args.max_rho )

if __name__ == "__main__":
    main(read_params())
