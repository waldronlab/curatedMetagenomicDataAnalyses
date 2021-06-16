#!/usr/bin/env python


import pandas as pd
import numpy as np
import os, argparse as ap
import sys
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
import itertools as it
from num2words import num2words



def read_params( args ):
    p = ap.ArgumentParser(description=\
        "", \
        formatter_class=ap.RawTextHelpFormatter)
    add = p.add_argument
 
    add("metaanalysis", type=str, nargs="+", help="A meta-analysis with features as rows and fields as columns")
 
    add("-n", "--narrowed", action="store_true", help="")

    add("-a", "--relative_abundances", type=str, nargs="+", help="A file with features as rows and samples as columns")

    add("-i", "--imp", type=int, default=30, help="Tries to plot the first *imp*ortant features. ")

    add("-how", "--how", type=str, default="first", choices=["first", "union", "intersection"], help=\
        "How to chose the [--imp] features to plot, default : intersection. NOTE: "
        "meaningfull only for multiple meta-analysis plot")

    add("-ps", "--positive_direction", type=str, default="Positive side", help=\
        "Label of the +1 group (with arrow)")

    add("-ns", "--negative_direction", type=str, default="Negative side", help=\
        "Label of the -1 group (with arrow)")

    add("-x", "--x_axis", type=str, default="Effect-size")

    add("-y", "--y_axis", type=str, default="")

    add("-es", "--e_suff", type=str, default=["_es"], nargs="+")
 
    add("-qs", "--q_suff", type=str, default=["_Q"], nargs="+")

    add("-pe", "--prevalence", type=str, default=[""], nargs="+")

    add("-mp", "--min_prevalence", type=float, default=[0.01], nargs="+")

    add("-mna", "--min_ab", type=float, default=[""], nargs="+")

    add("-ms", "--min_studies", type=int, nargs="+", default=[4], help=\
        "Minimum number of meta-analysis in which a feature must be present "
        "with a real effect-size to be considerable for the hierarchical meta-analysis.")

    add("-m", "--markers", action="store_true", help=\
        "Activating this flag, each single-analysis/population is "
        "is represented with a different marker/symbol.")

    add("-o", "--outfile", type=str, default="", help="Uses the input name of the first meta-analysis if not specified. ")
    
    add("-re", "--random_effect", type=str, nargs="+", default=["Effect-size"], help=\
        "Random/fixed effect model effect-size")

    add("-ci", "--confint", type=str, nargs="+", default=["Conf-int"], help=\
        "Random/fixed effect model confidence intervals")

    add("-rq", "--random_effect_q", type=str, nargs="+", default=["FDR-Q"], help=\
        "Random/fixed effect model Q-value")

    add("-cr", "--color_red", type=str, nargs="+", default=["goldenrod"], help=\
        "Color of the significant results and of the two important lines marking the desired thresholds")

    add("-cb", "--color_blue", type=str, nargs="+", default=["dodgerblue"], help=\
        "Color the non-significant results and of the zero-marking line ")

    add("-cbl", "--color_black", type=str, nargs="+", default=["black"], help=\
        "Color of the random/fixed effect diamond ")

    add("-dmr", "--diam_marker", type=str, nargs="+", default=["D"], help=\
        "Marker of the random/fixed effect [default=Diamond]")

    add("-il", "--important_lines", type=float, nargs="+", default=[0.3], help=\
        "Two lines are draw at this level to mark good results [default = 0.3]")

    add("-as", "--a_single", type=float, default=0.2, \
        help="Significant alpha for the single analyses/populations")

    add("-ar", "--a_random", type=float, default=0.05, \
        help="Significant alpha for the random/fixed effect model")

    add("--dotsize", type=int, default=9)

    add("--diamsize", type=int, default=17)

    add("--neg_max_rho", type=float, default=0.80)

    add("--pos_max_rho", type=float, default=0.80)

    add("--legloc", type=str, default="best")

    print(vars(p.parse_args()))

    return p.parse_args()





def parse_args_table(args, ij):

    name = num2words(ij, ordinal=True)

    take_a_param = lambda param, ij : param[ij] if len(param)>1 else param[0]

    if isinstance( args, dict ):
        from attrdict import AttrDict
        args = AttrDict(args)

    meta_analysis = pd.read_csv( take_a_param( args.metaanalysis, ij ), sep="\t", header=0, index_col=0, low_memory=False, engine="c")

    meta_analysis.fillna( "NA", inplace=True )

    if (not take_a_param( args.relative_abundances, ij)) and (not take_a_param( args.prevalence, ij )):
        raise Warning("You must specifiy one of prevalence OR relative_abundances (to now the minimum prevalence)")

    if take_a_param( args.relative_abundances, ij) and take_a_param( args.prevalence, ij ):
        filter_on = "PREVALENCE"
        raise Warning("You set both prevalence AND relative_abundances: the first one will be used")

    else:
        filter_on = "PREVALENCE" if take_a_param( args.prevalence, ij ) else "ABUNDANCE"

    min_prevalence = dict( )
    features = meta_analysis.index.tolist( )

    if take_a_param( args.relative_abundances, ij ):
        abundances = pd.read_csv( take_a_param( args.relative_abundances, ij), sep="\t", header=0, index_col=0, low_memory=False, engine="c")
        prevalences = dict([ (feat, (np.count_nonzero(abundances.loc[feat, :].values.astype(float))/float(abundances.shape[1]) \
            if feat in abundances.index.tolist() else 1.0) ) for feat in features ])
 
    if take_a_param( args.prevalence, ij ): # and not take_a_param( args.min_ab, ij )):  
        prevalences = dict([ (feat, meta_analysis.loc[feat, take_a_param( args.prevalence, ij )]) for feat in features])

    single_effects = [ c for c in meta_analysis.columns.tolist() if ((c.endswith( take_a_param( args.e_suff, ij ) )) and (c != take_a_param( args.random_effect, ij ) ))  ]
    single_effect_Qs = [ c for c in meta_analysis.columns.tolist() if ((c.endswith( take_a_param( args.q_suff, ij ) )) and (c != take_a_param( args.random_effect_q, ij ) ))  ]
 
    results = meta_analysis.loc[ meta_analysis[ take_a_param( args.random_effect_q, ij ) ]<\
            (args.a_random if args.how=="first" else 1.0), : ]
    results = results.loc[ [ i for i in results.index.tolist() if ( meta_analysis.loc[ i, single_effects ].tolist().count("NA")>= take_a_param( args.min_studies, ij) ) ], :]
    
    if filter_on == "PREVALENCE":
        results = results.loc[ [ i for i in results.index.tolist() if (prevalences[i]>=take_a_param( args.min_prevalence, ij))], : ]
    else:
        results = results.loc[ [ i for i in results.index.tolist() if \
            (not i in abundances.index.tolist()) or (np.mean(abundances.loc[i].values.astype(float))>=take_a_param( args.min_prevalence, ij)) \
            ], : ]
    
    if args.narrowed:
        results = results.loc[[i for i in results.index.tolist() if \
            (float(results.loc[i, take_a_param( args.confint, ij) ].split(";")[0])*float(results.loc[i, take_a_param( args.confint, ij)].split(";")[1])) \
            >0.0], :]

    results["abs"] = np.abs( results[ take_a_param( args.random_effect, ij) ].values.astype(float) )
    results.sort_values( "abs", inplace=True, ascending=False )

    if args.how=="first":
        results = results[:args.imp]

    positive = [i for i in results.index.tolist() if ( float(results.loc[i, take_a_param( args.random_effect, ij)] >0.) ) ]
    negative = [i for i in results.index.tolist() if ( float(results.loc[i, take_a_param( args.random_effect, ij)] <0.) ) ]
    res_pos = results.loc[ positive, : ].sort_values( take_a_param( args.random_effect, ij), ascending=False )
    res_neg = results.loc[ negative, : ].sort_values( take_a_param( args.random_effect, ij), ascending=True )
    results = res_pos.append(res_neg)

    return [results, prevalences, single_effects, single_effect_Qs, results.index.tolist( ), name]
            





def select_features(how, results):
    if how == "first":
        return results[0].index.tolist()
    elif how == "union":
        return set().union(*results)
    elif how == "intersection":
        return set().intersection(*results)
 




 
def build_long_frame(args, analysis, ij, result_features, all_markers):
 
    results, prevalences, single_effects, single_effect_Qs, Feats, name_of_this_one = tuple(analysis)
    take_a_param = lambda param, ij : param[ij] if len(param)>1 else param[0]

    print(name_of_this_one)
    print(Feats)
    print(results)

    lgf_singles = pd.DataFrame(\
        {
        "Feature": list( it.chain.from_iterable( [ [ feature for e in range(len(single_effects)) if \
            str(results.loc[feature, single_effects[e]])!="NA"] for feature in Feats ] ) ) , \
        "effect-size": list( it.chain.from_iterable( [ [ float(results.loc[ feature, e ])  for e in single_effects if \
            str(results.loc[feature, e])!="NA"] for feature in Feats ] ) ), \
        "population": list( it.chain.from_iterable( [ [ e.replace(take_a_param( args.e_suff, ij), "")  for e in single_effects if \
            str(results.loc[feature, e])!="NA"] for feature in Feats ] ) ),
        "q-values": list( it.chain.from_iterable( [ [ float(results.loc[ feature, e ])  for e in single_effect_Qs if \
            str(results.loc[feature, e])!="NA"] for feature in Feats ] ) ),
        } \
        )
    
    if args.markers:
        if isinstance(all_markers, int):
            all_markers = [".", ",", "v", "^", "<", ">", "1", "2", "3", "8", "s", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d", "|"]
        pop_to_marker = {}
        for pop in lgf_singles["population"].unique().tolist():
            if not pop in pop_to_marker:
                pop_to_marker[pop] = all_markers[-1]
                all_markers = all_markers[:-1]
 
        lgf_singles["marker_of_population"] = [pop_to_marker[pop] for pop in lgf_singles["population"].tolist()]

    remaining_markers = all_markers

    lgf_singles.index = lgf_singles["Feature"].tolist()
    lgf_singles["SG"] = [ ("yes" if (q<args.a_single) else "no") for q in lgf_singles["q-values"].tolist() ]
 
    confint = [ results.loc[ feature, take_a_param( args.confint, ij) ].split(";") for feature in Feats ]

    lgf_effect = pd.DataFrame(\
        { "Feature": Feats, \
          "random-effect": [ float(results.loc[ feature, take_a_param( args.random_effect, ij) ]) for feature in Feats ], \
          "Q-values": [ float(results.loc[ feature, take_a_param( args.random_effect_q, ij) ]) for feature in Feats ], \
          "CI_low":  [ float(ci[0]) for ci in confint ], \
          "CI_upp":  [ float(ci[1]) for ci in confint ], \
        } \
        )

    lgf_effect.index = lgf_effect["Feature"].tolist()
    lgf_effect["SG"] = [ ("yes" if (q<args.a_random) else "no") for q in lgf_effect["Q-values"].tolist() ]
    
    return lgf_singles, lgf_effect, name_of_this_one, remaining_markers





def draw_figure(args):

    analyses = [\
        parse_args_table( args, ij ) for ij in range(len( args.random_effect ))
        ]

    result_features = select_features( args.how, [ an[0] for an in analyses ] )
 
    take_a_param = lambda param, ij : param[ij] if len(param)>1 else param[0]
 
    labels, colors, markers = [], [], []

    # [results, prevalences, single_effects, single_effect_Qs, results.index.tolist( )]

    fig = plt.figure(figsize=(22, 16))

    gs = gridspec.GridSpec( 2,1, height_ratios=[1]+[len(result_features)] )

    ax = plt.subplot( gs[ 1,0 ] )

    ax_arrows = plt.subplot( gs[ 0,0 ])
  
    ####all_markers = [".", ",", "o", "v", "^", "<", ">", "1", "2", "3", "8", "s", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d", "|"]
 

    for ij,analysis in enumerate( analyses ):

        if not ij: 
            long_form_frame_singles, long_form_frame_effect, name_of, remaining_markers = build_long_frame( args, analysis, ij, result_features, -1 )
        else:
            long_form_frame_singles, long_form_frame_effect, name_of, remaining_markers = build_long_frame( args, analysis, ij, result_features, remaining_markers )
     

        if not args.markers:

            strip_one = sns.scatterplot(\
                x="effect-size", y="Feature", \
                data=long_form_frame_singles, size="SG", \
                sizes={"no": args.dotsize*8, "yes": args.dotsize*16}, \
                palette={"yes": take_a_param(args.color_red, ij), \
                "no": take_a_param(args.color_blue, ij)}, \
                edgecolor="black", ax=ax, hue="SG")

        else:

            marker_of_this_population = dict([(p,m)\
                for p,m in zip( long_form_frame_singles["population"].tolist(), long_form_frame_singles["marker_of_population"].tolist() )])

            for p,m in marker_of_population.items():
                labels += [ p ]
                markers += [ m ]
                markers += [ "darkgrey" ]

            strip_one = sns.scatterplot(\
                x="effect-size", y="Feature", \
                data=long_form_frame_singles, size="SG", \
                sizes={"no": args.dotsize*8, "yes": args.dotsize*16}, \
                palette={"yes": take_a_param(args.color_red, ij), \
                "no": take_a_param(args.color_blue, ij)}, \
                edgecolor="black", ax=ax, hue="SG", \
                style="style", markers=marker_of_population)

        strip_RE = sns.stripplot(\
            x="random-effect", y="Feature", data=long_form_frame_effect, \
            marker=take_a_param(args.diam_marker, ij), ax=ax, color=take_a_param(args.color_black, ij), \
            size=args.diamsize, palette={"no": "goldenrod", "yes": "black"}, \
            hue="SG")
 
        for e,feat in enumerate(long_form_frame_effect["Feature"].tolist()):
            lw, upp = float(long_form_frame_effect.loc[feat, "CI_low"]), float(long_form_frame_effect.loc[feat, "CI_upp"])
            ax.plot([ lw, lw, lw, upp, upp, upp ], [ e-0.01, e+0.01, e, e, e+0.01, e-0.01 ], \
            c=take_a_param(args.color_red, ij), linewidth=3.6)

        blue_line = ax.axvline(ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], linewidth=3.5, \
            color=args.color_blue[0], linestyle="--", alpha=1.0)

        red_line_1 = ax.axvline(x=-take_a_param(args.important_lines, ij), ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], \
            linewidth=2.0, color=take_a_param(args.color_red, ij), linestyle="--", alpha=1.0)

        red_line_2 = ax.axvline(x=take_a_param(args.important_lines, ij), ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], \
            linewidth=2.0, color=take_a_param(args.color_red, ij), linestyle="--", alpha=1.0)
  

        if len(long_form_frame_effect["SG"].unique().tolist())==2:
            labels += [ name_of, name_of ]
            markers += [ take_a_param(args.diam_marker, ij), take_a_param(args.diam_marker, ij) ]
            colors += [ take_a_param(args.color_black, ij), "yellow" ]

        else:
            labels += [ name_of ]
            markers += [ take_a_param(args.diam_marker, ij) ]
            colors += [ take_a_param(args.color_black, ij) ]
 
        labels += [ "Non significant", "Significant" ]
        markers += [ "o", "o" ]
        colors += take_a_param(args.color_blue, ij), take_a_param(args.color_red, ij)


    ax_arrows.set_xlim([-args.neg_max_rho, args.pos_max_rho])
    ax_arrows.set_ylim([-1.0,1.0])

    leg_pos = ax_arrows.arrow(x=0.02, y=0.5, dx=args.pos_max_rho*0.80, dy=0.0, color=args.color_blue[0], \
        shape="full", width=.075, length_includes_head=True, linestyle="-", \
        head_width=0.5, head_length=0.07)

    leg_neg = ax_arrows.arrow(x=-0.02, y=0.5, dx=-args.neg_max_rho*0.80, dy=0.0, color=args.color_blue[0], \
        shape="full", width=.075, length_includes_head=True, linestyle="-", \
        head_width=0.5, head_length=0.07)

    ax_arrows.annotate(args.positive_direction, xy=(args.pos_max_rho*0.80 - 0.01*len(args.positive_direction) - 0.015, -1.0), fontsize=16)
    ax_arrows.annotate(args.negative_direction, xy=(-args.neg_max_rho*0.80, -1.0), fontsize=16)

    ax_arrows.axis("off")

    leg_handles = [(mlines.Line2D([], [], color=color, marker=marker, linestyle='None', alpha=1.0, \
        markersize=args.dotsize, label=label)) for label,color,marker in zip(labels, colors, markers)]
 
    legend = ax.legend(handles=leg_handles, fontsize=16, loc=args.legloc, \
        frameon=True, markerfirst=True, ncol=1)

    ax.set_xlim([-args.neg_max_rho, args.pos_max_rho])
 
    ax.set_xlabel(args.x_axis, fontsize=16)

    ax.set_ylabel(args.y_axis, fontsize=16)

    if not args.outfile:
        args.outfile = args.metaanalysis[ 0 ].replace( ".tsv", "" )

    [plt.savefig("%s.%s" %(args.outfile, fmt), dpi=200) for fmt in ["svg", "png"]]




if __name__ == "__main__":

    draw_figure( read_params( sys.argv ) )
