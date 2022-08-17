#!/usr/bin/env python


import pandas as pd
import numpy as np
import os, argparse as ap
import sys
import matplotlib
#matplotlib.use('Agg')
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
import warnings
import matplotlib.ticker as mticker
from num2words import num2words


sns.set_style("whitegrid")
sns.set_context("paper", rc={"font.size":24,"axes.titlesize":24,"axes.labelsize":24, "axes.ticklabelssize": 24})
## sns.set(font="Arial")


def read_params( args ):
    p = ap.ArgumentParser(description=\
        "", \
        formatter_class=ap.RawTextHelpFormatter)
    add = p.add_argument
 
    add("metaanalysis", type=str, nargs="+", help="A meta-analysis with features as rows and fields as columns")

    add("--names", type=str, nargs="+", default=[""], help="Toy can use this to specify names that will be used in the legend."
	"Is useful when p.otting multiple analyses together. Default [None] will switch automatically to naming "
        "analyses with an ordinal number")
 
    add("-n", "--narrowed", action="store_true", help="Only considers confidence intervals which don't go over zero.")

    add("-b", "--boxes", action="store_true", help="Draws also boxplots, of the distribution of the populations.")

    add("-a", "--relative_abundances", type=str, default=[""], nargs="+", help="A file with features as rows and samples as columns")

    add("-i", "--imp", type=int, default=30, help="Tries to plot the first *imp*ortant features. ")

    add("-how", "--how", type=str, default="first", choices=["first", "union", "intersection", "joint_top"], help=\
        "How to chose the [--imp] features to plot, default : intersection. NOTE: "
        "meaningfull only for multiple meta-analysis plot")

    add("-ps", "--positive_direction", type=str, default="Positive side", help=\
        "Label of the +1 group (with arrow)")

    add("-ns", "--negative_direction", type=str, default="Negative side", help=\
        "Label of the -1 group (with arrow)")

    add("-x", "--x_axis", type=str, default="Effect-size")

    add("-y", "--y_axis", type=str, default="")

    add("-t", "--title", type=str, default="")

    add("-es", "--e_suff", type=str, default=["_es"], nargs="+")
 
    add("-qs", "--q_suff", type=str, default=["_Q"], nargs="+")

    add("-ses", "--se_suff", type=str, default=["_SE"], nargs="+")
 
    #add("-pe", "--prevalence", type=str, default=[""], nargs="+")

    add("-mp", "--min_prevalence", type=float, default=[0.00], nargs="+")

    add("-mna", "--min_ab", type=float, default=[0.0], nargs="+")

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

    add("-rs", "--random_effect_se", type=str, nargs="+", default=["RE_StdErro"], help=\
	"Random/fixed effect model SERR")

    add("-cr", "--color_red", type=str, nargs="+", default=["goldenrod"], help=\
        "Color of the significant results and of the two important lines marking the desired thresholds")

    add("-cb", "--color_blue", type=str, nargs="+", default=["dodgerblue"], help=\
        "Color the non-significant results and of the zero-marking line ")

    add("-cbl", "--color_black", type=str, nargs="+", default=["black"], help=\
        "Color of the random/fixed effect diamond ")

    add("-pmr", "--pop_marker", type=str, nargs="+", default=["o"], help=\
        "Marker of the single population-effect [default=Circle]")

    add("-dmr", "--diam_marker", type=str, nargs="+", default=["D"], help=\
        "Marker of the random/fixed effect [default=Diamond]")

    add("-il", "--important_lines", type=float, nargs="+", default=[0.3], help=\
        "Two lines are draw at this level to mark good results [default = 0.3]")

    add("-as", "--a_single", type=float, default=0.2, \
        help="Significant alpha for the single analyses/populations")

    add("-ar", "--a_random", type=float, default=0.05, \
        help="Significant alpha for the random/fixed effect model")

    add("--dotsize", type=int, default=9, help="Default = 9")

    add("--diamsize", type=int, default=17, help="Default = 17")

    add("--neg_max_rho", type=float, default=0.80)

    add("--pos_max_rho", type=float, default=0.80)

    add("--legloc", type=str, default="best")

    add("--check_genera", action="store_true", help="meaningfull only with --how==joint_top; it searches that if genera and species are present genera are"
            " selected only in case they are different at the 3rd digimal digit")
 
    add("--shrink_names", type=str, default="|", help="Apply .split(\"<shrink-names>\") to the Y-axis (feature) names of the figure")

    add("--size_on_err", action="store_true", help="Size the dots on the SE of the population")

    print(vars(p.parse_args()))

    return vars(p.parse_args())





def parse_args_table(args, ij):
 
    take_a_param = lambda param, ij : param[ij] if len(param)>1 else param[0]

    if args["names"][0] and (len(args["names"])!=len(args["metaanalysis"])):
        warnings.warn("You specify names of the analysis (which is used for the legend), "
	    "but the number of names specified don't match the number of analysis. "
	    "Switching to ordinal names.")
    
    name = num2words(ij, ordinal=True) if (not take_a_param( args["names"], ij )) else take_a_param( args["names"], ij )

    #if isinstance( args, dict ):
    #    from attrdict import AttrDict
    #    args = AttrDict(args)

    meta_analysis = pd.read_csv( take_a_param( args["metaanalysis"], ij ), sep="\t", header=0, index_col=0, low_memory=False, engine="c")

    meta_analysis.fillna( "NA", inplace=True )

    filter_on = "NO_FILTER"
 
    if (take_a_param( args["relative_abundances"], ij)): # and (not take_a_param( args["prevalence"], ij )):
        #warnings.warn("You didn't specifiy neither abundances tables nor prevalence. This means that all the feats found will be retained.")
        #filter_on = "NO_FILTER"

        if take_a_param( args["relative_abundances"], ij) and take_a_param( args["min_prevalence"], ij ):
            filter_on = "PREVALENCE"

            #warnings.warn("You set both prevalence AND relative_abundances: the first one will be used")
        elif take_a_param( args["relative_abundances"], ij) and (take_a_param( args["min_ab"], ij ) > 0.0):
            filter_on = "BOTH" if take_a_param( args["min_prevalence"], ij ) else "ABUNDANCE"

    min_prevalence = dict( )
    features = meta_analysis.index.tolist( )

    if take_a_param( args["relative_abundances"], ij ):
        abundances = pd.read_csv( take_a_param( args["relative_abundances"], ij), sep="\t", header=0, index_col=0, low_memory=False, engine="c")
        prevalences = dict([ (feat, (np.count_nonzero(abundances.loc[feat, :].values.astype(float))/float(abundances.shape[1]) \
            if feat in abundances.index.tolist() else 1.0) ) for feat in features ])
  
    #if take_a_param( args["prevalence"], ij ): # and not take_a_param( args.min_ab, ij )):  
    #    prevalences = dict([ (feat, meta_analysis.loc[feat, take_a_param( args["prevalence"], ij )]) for feat in features])

    single_effects = [ c for c in meta_analysis.columns.tolist() if \
        ((c.endswith( take_a_param( args["e_suff"], ij ) )) and (c != take_a_param( args["random_effect"], ij ) ))  ]
 
    single_effect_Qs = [ c for c in meta_analysis.columns.tolist() if \
        ((c.endswith( take_a_param( args["q_suff"], ij ) )) and (c != take_a_param( args["random_effect_q"], ij ) ))  ]
  
    single_std_errors = [ c for c in meta_analysis.columns.tolist() if \
	((c.endswith( take_a_param( args["se_suff"], ij) )) and (c != take_a_param( args["random_effect_se"], ij) ))  ]
 
    if take_a_param( args["random_effect"], ij ) in single_effects:
        single_effects = single_effects.remove( take_a_param( args["random_effect"], ij ) )

    if take_a_param( args["random_effect_q"], ij ) in single_effect_Qs:
        single_effect_Qs = single_effect_Qs.remove( take_a_param( args["random_effect_q"], ij ) )
        
    results = meta_analysis.loc[ meta_analysis[ take_a_param( args["random_effect_q"], ij ) ]<\
            (args["a_random"] if (args["how"] in ["first", "joint_top"]) else 1.0), : ]

    results = results.loc[ [ i for i in results.index.tolist() if ((len(single_effects) - results.loc[ i, single_effects ].tolist().count("NA")) \
        >= take_a_param( args["min_studies"], ij)) ], :]

    if filter_on != "NO_FILTER":

        if filter_on == "PREVALENCE":
            results = results.loc[ [ i for i in results.index.tolist() if (prevalences[i]>=take_a_param( args["min_prevalence"], ij))], : ]

        elif filter_on == "ABUNDANCE":
            results = results.loc[ [ i for i in results.index.tolist() if (np.mean(abundances.loc[i].values.astype(float))>=take_a_param( args["min_ab"], ij)) ], : ]

        else:
            results = results.loc[ [ i for i in results.index.tolist() if (prevalences[i]>=take_a_param( args["min_prevalence"], ij))], : ]
            results = results.loc[ [ i for i in results.index.tolist() if (np.mean(abundances.loc[i].values.astype(float))>=take_a_param( args["min_ab"], ij)) ], : ]
                 
            #    (not i in abundances.index.tolist()) or (np.mean(abundances.loc[i].values.astype(float))>=take_a_param( args["min_prevalence"], ij)) \
            #    ], : ]
    
    if args["narrowed"]:
        results = results.loc[[i for i in results.index.tolist() if \
            (float(results.loc[i, take_a_param( args["confint"], ij) ].split(";")[0])*float(results.loc[i, take_a_param( args["confint"], ij)].split(";")[1])) \
            >0.0], :]

    results["abs"] = np.abs( results[ take_a_param( args["random_effect"], ij) ].values.astype(float) )
    results.sort_values( "abs", inplace=True, ascending=False )

    if args["how"] in ["first", "joint_top"]:
        results = results[:args["imp"]]

    positive = [i for i in results.index.tolist() if ( float(results.loc[i, take_a_param( args["random_effect"], ij)] >0.) ) ]
    negative = [i for i in results.index.tolist() if ( float(results.loc[i, take_a_param( args["random_effect"], ij)] <0.) ) ]
    res_pos = results.loc[ positive, : ].sort_values( take_a_param( args["random_effect"], ij), ascending=False )
    res_neg = results.loc[ negative, : ].sort_values( take_a_param( args["random_effect"], ij), ascending=True )
    results = res_pos.append(res_neg)

    if args["shrink_names"]:
       results.index = [i.split(args["shrink_names"])[-1] for i in results.index ]

       #ax.set_yticklabels([str(x.get_text()).split(args["shrink_names"])[-1] for x in ax

    return [results, single_effects, single_effect_Qs, single_std_errors, results.index.tolist( ), name]
            


def select_features(how, results):
    if how == "first":
        return results[0].index.tolist()
    elif how == "union":
        return set().union(*results)
    elif how == "intersection":
        return set().intersection(*results)


  
def select_joint_top_features( args, results ):
    take_a_param = lambda param, ij : param[ij] if len(param)>1 else param[0]
    effects = sorted(
            list(it.chain.from_iterable([ [(k,abs(v),v) for k,v in curr_res.loc[ :, take_a_param( args["random_effect"], i ) ].to_dict().items()] for i,curr_res in enumerate(results) ])), 
            key=lambda a : a[1], reverse=True)

    to_rem = []

    if args["check_genera"]:
        for name in effects:
            e = name[0]
            if \
            (e.startswith("g__")) and (len([ x for x in effects if (x[0].startswith("s__") and x[0].replace("s__", "g__").startswith(e)) ]) == 1):
            #(any([(x[0].startswith(e.replace("g__", "s__"))) for x in effects if (x[0]!=e) ])):  #and \
            #(len([ x for x in effects if (x[0].startswith("s__") and x[0].replace("s__", "g__").startswith(e)) ]) == 1):
 
            #    effect = str(name[1])
            #    spc_l = [ (x[0]) for x in effects if ( (x[0]!=e) and (x[0].startswith(e.replace("g__", "s__"))) ) ]
            #    spc_e = [ (x[1]) for x in effects if ( (x[0]!=e) and (x[0].startswith(e.replace("g__", "s__"))) ) ]
            
            #    mapp = dict([(x,str(y)) for x,y in zip(spc_l, spc_e)])
            #    for x in mapp:
            #        if mapp[x].split(".")[1][:3] == effect.split(".")[1][:3]:
                        to_rem.append( e )

    effects = [u for u in effects if (not u[0] in to_rem)]
    result_feats = [w[0] for w in sorted([x for x in effects[:args["imp"]]], key=lambda a : a[2])]

    #if args["shrink_names"]:
     #   result_feats = [x.split(args["shrink_names"])[-1] for x in result_feats]

    return result_feats


    

 
def build_long_frame(args, analysis, ij, result_features, all_markers):
 
    results, single_effects, single_effect_Qs, single_std_errors, Feats__, name_of_this_one = tuple(analysis)

    Feats = result_features #[feat for feat in result_features if feat in results.index]

    take_a_param = lambda param, ij : param[ij] if len(param)>1 else param[0]

    get = lambda df, index, column : np.nan if (not index in df.index.tolist()) else df.loc[index, column]

    #print( len(list( it.chain.from_iterable( [ [ feature for e in range(len(single_effects)) if \
    #        str(get( results, feature, single_effects[e] ))!="NA"] for feature in Feats ] ) )) )

    #print(len(  list( it.chain.from_iterable( [ [ float(get(results, feature, e ))  for e in single_effects if \
    #        str(get(results, feature, e))!="NA"] for feature in Feats ] ) )    ))
    #print(len(  list( it.chain.from_iterable( [ [ e.replace(take_a_param( args["e_suff"], ij), "")  for e in single_effects if \
    #        str(get(results, feature, e))!="NA"] for feature in Feats ] ) )    ))
    #print(len(  list( it.chain.from_iterable( [ [ float(get(results, feature, e ))  for e in single_effect_Qs if \
    #        str(get(results, feature, e))!="NA"] for feature in Feats ] ) )    ))

    lgf_singles = pd.DataFrame(\
        {
        "Feature": list( it.chain.from_iterable( [ [ feature for e in range(len(single_effects)) if \
            str(get(results, feature, single_effects[e]))!="NA"] for feature in Feats ] ) ) , \
        "effect-size": list( it.chain.from_iterable( [ [ float(get(results, feature, e ))  for e in single_effects if \
            str(get(results, feature, e))!="NA"] for feature in Feats ] ) ), \
        "population": list( it.chain.from_iterable( [ [ e.replace(take_a_param( args["e_suff"], ij), "")  for e in single_effects if \
            str(get(results, feature, e))!="NA"] for feature in Feats ] ) ),
        "q-values": list( it.chain.from_iterable( [ [ float(get( results, feature, e ))  for e in single_effect_Qs if \
            str(get(results, feature, e))!="NA"] for feature in Feats ] ) ),
	"std-errors": list( it.chain.from_iterable([ [ float(get( results, feature, e ))  for e in single_std_errors if \
	    str(get(results, feature, e))!="NA"] for feature in Feats ] ))
        } )
    
    if args["markers"]:
        if isinstance(all_markers, int):
            all_markers = [',', '.', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
            # [".", ",", "v", "^", "<", ">", "s", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d", "8"]
        pop_to_marker = {}
        for pop in lgf_singles["population"].unique().tolist():
            if not pop in pop_to_marker:
                pop_to_marker[pop] = all_markers[-1]
                all_markers = all_markers[:-1]
 
        lgf_singles["marker_of_population"] = [pop_to_marker[pop] for pop in lgf_singles["population"].tolist()]

    remaining_markers = all_markers

    lgf_singles.index = lgf_singles["Feature"].tolist()
    lgf_singles["SG"] = [ ("yes" if (q<args["a_single"]) else "no") for q in lgf_singles["q-values"].tolist() ]
 
    ## confint = [ get(results, feature, take_a_param( args["confint"], ij) ).split(";") for feature in Feats ]
  
    confint = [ ([np.nan, np.nan] if (not feature in results.index.tolist()) else \
	results.loc[feature, take_a_param( args["confint"], ij)].split(";")) \
	for feature in Feats ]

    lgf_effect = pd.DataFrame(\
        { "Feature": Feats, \
          "random-effect": [ float(get(results, feature, take_a_param( args["random_effect"], ij) )) for feature in Feats ], \
          "Q-values": [ float(get( results, feature, take_a_param( args["random_effect_q"], ij) )) for feature in Feats ], \
          "CI_low":  [ float(ci[0]) for ci in confint ], \
          "CI_upp":  [ float(ci[1]) for ci in confint ], \
        } \
        )

    lgf_effect.index = lgf_effect["Feature"].tolist()
    lgf_effect["SG"] = [ ("yes" if (q<args["a_random"]) else "no") for q in lgf_effect["Q-values"].tolist() ]
    
    return lgf_singles, lgf_effect, name_of_this_one, remaining_markers





def draw_figure(args, show):

    if not show:
        matplotlib.use('Agg')

    analyses = [ \
        parse_args_table( args, ij ) for ij in range(len( args["metaanalysis"] ))
        ]

    if args["how"] != "joint_top":
        result_features = select_features( args["how"], [ an[0] for an in analyses ] )
    else:
        result_features = select_joint_top_features( args, [ an[0] for an in analyses ] )
 
    take_a_param = lambda param, ij : param[ij] if len(param)>1 else param[0]
 
    labels, colors, markers = [], [], []

    # [results, prevalences, single_effects, single_effect_Qs, results.index.tolist( )]

    fig = plt.figure(figsize=(14, 16))

    gs = gridspec.GridSpec( 2,1, height_ratios=[1]+[len(result_features)] )

    ax = plt.subplot( gs[ 1,0 ] )

    ax_arrows = plt.subplot( gs[ 0,0 ])
  

    for ij,analysis in enumerate( analyses ):

        if not ij: 
            long_form_frame_singles, long_form_frame_effect, name_of, remaining_markers = build_long_frame(args, analysis, ij, result_features, -1 )
        else:
            long_form_frame_singles, long_form_frame_effect, name_of, remaining_markers = build_long_frame(args, analysis, ij, result_features, remaining_markers )


        #long_form_frame_effect.to_csv("mia_nonna_%i.tsv" %ij, sep="\t", header=True, index=True)

        #long_form_frame_singles.index = [(str(i)+k) for i,k in enumerate(long_form_frame_singles.index)]
        #long_form_frame_singles.drop([i for i in long_form_frame_singles.index if \
        #    (np.sum(long_form_frame_singles.loc[i, ["effect-size", "q-values"]].values.astype(float))==0.)], inplace=True)
        #long_form_frame_singles.index = [k[1:] for k in long_form_frame_singles.index]

        #print(long_form_frame_effect, "BEFORE")

        #long_form_frame_effect.index = [(str(i)+k) for i,k in enumerate(long_form_frame_effect.index)]
        #long_form_frame_effect.drop([i for i in long_form_frame_effect.index if \
	#    (np.sum(long_form_frame_effect.loc[i, ["random-effect", "CI_low", "CI_upp"]].values.astype(float))==0.)], inplace=True)
        #long_form_frame_effect.index = [k[1:] for k in long_form_frame_effect.index]

        #rint(long_form_frame_effect, "AFTER")


        if not args["markers"]:
            long_form_frame_singles["dot_is"] = "the_population_marker" #take_a_param(args["pop_marker"], ij)

            if not args["size_on_err"]:
                strip_one = sns.scatterplot(\
                    x="effect-size", y="Feature", \
                    data=long_form_frame_singles, size="SG", \
                    sizes={"no": args["dotsize"]*8, "yes": args["dotsize"]*16}, \
                    palette={"yes": take_a_param(args["color_red"], ij), \
                    "no": take_a_param(args["color_blue"], ij)}, \
                    edgecolor="black", ax=ax, hue="SG", \
                    style="dot_is", \
                    markers={"the_population_marker": take_a_param(args["pop_marker"], ij)})
            else:
                strip_one = sns.scatterplot(\
                    x="effect-size", y="Feature", \
                    data=long_form_frame_singles, size="std-errors", \
                    palette={"yes": take_a_param(args["color_red"], ij), \
                    "no": take_a_param(args["color_blue"], ij)}, \
                    edgecolor="black", ax=ax, hue="SG", \
                    markers={"the_population_marker": take_a_param(args["pop_marker"], ij)})

        else:
            marker_of_this_population = dict([(p,m)\
                for p,m in zip(\
                long_form_frame_singles["population"].tolist(), \
                long_form_frame_singles["marker_of_population"].tolist() )]\
            )

            for p,m in marker_of_this_population.items():
                labels += [ p ]
                markers += [ m ]
                colors += [ "darkgrey" ]

            #if not args["size_on_err"]:
            #    print("SONO PALESEMENTE AFFANCULOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO")
 
            strip_one = sns.scatterplot(\
                    x="effect-size", y="Feature", \
                    data=long_form_frame_singles, size="SG", \
                    sizes={"no": args["dotsize"]*8, "yes": args["dotsize"]*16}, \
                    palette={"yes": take_a_param(args["color_red"], ij), \
                    "no": take_a_param(args["color_blue"], ij)}, \
                    edgecolor="black", ax=ax, hue="SG", \
                    style="population", markers=marker_of_this_population)

            #else:
            #    print(" SONOEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa")
            #strip_one = sns.scatterplot(\
            #        x="effect-size", y="Feature", \
            #        data=long_form_frame_singles, size="std-errors", \
            #        palette={"yes": take_a_param(args["color_red"], ij), \
            #        "no": take_a_param(args["color_blue"], ij)}, \
            #        edgecolor="black", ax=ax, hue="SG", \
            #        markers=marker_of_this_population)

        if args["boxes"]:
            strip_b = sns.boxplot(\
                x="effect-size", y="Feature", \
                data=long_form_frame_singles, \
                saturation=0.25, color="palegoldenrod", ax=ax, width=0.45)

        strip_RE = sns.stripplot(\
            x="random-effect", y="Feature", data=long_form_frame_effect, \
            marker=take_a_param(args["diam_marker"], ij), ax=ax, color=take_a_param(args["color_black"], ij), \
            size=args["diamsize"], palette={"no": "goldenrod", "yes": "black"}, \
            hue="SG")
 
        for e,feat in enumerate(long_form_frame_effect["Feature"].tolist()):
            lw, upp = float(long_form_frame_effect.loc[feat, "CI_low"]), float(long_form_frame_effect.loc[feat, "CI_upp"])
            ax.plot([ lw, lw, lw, upp, upp, upp ], [ e-0.01, e+0.01, e, e, e+0.01, e-0.01 ], \
            c=take_a_param(args["color_red"], ij), linewidth=3.6)

        blue_line = ax.axvline(ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], linewidth=3.5, \
            color=args["color_blue"][0], linestyle="--", alpha=1.0)

        red_line_1 = ax.axvline(x=-take_a_param(args["important_lines"], ij), ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], \
            linewidth=2.0, color=take_a_param(args["color_red"], ij), linestyle="--", alpha=1.0)

        red_line_2 = ax.axvline(x=take_a_param(args["important_lines"], ij), ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], \
            linewidth=2.0, color=take_a_param(args["color_red"], ij), linestyle="--", alpha=1.0)
  

        if len(long_form_frame_effect["SG"].unique().tolist())==2:
            labels += [ name_of, name_of ]
            markers += [ take_a_param(args["diam_marker"], ij), take_a_param(args["diam_marker"], ij) ]
            colors += [ take_a_param(args["color_black"], ij), "yellow" ]

        else:
            labels += [ name_of ]
            markers += [ take_a_param(args["diam_marker"], ij) ]
            colors += [ take_a_param(args["color_black"], ij) ]
 
        labels += [ ("%s-pop non sign." %name_of), ("%s-pop Sign." %name_of) ]
        markers += [ take_a_param(args["pop_marker"], ij), take_a_param(args["pop_marker"], ij) ]
        colors += take_a_param(args["color_blue"], ij), take_a_param(args["color_red"], ij)


    ax_arrows.set_xlim([-args["neg_max_rho"], args["pos_max_rho"]])
    ax_arrows.set_ylim([-1.0,1.0])

    leg_pos = ax_arrows.arrow(x=0.02, y=0.5, dx=args["pos_max_rho"]*0.80, dy=0.0, color=args["color_blue"][0], \
        shape="full", width=.075, length_includes_head=True, linestyle="-", \
        head_width=0.5, head_length=0.07)

    leg_neg = ax_arrows.arrow(x=-0.02, y=0.5, dx=-args["neg_max_rho"]*0.80, dy=0.0, color=args["color_blue"][0], \
        shape="full", width=.075, length_includes_head=True, linestyle="-", \
        head_width=0.5, head_length=0.07)


    #print(args["positive_direction"], len(args["positive_direction"]))
    #print(1.0*float(len(args["positive_direction"])))

    ax_arrows.annotate(args["positive_direction"], xy=(args["pos_max_rho"]*0.80 - 0.01*len(args["positive_direction"]) - 0.015, -1.0), fontsize=24)
    ax_arrows.annotate(args["negative_direction"], xy=(-args["neg_max_rho"]*0.80, -1.0), fontsize=24)

    ax_arrows.axis("off")
 
    labels_, colors_, markers_ = [], [], []
    duplicates_question = [] 
    for label,color,marker in zip(labels, colors, markers):
        if not (label,color,marker) in duplicates_question:
            duplicates_question += [tuple((label,color,marker))]
    #print(duplicates_question, "duplicated")

    for ob in duplicates_question:
        labels_ += [ob[0]]
        colors_ += [ob[1]]
        markers_ += [ob[2]]

    #print(labels_)
    #print(colors_)
    #print(markers_)

    leg_handles = [(mlines.Line2D([], [], color=color, marker=marker, linestyle='None', alpha=1.0, \
        markersize=args["dotsize"], label=label)) for label,color,marker in zip(labels_, colors_, markers_) ]
 
    legend = ax.legend(handles=leg_handles, fontsize=24, loc=args["legloc"], \
        frameon=True, markerfirst=True, ncol=1)

    ax.set_xlim([-args["neg_max_rho"], args["pos_max_rho"]])
 
    ax.set_xlabel(args["x_axis"], fontsize=24)

    ax.set_ylabel(args["y_axis"], fontsize=24)

    if not args["outfile"]:
        args["outfile"] = args["metaanalysis"][ 0 ].replace( ".tsv", "" )

    ax_arrows.set_title(args["title"], fontsize=24)

    xticks_loc = ax.get_xticks().tolist()
    ax.xaxis.set_major_locator(mticker.FixedLocator(xticks_loc))
    ax.set_xticklabels( [ ("%.2f" %n) for n in xticks_loc ], fontsize=24 )

    if args["shrink_names"]: 
 
       #ax.set_yticklabels([str(x.get_text()).split(args["shrink_names"])[-1] for x in ax.get_yticklabels()], fontsize=24)
     
       for lab in ax.get_yticklabels():
           lab.set_style("italic")

    if not show:
        [plt.savefig("%s.%s" %(args["outfile"], fmt), dpi=200) for fmt in ["svg", "png"]]
    else:
        plt.show()




if __name__ == "__main__":

    draw_figure( read_params( sys.argv ), show=False )
