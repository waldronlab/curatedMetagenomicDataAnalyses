#!/usr/bin/env python

import pandas as pd
import numpy as np
import os, sys
import argparse as ap
from scipy import stats as sts
from scipy.spatial.distance import pdist
from skbio.stats import distance as dts


def read_params():
    p = ap.ArgumentParser()
    add = p.add_argument
    add("dataset", type=str, help="(K (features) + M (metadata)) * N (samples) table. "
        "Metadata are detected with the FeatID arg.")
    add("groups", type=str, help="Determines on what the different of the groups is computed. "
        "NaN are excluded.")
    add("-si", "--study_id", type=str, default="study_identifier", help="Study identifier")
    add("-z", "--feat_id", type=str, default="s__", help="A substring present in all the features. "
        "Default (s__) is the substring present in all metaphlan 3 microbial species.")
    add("-a", "--analysis", type=str, default="permanova", choices=["permanova", "anosim"])
    add("-m", "--metric", type=str, default="euclidean")
    add("-s", "--significance", type=float, default=0.05)
    return p.parse_args()


def get_input(args):
    I = pd.read_csv(args.dataset, sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")
    feats = [i for i in I.index.tolist() if (args.feat_id in i)]
    I = I.loc[ :, I.loc[args.groups]!="NA" ]
    I = I.loc[ [args.study_id, args.groups] + feats ].astype(float, errors="ignore")
    return I


def store_analyze(data_dic, analy, metric, study_id, groups):
    for study in data_dic:
        data_dic[study]["DISTMAT"] = dts.DistanceMatrix( pdist( data_dic[study]["DATA"].drop([ study_id, groups ]).T, metric=metric ) )
        data_dic[study]["STATS"] = analy( data_dic[study]["DISTMAT"], grouping=data_dic[study]["DATA"].loc[groups] )
        data_dic[study]["PVALUE"] = data_dic[study]["STATS"].loc["p-value"]
    return data_dic


def main():
    args = read_params()

    if args.analysis == "permanova":
        func = dts.permanova
        hie = "BINOM"

    elif args.analysis == "anosim":
        func = dts.anosim
        hie = "FISHER"

    indat = get_input(args)
    indat_di = dict([ (dat, {"DATA": indat.loc[ :, indat.loc[ args.study_id ]==dat  ]}) for dat in indat.loc[ args.study_id ].unique()  ])

    analyses = store_analyze( indat_di, func, args.metric, args.study_id, args.groups )
 
    for k in analyses:
        print("\n============\n", k, analyses[k]["STATS"])

    print("\n=========\n")
    print("BINOM P: ", sts.binom_test( \
        len([ analyses[study]["PVALUE"] for study in indat.loc[ args.study_id ].unique() \
        if (analyses[study]["PVALUE"]<args.significance) ]), \
        len([ analyses[study]["PVALUE"] for study in indat.loc[ args.study_id ].unique()]), \
        p=args.significance
        ))



if __name__ == "__main__":
    main()
