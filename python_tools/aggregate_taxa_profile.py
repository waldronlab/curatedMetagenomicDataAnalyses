#!/usr/bin/env python


import pandas as pd
import numpy as np
import argparse as ap
import sys, os


def read_params():
    p = ap.ArgumentParser()
    add = p.add_argument
    add("metaphlan_profile", type=str, help=\
        "This must be a metaphlan-species table to aggregate. "
        "It can have the metadata above, but MUST HAVE SAMPLES in COLUMNS.")
    add("-l", "--level", type=str, default="g", choices=\
        ["s", "g", "f", "o", "c", "p", "k", "S", "G", "F", "O", "C", "P", "K"], \
        help="The level to which you want to aggregate. Default is g. "
        "Other levels are, in order: k, c, p, o, f")
    add("-o", "--outfile", type=str, default="", help=\
        "If not specified, it uses the name of the input and attacches the name of the level")
    return p.parse_args()


def aggregate_table(table_name, level):

    aggregation_levels_dic = {"s__": 6, "g__": 5, "f__": 4, "o__": 3, "c__": 2, "p__": 1, "k__": 0}
    table = pd.read_csv(table_name, sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")
    N = table.shape[1]

    species = [j for j in table.index.tolist() if ("s__" in j)]
    aggregated_table = table.copy().loc[~species]

    desired_level = list(  set(  [ i for i in species if (len(i.split("|"))==(aggregation_levels_dic[ level ]+1))  ]  ) )
    desired_level_ab = dict([(taxon, np.zeros(N)  ) for taxon in desired_level])

    for taxon in species:
        if len( taxon.split( "|" ) ) > ( aggregation_levels_dic[ level ] +1 ):
            desired_level_ab[ taxon.split( "|" )[ aggregation_levels_dic[ level ] ] ] += table.loc[ taxon  ].values.astype( float )

    for aggregated_taxon in desired_level_ab:
        aggregated_table.loc[ aggregated_taxon ] = desired_level_ab[ aggregated_taxon ]

    return aggregated_table

    
def main():
    args = read_params()
    args.level = args.level.lower()

    if args.level == "s":
        raise("The level required (species) is the standard level of the table, "
            "and therefore can't be aggregated")

    if not args.outfile:
        args.outfile = args.metaphlan_profile.split(".")[0] + args.level + "_level.tsv"

    aggregated_table = aggregate_table(args.metaphlan_profile, args.level + "__")
    aggregated_table.to_csv(args.outfile, sep="\t", header=True, index=True)


if __name__ == "__main__":
    main()
