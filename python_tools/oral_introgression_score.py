#!/usr/bin/env python


import pandas as pd
import numpy as np
import argparse as ap
import sys, os


def read_params():
    p = ap.ArgumentParser(description=\
        "This program computes an oral-introgression score given "
        "a query, **gut** metagenomic dataset and 1 of  "
        "\n - an oral-only instance of curatedMetagenomicData 3"
        "\n - a list of oral-species."
        "\n The easiest solution is to use a precomputed list available in the directory where this program is stored."
        "\n The choice allows to use different list, and definitely theorize different scores."
        "\n If both choices or no-choices are specified, the program throws an error. "
        "The goal of the program is, besides, to compute such list."
        "\n The list of oral microbial species is used to estimate the introgression "
        "of oral species in the gut of the query samples. The Flag -c --entropy activates a "
        "a score computation modality which uses Shannon entropy (slightly more complicated "
        "but slightly more accurate according to our first tests)."
        "\n **Usage**"
        "\n a) python oral_introgression_score.py my_metagenomic_samples_toy_a.tsv --oral_cmd oral_cmd_exp_a.tsv"
        "\n After running this, you should find two files, called: "
        "\n my_metagenomic_samples_toy_a_oral_score.tsv & oral_cmd_exp_a_oral_species.tsv"
        "\n b) python oral_introgression_score.py my_metagenomic_samples_toy_b.tsv --transpose_input -os oral_cmd_exp_a_oral_species.tsv -c"
        formatter_class=ap.RawTextHelpFormatter)
    add = p.add_argument
    add("query_sample", type=str, \
        help="A table with Sample IDs as columns and index-names as microbial species")
    add("-ti", "--transpose_input", action="store_true", \
        help="With this flag, table is expected to have samples as index and microbial species as column-names")
    add("-z", "--featid", type=str, default="", \
        help="If activated, knows if a row/column is a species or not "
        "(default: thinks that all rows except the first are species)")
    add("-ocmd", "--oral_cmd", type=str, default="", \
        help="If passed, compute the list of oral features from a custom "
        "database, which follows the same defaults of the input query.")
    add("-os", "--oral_species", type=str, default="", \
        help="Compute the score for a [custom] list of species")
    add("-c", "--entropy", action="action_true", \
        help="Uses the entropy score (slightly more complicated, according to the first tests more accurate)")
    return p.parse_args()


def take_input(args):
    input_table = pd.read_csv(args.query_sample, sep="\t", header=0, index_col=0, low_memory=False, engine="c")
    if args.transpose_input:
        input_table = input_table.T
    query_samples = input_table.columns
    if args.featid:
        query_features = [ i for i in input_table.index if args.featid in i ]
    else:
        query_features = input_table.index.tolist()
    return input_table.loc[query_features]


def evaluate_oral_species(args):
    if (not args.oral_cmd) and (not args.oral_species):
        raise ("One of oral_cmd and oral_species must be activated. Exiting")
    if args.oral_cmd and args.oral_species:
        raise ("Arguments oral_cmd and oral_species are meant to be mutually exclusive."
            "Specifying the first serves to evaluate and save an istance of the second.")
    if args.oral_cmd:
        oral_cmd = pd.read_csv(args.oral_cmd, sep="\t", header=0, index_col=0, low_memory=False, engine="c")
        if args.transpose_inputs:
            oral_cmd = oral_cmd.T
        if args.featid:
            oral_features = [ i for i in input_table.index if args.featid in i ]
        else:
            oral_features = input_table.index.tolist()
        oral_species_list = [i for i in oral_features if \
            ((np.mean(oral_cmd.loc[i].values.astype(float))>0.01) and \
            ((np.count_nonzero(oral_cmd.loc[i].values.astype(float))/float(len(oral_features)))>0.05))]
        with open(args.oral_cmd.split(".")[0] + "_oral_species.tsv", "w") as opl:
            [opl.write(species + "\n") for species in oral_species_list]
    elif args.oral_species:
        with open(args.oral_species) as osp:
            oral_species_list = [line.rstrip() for line in osp.readlines()]
    return oral_species_list


def compute_oral_score(samples_frame, oral_species, simple=True):
    oral_species = [o for o in oral_species if (o in samples_frame.index)]
    samples_frame = samples_frame.loc[ oral_species, : ]
    if simple:
        give_score = lambda array : np.sum(array)
        score_name = "Oral-totabundance-score"
    else:
        give_score = lambda array : -np.sum([ ( (p/100.) * ( np.log( p/100 ) if (p>0.) else 0.  )) for p in array ])
        score_name = "Oral-entropy-score"
    ScoreS = pd.DataFrame({\
        "sample_id": samples_frame.columns.tolist(), \
        score_name: [ give_score( samples_frame.loc[:, sample] ) for sample in samples_frame.columns]\
        })
    return ScoreS


def main():
    args = read_params()
    query = take_input(args)
    oral_species = evaluate_oral_species(args)
    scores = compute_oral_score(query, oral_species, simple=(not args.entropy))
    scores.to_csv(args.query_sample.split(".")[0] + "_oral_score.tsv", sep="\t", header=True, index=True)
    

if __name__ == "__main__":
    main()
    print("Have a nice day!")
