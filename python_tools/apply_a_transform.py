#!/usr/bin/env python

import argparse as ap
import numpy as np
import pandas as pd

def read_params():
    p = ap.ArgumentParser(description=\
	"This script takes a metadata + data matrix and apply a transformation "
	"to the data (CLR or arcsin-square-root (default)). It checks first if "
	"if one of these two might have been previously applied. \n"
        "NOTE: the std. featured identifier is s__. So, if you are working "
	"with ANY OTHER feature type, different from species-level relative "
	"abundances, remember to specify a different feature identifier.", \
	formatter_class=ap.RawTextHelpFormatter)
    add = p.add_argument
    add("datatable", type=str, help="metadata + data table. Data are expected with fields name "
	"on the index and samples as column names.")
    add("-a", "--apply_t", type=str, default="arcsin", choices=["CLR", "arcsin"], help=\
	"The script makes a check to test wether a transformation has been applied previously.")
    add("-z", "--feat_id", type=str, default="s__", choices=["k__", "p__", "o__", "c__", \
	"f__", "g__", "s__", "PWY", "UniRef90", "K"], help="Substring of any type "
	"uniquely identifing the features with respect to the metadata.")
    add("-of", "--outfile", type=str, default="", help=)
    return p.parse_args()

def apply_transform_on_data(dataframe, T):
    samples = dataframe.columns.tolist()
    features = dataframe.index.tolist()
    sums = np.sum(dataframe.values.astype(float), axis=0)
    if np.isclose(np.sum(sums), 100*len(samples)):
        sys.stdout.write("\n Relative abundances (summing to 100) detected. Proceeding...")
        if T.lower() == "clr":
            for sample in samples:
                dataframe[sample] = cps.clr(cps.multiplicative_replacement(dataframe.loc[features, sample].values.astype(float)))
        elif T.lower() == "arcsin":
            for sample in samples:
                dataframe[sample] = np.arcsin(np.sqrt(dataframe[sample].values.astype(float)/100.))
        else:
            raise NotImplementedError("Sorry: the transformation %s is not available yet." %T)
    elif np.isclose(np.sum(sums), len(samples)):
        sys.stdout.write("\n Relative abundances (summing to 1) detected. Proceeding...")
        if T.lower() == "clr":
            for sample in samples:
                dataframe[sample] = cps.clr(cps.multiplicative_replacement(dataframe.loc[features, sample].values.astype(float)*100.))
        elif T.lower() == "arcsin":
            for sample in samples:
                dataframe[sample] = np.arcsin(np.sqrt(dataframe[sample].values.astype(float)))
        else:
            raise NotImplementedError("Sorry: the transformation %s is not available yet." %T)
    elif np.isclose(np.sum(sums), 0.0):
        raise ValueError("\n Data summing to 0.0 with tol.: CLR already performed? Exiting...")
    elif np.isclose( np.sum([((np.sin(dataframe.loc[:, sample].values.astype(float)))**2)*100 \
        for sample in samples]), 100*len(samples)):
        raise ValueError("\n Data summing to 100 after inverse of arcsin-sqrt. Exiting...")
    return dataframe

def main(args):
    I = pd.read_csv(args.datatable, sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")
    feats = [i for i in I.index.tolist() if (args.feat_id in i)]
    metadata = [m for m in I.index.tolist() if (not m in feats)]
    dataframe = apply_transform_on_data(I.loc[feats, :], args.apply_t)
    metadata_dataframe = I.loc[metadata, :].append(dataframe)
    metadata_dataframe.index.name = "sampleID"
    if not args.outfile:
        args.outfile = args.datatable.replace(".tsv", "_%s.tsv" %args.apply_t)
    metadata_dataframe.to_csv(args.outfile, sep="\t", header=True, index=True)

if __name__ == "__main__":
    main(read_params())
