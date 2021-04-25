#!/usr/bin/env python

import argparse as ap
import numpy as np
import pandas as pd
import glob, os, sys

def read_params():
    ### we need to write somewgere that, if more tables
    ### for rhe same dataset are available it gets all of them
    p = ap.ArgumentParser(\
        description=\
        "Utility for retrieving cMD 3 data given "
        "a metadata table and a folder with profiles in it.", \
        formatter_class=ap.RawTextHelpFormatter)
    add = p.add_argument
    add("metadata", type=str, help="Metadata table of your analysis. "
        "The program guesses if it's row-wise or column-wise.")
    add("data_directory", type=str, help="Dir-prefix for the tables")
    add("-z", "--feat_id", type=str, help="Feature-identifier: if set, will use just this one. "
        "Otherwise uses a series of defaults deduced by the name of the file. "
        "Usefull mainly if you want genera (g__) as the default in relative_abundances "
        "is species (s__)")
    add("-re", "--release", type=str, default="2021-03-31") # 2021-03-31.BrooksB_2017.relative_abundance.tsv 
    add("-si", "--studyID", type=str, default="dataset_name")
    add("-of", "--outfile", type=str, default="")
    return p.parse_args()

def main(args):
    data_types = {"relative_abundance": "taxon", "pathway_abundance": "pathway", "gene_families": "genefamily"}
    default_feat_ids = {"taxon": "s__", "pathway": "PWY", "genefamily": "UniRef90", "KEGG": "K", "ECnumber": "EC"}
    metadata = pd.read_csv(args.metadata, sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")
    data = []
    data_type_merged = []
    featuress = set()
    if "number_reads" in metadata.index.tolist():
        metadata = metadata.T
    for dat in metadata[args.studyID].unique():
        metadata_this = metadata.loc[metadata[args.studyID]==dat]
        study = dat
        for datatable in glob.glob(os.path.join(args.data_directory, ".".join([args.release, study, "*", "tsv"]))):
            data_type = data_types[ os.path.basename(datatable).split(".")[2] ]
            feat_id = default_feat_ids[ data_type ]
            if (not data_type in data_type_merged):
                data_type_merged += [data_type]
            datatab = pd.read_csv(datatable, sep="\t", header=0, index_col=0, low_memory=False, engine="c")
            Features = [j for j in datatab.index.tolist() if (feat_id in j)]
            featuress = featuress | set(Features)
            Samples = metadata_this.index.tolist()
            datatab.columns = [(study + "_" + c) for c in datatab.columns.tolist()]
            ### HERE WE MUST INSERT A COUPLE OF CHECKS
            metadata_with_data = metadata_this.T.append( datatab.loc[ Features, Samples ] )
            if not len(data):
                data = metadata_with_data
            else:
                data = data.merge(metadata_with_data, right_index=True, left_index=True, how="outer")
            sys.stdout.write("\n Metadata+data addition of study %s, of shape: " %study)
            sys.stdout.write(" ".join(map(str, list(metadata_with_data.shape))))
    if (not args.outfile):
        args.outfile = args.metadata.replace(".tsv", "_".join(data_type_merged) + ".tsv")
    sys.stdout.write("\n Final table to save...")
    sys.stdout.write("\n Number features = %i" %len(list(featuress)))
    sys.stdout.write("\n Number samples = %i" %len(data.columns.tolist()))
    sys.stdout.write("\n Tot.number metadata = %i" %(int(data.shape[0]) - len(list(featuress))))
    sys.stdout.write("\n ... saving in %s.\n" %args.outfile)
    data.to_csv(args.outfile, sep="\t", header=True, index=True)

if __name__ == "__main__":
    main(read_params())    
