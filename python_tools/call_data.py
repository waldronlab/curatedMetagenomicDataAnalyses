#!/usr/bin/env python

import argparse as ap
import numpy as np
import pandas as pd
import glob, os, sys
from skbio.stats import composition as cps

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
    add("-shr", "--shrink", action="store_true", help=\
	"In the case of any taxonomic feature, it just uses that level (gives you shorter names of feaures)")
    add("-sta", "--stratified", action="store_true", help=\
	"In case if pathways, pathcoverages and genefamilies, it returns also contributors for each feature [default: False]")
    add("-of", "--outfile", type=str, default="", help=\
	"If not set, will be a dump of metadata + the data-types inferred from the profiles")
    add("-a", "--apply_t", type=str, default="arcsin", choices=["CLR", "arcsin"], help=\
	"You can specify one of two transformation you want (CLR or arcsin squareroot). "
	"Note that if you already applied it before you musn't apply it twice (eve tough the program "
	"outputs a warning).")
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
    data_types = {"relative_abundance": "taxon", "pathway_abundance": "pathway", "gene_families": "genefamily", "pathway_coverage": "coverage"}
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
            
            if (data_type == "taxon") and (args.shrink):
                feature_mapper = dict([ (ft,ft.split("|")[-1]) for ft in Features if ft.split("|")[-1].startswith(feat_id) ])
                datatab.rename(index=feature_mapper, inplace=True)
                new_Features = [feature_mapper[ft] for ft in feature_mapper]
                datatab.drop([ft for ft in Features if (not ft in feature_mapper)], inplace=True, axis=0)
                Features = new_Features

            if (data_type in ["pathway", "genefamily", "coverage"]) and (not args.stratified):
                new_Features = [ft for ft in Features if (not "|" in ft)]
                datatab.drop([ft for ft in Features if (not ft in new_Features)], inplace=True, axis=0)
                Features = new_Features

            featuress = featuress | set(Features)
            Samples = metadata_this.index.tolist()
            datatab.columns = [(study + "_" + c) for c in datatab.columns.tolist()]
            ## Here I apply the transformation
            datatab = apply_transform_on_data(datatab.loc[ Features, Samples ], args.apply_t)
            ### HERE WE MUST INSERT A COUPLE OF CHECKS
            metadata_with_data = metadata_this.T.append( datatab )
            if not len(data):
                data = metadata_with_data
            else:
                data = data.merge(metadata_with_data, right_index=True, left_index=True, how="outer").fillna(0.0)
            sys.stdout.write("\n Metadata+data addition of study %s, of shape: " %study)
            sys.stdout.write(" ".join(map(str, list(metadata_with_data.shape))))

    data = data.loc[[i for i in data.index.tolist() if (not i in featuress)] + list(featuress)]
    data.index.name = "sampleID"

    if (not args.outfile):
        args.outfile = args.metadata.replace(".tsv", "_" + "_".join(data_type_merged) + ".tsv")

    sys.stdout.write("\n Final table to save...")
    sys.stdout.write("\n Number features = %i" %len(list(featuress)))
    sys.stdout.write("\n Number samples = %i" %len(data.columns.tolist()))
    sys.stdout.write("\n Tot.number metadata = %i" %(int(data.shape[0]) - len(list(featuress))))
    sys.stdout.write("\n ... saving in %s.\n" %args.outfile)
    data.to_csv(args.outfile, sep="\t", header=True, index=True)

if __name__ == "__main__":
    main(read_params())
