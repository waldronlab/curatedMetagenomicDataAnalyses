#!/usr/bin/env python

import argparse as ap
import numpy as np
import pandas as pd
import glob, os, sys
from scipy import stats as sts
 
def read_params():
    p = ap.ArgumentParser(description=\
        "Utility for creating a metadata table to perform a meta-analysis. \n"
        "This program has two main scopes: a) if you specify as the input_folder argument \n"
        "a combinedMetadata table, it will merge the population queried starting from that \n"
        "and therefore will store no-data. In this case, data will be \n"
        "retrievable using another script: '[python] call_data.py' which takes \n"
        "a metadata table and a type of data as input and ensures "
        "you get your data. The other usage for this program is "
        "\nto pass as a first argument a folder which contains metadata + data. "
        "\nIn this second settings, the program will subset the load of \n"
        "data i input_folder for you, returning a metadata + data table directly."
	"\n\n*** Example usages *** \n\n"
        "The following 2 examples return \nthe datasets for 2 analyses "
        "described in the paper \n'Manghi, Schiffer, et al, 2021 in press.'\n"
        ""
        "\nExample a):"
        "\n==> returns 4207 samples used to "
        "study the microbiome related differences between "
        "males & females\n\n python meta_analysis_data.py "
        "combinedMetadata.tsv sex_contrast_meta_analysis.tsv "
        " -min age:16 \\\n --cat study_condition:control "
        "body_site:stool \\\n -mp gender:25 \\\n -cf BMI \\\n -mm gender:40 "
        " --verbose \n"
        ""
        "\nExample b)"
        "\n==> return 3385 samples used to "
        "study microbiome relationships with ageing. "
        "\n(Imagine to run this having metadata + data rel.abuns. "
        "previously stored in rel_abundances/):\n "
        "\n python meta_analysis_data.py rel_abundances/ "
        "ageing_dataset_ma.tsv \\\n -min age:16 \\\n --cat study_condition:control "
        "body_site:stool \\\n -iqr age:15 \\\n -cf BMI \\\n --debug\n\n" 
	"\n The flags used in these two commands"
        "serve for:"
        "\n -min age:16 populations will have minimum 16 year"
        "\n -cat body_site:stool: populations'll have body_site==stool"
        "\n -mp gender:25 least abun. class inside 'gender' will"
        "cover at least 25 percent of the total"
        "\n -iqr age:15: in each population, age will have an IQR of at least 15"
        "\n"
        "\nThere is no limit to each argument: e.g. you can select as mimum of maximum as you want.\n"
        "There other are parameters which are extremely usefull:\n "
        "\n--search a:b:c:d will select only samples and population having b,c or d under ciolumn a"
        "\n--minmin set a len minium required for the population ('-mm 40') or a minimum for the least abundant class ('-mm gender:40')"
        "\n"
        "--multiple. If negative (-1, default) allows just baseline samples. If an int, "
        "uses all samples inside the 'multiple' upper bound."
        "\n\n"
        "\n** IMPORTANT **\n"
        "\n\n"
        "When a flag is activated, it will also exclude the "
        "\nNaN values in the column it refers to: e.g.\n "
        "--min age:16 will exclude samples younger than 16 AND samples having age==NaN.\n"
        "Therefore: \n"
        "\n-cf/--cfd (confounder) flag simply ensure that the data are non-NaN in the column chosen.\n"
	"\n You call any argument as many times as you want BUT MULTIPLE", \
        formatter_class=ap.RawTextHelpFormatter)
    add = p.add_argument
    add("input_folder", type=str, help=\
        "Folder where the tables of metadata "
        "+ abundances are stored\n "
	"tables inside this folder will be "
        "filtered for selection \n(it can have more "
        "tables than) the ones you need.\n "
        "An alternative way is to specify "
        "under this parameter\n "
        "the combinedMetadata table, a table which "
        "contains \nall the metadata in this package.\n "
        "This last option is automatically "
        "activated when the first argument ends in "
        "'.tsv'. \n"
        "Consider, now, that if you specify this option, "
        " your dataset will be tailored \nas "
        "well for your analysis, but it will\n "
        "be constituted solely by metadata.")
    add("output_dataset", type=str, help="The name of the table with the population "
	"and relative abundances for your meta-analysis")
    #add("-z", "--featid", type=str, help="This keyword allows to distinguish features from metadata, a brief examples:"
    #	"\n 's__' is for species\n'g__' is for genera\n'PWY' is for pathways\n'K' if for KEGG-profiles")
    add("-min", "--min", type=str, nargs="+", default=[], help="Use ':' ==> age:16 means minimum 16 in the columns age")
    add("-max", "--max", type=str, nargs="+", default=[], help="Use ':' ==> BMI:45 means maximum 45 in the columns BMI")
    add("-iqr", "--iqr", type=str, nargs="+", default=[], help="Use ':' ==> age:10 means asking at least 10 as IQR for age in this population")
    add("-mul", "--multiple", type=int, default=-1, help="Upper bound of days from the baseline allowed (default -1 uses only baselines, while 0 uses all)")
    add("-cat", "--cat", type=str, nargs="+", default=[], help="Use ':' ==> study_condition:control means asking just controls")
    add("-ex", "--exclude", type=str, nargs="+", default=[], help="Works like --cat, but for exclusion (useful mixed with --search)")
    add("-search", "--search", type=str, nargs="+", default=[], help="Use ':' ==> search:study_condition:a:b:c:etc asks for a or b or c or etc in column 'study_condition'")
    add("-bn", "--binary", type=str, nargs="+", default=[], help="Use ':' (-bn cls:pos:prevpos:neg) ==> under column cls, substitute everything == prevpos with pos,  "
	"and everything else with neg (useful mixed with ex or cat)")
    add("-mp", "--min_perc", type=str, nargs="+", default=[], help="Use ':' ==> column-name:25 searches at least min 25%% in the less abundance class under column-name")
    add("-mm", "--minmin", type=str, default="0", help="Can have 2 shapes: gender:40 requires a minimum of 40 individuals in the less abundance sex "
	", while just 40 requires a minimum of 40 individuals in total (This keyword is very important because many dataset are quite small for a meta-analysis)")
    add("-cf", "--cfd", type=str, nargs="+", default=[], help="Use words, such as: BMI age gender: this will simply exclude samples **without** these confounders")
    add("-si", "--study_identifier", type=str, default="study_name", help="Name of the column identifing Dataset [default: study_name]")
    add("-v", "--verbose", action="store_true", help="Tell you what is adding and how big it is")
    add("--debug", action="store_true", help="Similar to verbose, but prints more stuff")
    return vars(p.parse_args())


def handle_input(input_argument, studyID, verbose):
    data_2_tables = {}
    if input_argument.endswith(".tsv"):
        if os.path.isfile(input_argument):
            metmet = pd.read_csv(input_argument, sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")
            for study in metmet[studyID].unique().tolist():
                data_2_tables[study] = metmet.loc[metmet[studyID]==study, :].T
        else:
            raise FileNotFoundError("Sorry: you input a combinedMetadata.tsv table string which is not present in cur dir. Exiting")
    else:
        if not os.path.isdir(  input_argument  ):
            raise FileNotFoundError("The directory you set as the input folder is not there. Exiting.")
        elif os.path.isdir(  input_argument  ) and (not len(os.listdir( input_argument  ))):
            raise IndexError( "The directory you set up as a input folder exists, but is empty. Exiting.")
        else:
            for table in glob.glob(os.path.join(input_argument, "*")):
                data_2_tables[os.path.basename(table.replace(".tsv", ""))] = pd.read_csv(table, sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")
    return data_2_tables


def select(args):

    madata = []

    import numpy

    print(args)

    metadata_dict = handle_input(args["input_folder"], args["study_identifier"], args["verbose"])

    for name in metadata_dict:
        tabtab = metadata_dict[name] #pd.read_csv(table, sep="\t", header=0, index_col=0, low_memory=False, engine="c").fillna("NA")
        tabtab.index.name = "sample_id"
        tabtab.columns = [ (dt + "_" + c) for c,dt in zip(tabtab.columns.tolist(), tabtab.loc[args["study_identifier"]].tolist()) ]
        study = tabtab.loc[args["study_identifier"]].tolist()[0]

        print(study)



        if args["min"] and (not isinstance(tabtab, np.ndarray)):
            for ag in args["min"]:
                ag = ag.split(":")
                col = ag[0]
                min_ = int(ag[1])
                if col in tabtab.index.tolist():
                    tabtab = tabtab.loc[ :, tabtab.loc[col]!="NA" ]
                    tabtab = tabtab.loc[ :, tabtab.loc[col].astype(float).astype(int)>=min_ ]
                    if args["debug"]:
                        sys.stdout.write("\nMin ==> " + str(tabtab.shape) + "[" + study + "]") 
                else:
                    tabtab = np.array([])



        if args["max"] and (not isinstance(tabtab, np.ndarray)):
            for ag in args["max"]:
                ag = ag.split(":")
                col = ag[0]
                max_ = int(ag[1])
                if col in tabtab.index.tolist():
                    tabtab = tabtab.loc[ :, tabtab.loc[col]!="NA" ]
                    tabtab = tabtab.loc[ :, tabtab.loc[col]<=max_ ]
                    if args["debug"]:
                        sys.stdout.write("\nMAX ==> " + str(tabtab.shape) + "[" + study + "]")
                else:
                    tabtab = np.array([])


            
        if args["cat"] and (not isinstance(tabtab, np.ndarray)):
            for ct in args["cat"]:
                ct = ct.split(":")
                col = ct[0]
                ctgs = ct[1:]

                if col in tabtab.index.tolist():
                    tabtab = tabtab.loc[ :, tabtab.loc[col]!="NA" ]
                    tabtab = tabtab.loc[ :, tabtab.loc[col].isin(ctgs) ]
                    if args["debug"]:
                        sys.stdout.write("\nCATG ==> " + str(tabtab.shape) + "[" + study + "]")
                else:
                    tabtab = np.array([])



        if args["exclude"] and (not isinstance(tabtab, np.ndarray)):
            for ag in args["exclude"]:
                ag = ag.split(":")
                col = ag[0]
                exc = ag[1]
                if col in tabtab.index.tolist():
                    tabtab = tabtab.loc[ :, ~tabtab.loc[col].isin([exc])]


            
        if args["cfd"] and (not isinstance(tabtab, np.ndarray)):
            for confounder in args["cfd"]:
                if (not isinstance(tabtab, np.ndarray)) and (confounder in tabtab.index.tolist()):
                    tabtab = tabtab.loc[ :, tabtab.loc[confounder]!="NA" ]
                else:
                    tabtab = np.array([])



        if args["search"] and (not isinstance(tabtab, np.ndarray)):
            for ag in args["search"]:
                col = ag.split(":")[0]
                arrgs = ag.split(":")[1:]
                if col in tabtab.index.tolist():
                    #if (not "NA" in arrgs): 
                    tabtab = tabtab.loc[ :, tabtab.loc[col]!="NA" ]
                    tabtab = tabtab.loc[ :, tabtab.loc[col].isin(arrgs) ]
                    if args["debug"]:
                        sys.stdout.write("\nSEarch ==> " + str(tabtab.shape) + "[" + study + "]")
                else:
                    tabtab = np.array([])




        if args["multiple"] and (not isinstance(tabtab, np.ndarray)):
            if ("days_from_first_collection" in tabtab.index.tolist()) and \
		(len(tabtab.loc["days_from_first_collection"].unique()) >1):

                if args["multiple"] < 0.:
                    tabtab.loc["days_from_first_collection"] = [ (0.0 if str(n)=="NA" else float(n)) \
			for n in tabtab.loc["days_from_first_collection"].tolist() ]

                    tabtab = tabtab.loc[ :, tabtab.loc["days_from_first_collection"] == 0.0 ]

                    tabtab = tabtab.T.drop_duplicates(["subject_id"], keep="first").T

                    if args["debug"]:
                        sys.stdout.write("\nDROP ==> " + str(tabtab.shape) + "[" + study + "]")
                        sys.stdout.write("\nMULT ==> " + str(tabtab.shape) + "[" + study + "]")
                    
                else:
                    tabtab = tabtab.loc[ :, tabtab.loc["days_from_first_collection"] <= args["multiple"] ]
            else:
                if args["multiple"] < 0.:
                    tabtab = tabtab.T.drop_duplicates(["subject_id"], keep="first").T
                    if args["debug"]:
                        #exit(1)
                        sys.stdout.write("\nDROP ==> " + str(tabtab.shape) + "[" + study + "]")

                #else:
                #    tabtab = np.array([])

            if not tabtab.shape[1]:
                tabtab = np.array([])

        else:
            if args["debug"]:
                sys.stdout.write("\n!! NO-DROP ==> " + str(tabtab.shape) + "[" + study + "]")

            if (not isinstance(tabtab, np.ndarray)) and (not tabtab.shape[1]):
                tabtab = np.array([])


        if args["binary"] and (not isinstance(tabtab, np.ndarray)):
            for bn in args["binary"]:
                cls, pos, prevpos, neg = bn.split(":")
                if cls in tabtab.index.tolist():
                    #if neg != "NA":
                    #    tabtab = tabtab.loc[ :, tabtab.loc[cls]!="NA" ]
                    tabtab.loc[ cls ] = [( pos if (nm==prevpos) else neg ) for nm in tabtab.loc[ cls ].tolist() ]
                else:
                    tabtab = np.array([])



        if args["min_perc"] and (not isinstance(tabtab, np.ndarray)):
            for mp in args["min_perc"]:
                col, perc = mp.split(":")
                perc = float(perc)

                if col in tabtab.index.tolist():
                    tabtab = tabtab.loc[ :, tabtab.loc[col]!="NA" ]

                    if tabtab.shape[1]:
                        tabtab = tabtab \
                                if \
                            ((np.min([((tabtab.loc[:, tabtab.loc[col]==z].shape[1]/float(tabtab.shape[1]))*100) for z in tabtab.loc[col].unique()])>=perc) \
                            and (len(tabtab.loc[col].unique())>1)) \
                                else \
                            np.array([])
                    

                    else:
                        tabtab = np.array([])

                    if args["debug"]:
                        sys.stdout.write("\nMIN PERC ==> " + str(tabtab.shape) + "[" + study + "]")
       

                else:
                    tabtab = np.array([])



        if args["iqr"] and (not isinstance(tabtab, np.ndarray)):
            for iqr in args["iqr"]:
                iqr = iqr.split(":")
                col = iqr[0]
                num = float(iqr[1])
                if (not isinstance(tabtab, np.ndarray)) and (col in tabtab.index.tolist()):
                    tabtab = tabtab.loc[ :, tabtab.loc[col]!="NA" ]
                    tabtab = tabtab if (sts.iqr(tabtab.loc[ col ].values.astype(float)) >= num) else np.array([])
                    if args["debug"]:
                        sys.stdout.write("\nIQR ==> " + str(tabtab.shape) + "[" + study + "]")
                else:
                    tabtab = np.array([])


        
        if args["minmin"] and (not isinstance(tabtab, np.ndarray)):
            if len(args["minmin"].split(":")) == 2:
                minmin = args["minmin"].split(":")
                col = minmin[0]
                mm = int(minmin[1]) 
                if tabtab.shape[1]:
                    #### print(study, [ tabtab.loc[col, tabtab.loc[col]==u].shape[0] for u in tabtab.loc[col].unique() ] )
                    tabtab = tabtab \
                        if (((np.min([ tabtab.loc[col, tabtab.loc[col]==u].shape[0] for u in tabtab.loc[col].unique()]) >= mm)) \
                            and ( len(tabtab.loc[col].unique())>1  )) \
                        else np.array([])
                else:
                    tabtab = np.array([])

                if args["debug"]:
                    sys.stdout.write("\nMIN-MIN ==> " + str(tabtab.shape) + "[" + study + "]")
            else:
                mm = int(args["minmin"])
                tabtab = tabtab if (tabtab.shape[1]>=mm) else np.array([])



        if len(tabtab):
            if args["verbose"] or args["debug"]:
                sys.stdout.write("\nMerging population %s with %i samples - " %(tabtab.loc[args["study_identifier"]].tolist()[0], tabtab.shape[1]))
            if (not len(madata)):
                madata = tabtab
            else:
                madata = madata.merge( tabtab, left_index=True, right_index=True, how="outer" )


    if not isinstance(madata, list):
        madata.fillna(0.0, inplace=True)



    if len(madata):
        if args["verbose"] or args["debug"]:
            sys.stdout.write("\nYour dataset %s has %i samples! Does this make sense?\n" %(args["output_dataset"], madata.shape[1]))
        madata.to_csv( args["output_dataset"], sep="\t", header=True, index=True )
    else:
        raise IndexError("We are sorry: apparently your search is too detailed for this package, and we can't provide any sample with the desired characteristics Exiting")


def main():
    args = read_params()
    select(args)


if __name__=="__main__": 
    main()
