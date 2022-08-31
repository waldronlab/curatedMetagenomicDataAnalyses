#!/usr/bin/env python

import sys, os
import pandas as pd
import numpy as np
import itertools as it
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from meta_analyses import paule_mandel_tau
from meta_analyses import RE_meta_binary
from meta_analyses import RE_meta

class singleStudyEffect(object):
    def __init__(self, Rho_and_P, Name, Len, REG=True):
        self.effect, self.Pvalue = Rho_and_P
        self.accepted = (self.effect==self.effect)
        self.Name = Name
        self.Len = Len if REG else (Len[0] + Len[1])
        if REG:
            self.ncases = None
            self.ncontrols = None
        else:
            self.ncases = Len[1]
            self.ncontrols = Len[0]

class meta_analysis_with_linear_model(object):
    def __init__(self, dataTable, formula, studyid, feats, outfile, cls_or_reg, heterogeneity, pos=None, neg=None):

        #if not (("study_name" in dataTable.columns) or ("dataset_name" in dataTable.columns)):
        self.data = dataTable.T
        #else:
        #    self.data = dataTable

        self.feats = feats
        self.studyid = studyid
        self.formula = formula
        fm_covs = formula.split("+")
        self.predictor = fm_covs[0].strip()
        self.pos, self.neg = pos, neg
        self.covariates = [ c.strip() for c in fm_covs[1:] ]

        self.easy_names_covariates = dict([(c, c[2:-1]) for c in self.covariates if c.startswith("C(") ])

        for cv in self.covariates:
       
            if not cv.startswith("C"):
    
                print("NA" in self.data[cv].tolist())

                self.data[cv] = self.data[cv].values.astype(float)
            #print(cv)

            #if cv.startswith("log"):
            #    self.data[ cv.split("(")[1][:-1] ] = self.data[  cv.split("(")[1][:-1]  ].astype(float)   
 
        if self.neg:
            self.formula = formula.replace(self.predictor, "C(%s, Treatment('%s'))" %(self.predictor, self.neg))

        self.outfile = outfile
        self.het = heterogeneity
        self.studies = list(self.data[self.studyid].unique())
        self.N_studies = dict([(study, self.data.loc[self.data[self.studyid].isin([study]), :].shape[0]) for study in self.studies])
        self.cls_or_reg = cls_or_reg
 
        if (self.cls_or_reg in ["REG", "LOGC"]) or ((self.cls_or_reg == "UNS") and (not self.neg) and (not self.pos)):
            self.data[self.predictor] = self.data[self.predictor].values.astype(float)

        #if self.cls_or_reg == "LOG":
            



    def correlation_of_study(self, study, feature):
 
        ##print(feature, self.predictor, self.covariates, self.studyid)
        
        data_here = self.data.loc[self.data[self.studyid].isin([study]), [feature, self.predictor] \
            + [self.easy_names_covariates.get(c, c) for c in self.covariates]]
        data_here[feature] = data_here[feature].astype(float)

        formula = ('Q("%s") ~ ' %feature) + self.formula
        
        if self.cls_or_reg in ["LOG", "LOGC"]:
            #formula = ('%s ~ ' %feature) + self.formula
            data_here[feature] = data_here[feature].astype(float).astype(bool).astype(int)

        #print(formula)
        #print(data_here)

        if not self.cls_or_reg in ["LOG", "LOGC"]:
            md = smf.ols(formula, data=data_here)
            model_fit = md.fit()

        else:
            #print(data_here[feature]) ## = data_here[feature].astype(bool).astype(int)

            md = smf.logit(formula, data=data_here)
            model_fit = md.fit()

        #### print(model_fit.summary())

        if self.cls_or_reg == "UNS":
            if self.neg:
                b = model_fit.params.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos) ]
                p = model_fit.pvalues.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos) ]
                bse = model_fit.bse.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos) ]

            else:
                b = model_fit.params.loc[ "%s" %self.predictor ]
                p = model_fit.pvalues.loc[ "%s" %self.predictor ]
                bse = model_fit.bse.loc[ "%s" %self.predictor ]

            return b, bse, p, study

        elif self.cls_or_reg == "CLS":
            #predictor = self.predictor if (not "Treatment" in self.predictor) else (self.predictor.split(",")[0].replace("C(", "")) + ("[T.%s]" %self.pos)
            t = model_fit.tvalues.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos)]
            n1 = float(len(data_here.loc[(data_here[self.predictor]==self.neg)]))
            n2 = float(len(data_here.loc[(data_here[self.predictor]==self.pos)]))
            d = (t*(n1+n2))/float(np.sqrt(n1*n2)*np.sqrt(n1+n2-2))
            SEd = np.sqrt(((n1+n2-1)/float(n1+n2-3)) * ((4./float(n1+n2))*(1+((d**2.)/8.))))
            d_lw = d-(1.96*SEd)
            d_up = d+(1.96*SEd)
            return d, model_fit.pvalues.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos) ], (n1,n2), False

        elif self.cls_or_reg == "REG":
            t = model_fit.tvalues.loc[self.predictor]
            n = float(len(data_here))
            r = float(t) / np.sqrt(np.float((t**2.) + (n - 1.)))
            Zr = np.arctanh(r)
            SEr = 1/np.sqrt(n - 3)
            r_lw = Zr - (1.96*SEr)
            r_up = Zr + (1.96*SEr)
            return np.tanh(r), model_fit.pvalues.loc[self.predictor]

        elif self.cls_or_reg == "LOG":
            n1 = float(len(data_here.loc[(data_here[self.predictor]==self.neg)]))
            n2 = float(len(data_here.loc[(data_here[self.predictor]==self.pos)]))
            return model_fit.params.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos) ], \
		model_fit.pvalues.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos) ], \
		model_fit.bse.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos) ], (n1,n2)

        elif self.cls_or_reg == "LOGC":
            n = float(len(data_here))
            return model_fit.params.loc[ "%s" %(self.predictor) ], \
                model_fit.pvalues.loc[ "%s" %(self.predictor) ], \
                model_fit.bse.loc[ "%s" %(self.predictor) ], n
            




    def random_effect_regression_model(self):
        result = []

        for feature in self.feats:
            study2lens = {}
            singleStudiesClass = []

            if self.cls_or_reg == "REG":

                 singleStudiesRegs = [singleStudyEffect(self.correlation_of_study(study, feature), study, self.N_studies[study]) for study in self.studies]
                 considerable = sorted([e for e in singleStudiesRegs if e.accepted], key=lambda x : x.effect)

                 if len(considerable)>=(int(len(self.studies)//4) if (len(self.studies)>=5) else len(self.studies)):
                     re = RE_meta(\
                         [e.effect for e in considerable], \
                         [e.Pvalue for e in considerable], \
                         [e.Name for e in considerable], \
                         [e.Len for e in considerable], \
                         feature, het=self.het, REG=True)

                     sys.stdout.write("%s (%i studies) Random Effect Effect = %.3f [H: %.3f %.3f  %.3f]\n" \
                         %(feature, len(considerable), re.RE, re.t2_DL, re.t2_PM, re.I2))

                     if len(result) == 0: result = re.result
                     else: result = result.append(re.result)

            #elif self.cls_or_reg == "LOGC":
                

            elif self.cls_or_reg == "CLS":

                for study in self.studies:
                    cohenD, Pvalue, Lens, variance = self.correlation_of_study(study, feature)

                    singleStudiesClass += [singleStudyEffect((cohenD, Pvalue), study, Lens, False)]
                    study2lens[study] = {}
                    study2lens[study]["cases"] = Lens[1]
                    study2lens[study]["control"] = Lens[0]

                considerable = [e for e in singleStudiesClass if (e.accepted)]
                if len(considerable)>=(int(len(self.studies)//4) if (len(self.studies)>=5) else len(self.studies)):

                    re = RE_meta_binary(\
                        [e.effect for e in considerable], \
                        [e.Pvalue for e in considerable], \
                        [e.Name for e in considerable], \
                        [study2lens[e.Name]["cases"] for e in considerable], \
                        [study2lens[e.Name]["control"] for e in considerable], \
                        feature, "D", False, self.het)

                    sys.stdout.write("%s (%i studies) Random Effect Effect = %.3f [H: %.3f %.3f  %.3f]\n" \
                        %(feature, len(considerable), re.RE, re.t2_DL, re.t2_PM, re.I2))

                    if len(result) == 0:
                        result = re.result
                    else:
                        result = result.append(re.result)

            elif self.cls_or_reg in ["LOG", "LOGC"]:
                variances_considered = []

                for study in self.studies:
                    try:
                        LOR, Pvalue, variance, Lens = self.correlation_of_study(study, feature)
                        if self.cls_or_reg == "LOGC":
                            Lens = (Lens, 0)
                    except statsmodels.tools.sm_exceptions.PerfectSeparationError:
                        LOR, Pvalue, variance, Lens = np.nan, 1.0, np.nan, (0, 0)
                    except np.linalg.LinAlgError:
                        LOR, Pvalue, variance, Lens = np.nan, 1.0, np.nan, (0, 0)

                    if Pvalue != Pvalue:
                        LOR, Pvalue, variance, Lens = np.nan, 1.0, np.nan, (0, 0)

                    variances_considered += [ variance ]
                    
                    singleStudiesClass += [singleStudyEffect((LOR, Pvalue), study, Lens, False)]
                    study2lens[study] = {}
                    study2lens[study]["cases"] = Lens[1]
                    study2lens[study]["control"] = Lens[0]

                considerable = [(e) for e,v in zip(singleStudiesClass, variances_considered) if ((e.accepted) and (v**2. > 0.))]
                variances = [ (v**2.) for v,e in zip(variances_considered, singleStudiesClass) if ((e.accepted) and (v**2. > 0.))]

                if len(considerable)>=(int(len(self.studies)//4) if (len(self.studies)>=5) else len(self.studies)):

                    re = RE_meta_binary(\
                        [e.effect for e in considerable], \
                        [e.Pvalue for e in considerable], \
                        [e.Name for e in considerable], \
                        [study2lens[e.Name]["cases"] for e in considerable], \
                        [study2lens[e.Name]["control"] for e in considerable], \
                        feature, "precomputed", variances, self.het)

                    sys.stdout.write("%s (%i studies) Random Effect Effect = %.3f [H: %.3f %.3f  %.3f]\n" \
                        %(feature, len(considerable), re.RE, re.t2_DL, re.t2_PM, re.I2))

                    if len(result) == 0:
                        result = re.result
                    else:
                        result = result.append(re.result)

            elif self.cls_or_reg in ["UNS"]:
                variances_considered = []
                considerable = []

                for study in self.studies:
                    beta, betase, betapval, name = self.correlation_of_study(study, feature)

                    if betapval != betapval:
                        beta, betase, betapval, name = 0.0, 0.0, 1.0, "NULL"

                    considerable += [(beta, betase, betapval, name)]

                considerable = [e for e in considerable if (e[2] == e[2])]

                if len(considerable)>=(int(len(self.studies)//4) if (len(self.studies)>=5) else len(self.studies)):

                    variances = [(e[1]**2.) for e in considerable]
                    re = RE_meta_binary(\
                        [e[0] for e in considerable], \
                        [e[2] for e in considerable], \
                        [e[3] for e in considerable], \
                        [None for e in considerable], \
                        [None for e in considerable], \
                        feature, "precomputed", variances, self.het)

                    sys.stdout.write("%s (%i studies) Random Effect Effect = %.3f [H: %.3f %.3f  %.3f]\n" \
                        %(feature, len(considerable), re.RE, re.t2_DL, re.t2_PM, re.I2))

                    if len(result) == 0:
                        result = re.result
                    else:
                        result = result.append(re.result)

        _, FDR = fdrcorrection(result["RE_Pvalue"].values.astype(float), alpha=.05)
        
        result.insert((len(self.studies))+1, "RE_%s_Qvalue" %("Effect" if (self.cls_or_reg in ["CLS", "LOG", "LOGC", "UNS"]) else "Correlation"), FDR)
        result.fillna("NA", inplace=True)


        for S,_p_dist_col in zip(self.studies, [(s+"_Pvalue") for s in self.studies]):

            if _p_dist_col in set(result.columns.tolist()):
                feat_p_map = [[f,p] for f,p in zip(result.index.tolist(), result[_p_dist_col].tolist()) if p!="NA"]
                ps = [p[1] for p in feat_p_map]
                _,fdr = fdrcorrection(ps, alpha=0.05)
                littleFrame = pd.DataFrame({ S+"_Qvalue": fdr }, index=[f[0] for f in feat_p_map])
                result = result.join(littleFrame)
                result.fillna("NA", inplace=True)
            else:
                result[_p_dist_col] = "NA"
                result[S+"_Correlation"] = "NA"
                result[S+"_Qvalue"] = "NA"

        result.insert(0, "Feature", result.index.tolist())


        try:
            result = result[["Feature"] + \
                list(it.chain.from_iterable([[S+"_Correlation",  S+"_SE", S+"_Pvalue", S+"_Qvalue"] for S in self.studies])) + \
                ["RE_Correlation", "RE_Pvalue", "RE_Correlation_Qvalue", "RE_stdErr", "RE_conf_int", "Zscore", "Tau2_DL", "Tau2_PM", "I2"] \
                if (self.cls_or_reg=="REG") else \
                list(it.chain.from_iterable([[S+"_Effect",  S+"_SE", S+"_Pvalue", S+"_Qvalue"] for S in self.studies])) + \
                ["RE_Effect", "RE_Var", "RE_Pvalue", "RE_Effect_Qvalue", "RE_stdErr", "RE_conf_int", "Zscore", "Tau2_DL", "Tau2_PM", "I2"]]

        except KeyError:
            studies_found = [s for s in self.studies if ((s+"_SE") in result.columns)]
            result = result[["Feature"] + \
                list(it.chain.from_iterable([[S+"_Correlation",  S+"_SE", S+"_Pvalue", S+"_Qvalue"] for S in studies_found])) + \
                ["RE_Correlation", "RE_Pvalue", "RE_Correlation_Qvalue", "RE_stdErr", "RE_conf_int", "Zscore", "Tau2_DL", "Tau2_PM", "I2"] \
                if (self.cls_or_reg=="REG") else \
                list(it.chain.from_iterable([[S+"_Effect",  S+"_SE", S+"_Pvalue", S+"_Qvalue"] for S in studies_found])) + \
                ["RE_Effect", "RE_Var", "RE_Pvalue", "RE_Effect_Qvalue", "RE_stdErr", "RE_conf_int", "Zscore", "Tau2_DL", "Tau2_PM", "I2"]]
 
        #else:
        #    result = result[["Feature"] + \
        #        list(it.chain.from_iterable([[S+"_LOR",  S+"_SE", S+"_Pvalue", S+"_Qvalue"] for S in self.studies])) + \
        #       ["RE_LOR", "RE_Pvalue", "RE_LOR_Qvalue", "RE_stdErr", "RE_conf_int", "Zscore", "Tau2_DL", "Tau2_PM", "I2"]]

        

        result.to_csv(self.outfile, sep="\t", header=True, index=True)



class continuous_analysis_with_logistic_model(object):
    def __init(self, datatable, formula, feats, outfile, pos=None, neg=None):
        self.data = dataTable.T
        self.feats = feats
        fm_covs = formula.split("+")
        self.predictor = fm_covs[0].strip()
        self.pos, self.neg = pos, neg
        self.covariates = [c.strip() for c in fm_covs[1:]]

        ### self.easy_names_covariates = dict([(c, c[2:-1]) for c in self.covariates if c.startswith("C(") ])

        self.data[self.predictor] = self.data[self.predictor].values.astype(float)

        for cv in self.covariates:
            if not cv.startswith("C"):
                self.data[cv] = self.data[cv].astype(float)

        if self.neg:
            self.formula = formula.replace(self.predictor, "C(%s, Treatment('%s'))" %(self.predictor, self.neg))
        else:
            self.formula = formula

        self.outfile = outfile

        effects = []
        pvals = []
        Vars = []
        Stderrs = []

        for ft in self.feats:
            try:
                effect, StdErr, p_value, n1, n2 = self.logistic_of_study(ft)
                if p_value != p_value:
                    effect, StdErr, p_value, n1, n2 = np.nan, np.nan, 1.0, np.nan, np.nan
            except statsmodels.tools.sm_exceptions.PerfectSeparationError:
                effect, StdErr, p_value, n1, n2 = np.nan, np.nan, 1.0, np.nan, np.nan
            except np.linalg.LinAlgError:
                effect, StdErr, p_value, n1, n2 = np.nan, np.nan, 1.0, np.nan, np.nan

            effects += [effect]
            pvals += [p_value]
            Stderrs += [StdErr]
            Vars += [StdErr**2.]

        self.frame = pd.DataFrame(\
                {"RE_Effect": effects, "RE_Pvalue": pvals, "RE_stdErr": Stderrs, "RE_Var": Vars}, index=self.feats)

        _, FDR = fdrcorrection(np.nan_to_num( np.array(pvals, dtype=np.float64), nan=1.0), alpha=0.05)

        self.frame["Qvalue"] = FDR
        self.frame.index.name = "Feature"
        for ft in self.frame.index.tolist():
            print(ft, self.frame.loc[ft, "RE_Effect"], self.frame.loc[ft,"Qvalue"])

    def write_out(self):
        self.frame.to_csv(self.outfile, sep="\t", header=True, index=True)

    def logistic_of_study(self, feature):
        self.data[feature] = self.data[feature].astype(float).astype(bool).astype(int)

        #formula = ('Q("%s") ~ ' %feature) + self.formula
        formula = ('Q("%s") ~ ' %feature) + self.formula

        md = smf.logit(formula, data=self.data)
        model_fit = md.fit()
 
        ##predictor = self.predictor if (not "Treatment" in self.predictor) else (self.predictor.split(",")[0].replace("C(", "")) + ("[T.%s]" %self.pos)

        #t = model_fit.tvalues.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos)]
        #d = (t*(n1+n2))/float(np.sqrt(n1*n2)*np.sqrt(n1+n2-2))
        #SEd = np.sqrt(((n1+n2-1)/float(n1+n2-3)) * ((4./float(n1+n2))*(1+((d**2.)/8.))))
        #d_lw = d-(1.96*SEd)
        #d_up = d+(1.96*SEd)
        return model_fit.params.loc[ "%s" %(self.predictor)  ], \
                model_fit.bse.loc[ "%s" %(self.predictor)  ], \
                model_fit.pvalues.loc[ "%s" %(self.predictor)  ], L


class analysis_with_logistic_model(object):
    def __init__(self, dataTable, formula, feats, outfile, pos, neg):

        self.data = dataTable.T
        self.feats = feats
        fm_covs = formula.split("+")
        self.predictor = fm_covs[0].strip()
        self.pos, self.neg = pos, neg
        self.covariates = [c.strip() for c in fm_covs[1:]]
        ### self.easy_names_covariates = dict([(c, c[2:-1]) for c in self.covariates if c.startswith("C(") ])

        for cv in self.covariates:
            if not cv.startswith("C"):
                self.data[cv] = self.data[cv].astype(float)

        if self.neg:
            self.formula = formula.replace(self.predictor, "C(%s, Treatment('%s'))" %(self.predictor, self.neg))
        else:
            self.formula = formula #+ " + 1 "

        self.outfile = outfile

        ## main
        effects = []
        pvals = []
        Vars = []
        Stderrs = []

        for ft in self.feats:
            try:
                effect, StdErr, p_value, n1, n2 = self.logistic_of_study(ft)
                if p_value != p_value:
                    effect, StdErr, p_value, n1, n2 = np.nan, np.nan, 1.0, np.nan, np.nan
            except statsmodels.tools.sm_exceptions.PerfectSeparationError:
                effect, StdErr, p_value, n1, n2 = np.nan, np.nan, 1.0, np.nan, np.nan
            except np.linalg.LinAlgError:
                effect, StdErr, p_value, n1, n2 = np.nan, np.nan, 1.0, np.nan, np.nan

            effects += [effect]
            pvals += [p_value]
            Stderrs += [StdErr]
            Vars += [StdErr**2.]
          
        self.frame = pd.DataFrame(\
                {"RE_Effect": effects, "RE_Pvalue": pvals, "RE_stdErr": Stderrs, "RE_Var": Vars}, index=self.feats)

        _, FDR = fdrcorrection(np.nan_to_num( np.array(pvals, dtype=np.float64), nan=1.0), alpha=0.05)

        self.frame["Qvalue"] = FDR
        self.frame.index.name = "Feature"
        for ft in self.frame.index.tolist():
            print(ft, self.frame.loc[ft, "RE_Effect"], self.frame.loc[ft,"Qvalue"])

    def write_out(self):
        self.frame.to_csv(self.outfile, sep="\t", header=True, index=True)

    def logistic_of_study(self, feature):
        self.data[feature] = self.data[feature].astype(float).astype(bool).astype(int)

        #formula = ('Q("%s") ~ ' %feature) + self.formula
        formula = ('Q("%s") ~ ' %feature) + self.formula

        md = smf.logit(formula, data=self.data)
        model_fit = md.fit()
 
        print(model_fit.summary())
        ##predictor = self.predictor if (not "Treatment" in self.predictor) else (self.predictor.split(",")[0].replace("C(", "")) + ("[T.%s]" %self.pos)

        #t = model_fit.tvalues.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos)]
        n1 = float(len(self.data.loc[(self.data[ self.predictor ]==self.neg)]))
        n2 = float(len(self.data.loc[(self.data[ self.predictor ]==self.pos)]))
        #d = (t*(n1+n2))/float(np.sqrt(n1*n2)*np.sqrt(n1+n2-2))
        #SEd = np.sqrt(((n1+n2-1)/float(n1+n2-3)) * ((4./float(n1+n2))*(1+((d**2.)/8.))))
        #d_lw = d-(1.96*SEd)
        #d_up = d+(1.96*SEd)
        return model_fit.params.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos)  ], \
		model_fit.bse.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos)  ], \
		model_fit.pvalues.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos)  ], n1, n2



class analysis_with_linear_model(object):
    def __init__(self, dataTable, formula, feats, outfile, pos, neg, standardised=True):

        self.data = dataTable.T
        self.feats = feats
        fm_covs = formula.split("+")
        self.predictor = fm_covs[0].strip()
        self.pos, self.neg = pos, neg
        self.covariates = [c.strip() for c in fm_covs[1:]]
        ### self.easy_names_covariates = dict([(c, c[2:-1]) for c in self.covariates if c.startswith("C(") ])

        for cv in self.covariates:
            if not cv.startswith("C"):
                self.data[cv] = self.data[cv].astype(float)

        if self.neg:
            self.formula = formula.replace(self.predictor, "C(%s, Treatment('%s'))" %(self.predictor, self.neg))
        else:
            self.formula = formula #+ " + 1 "

        self.outfile = outfile
        
        ## main
        effects = []
        pvals = []
        Vars = []
        Stderrs = []

        for ft in self.feats:
            if standardised:
                effect, StdErr, p_value, n1, n2 = self.correlation_of_study(ft)
            else:
                effect, StdErr, p_value = self.unstandardised_correlation_of_study(ft)

            effects += [effect]
            pvals += [p_value]
            Stderrs += [StdErr]
            Vars += [StdErr**2.]
  
        self.frame = pd.DataFrame(\
                {"RE_Effect": effects, "RE_Pvalue": pvals, "RE_stdErr": Stderrs, "RE_Var": Vars}, index=self.feats)

        _, FDR = fdrcorrection(np.nan_to_num( np.array(pvals, dtype=np.float64), nan=1.0), alpha=0.05)

        self.frame["Qvalue"] = FDR
        self.frame.index.name = "Feature"
        for ft in self.frame.index.tolist():
            print(ft, self.frame.loc[ft, "RE_Effect"], self.frame.loc[ft,"Qvalue"])

    def write_out(self):
        self.frame.to_csv(self.outfile, sep="\t", header=True, index=True)

    def unstandardised_correlation_of_study(self, feature):
        self.data[feature] = self.data[feature].astype(float)
        formula = ('Q("%s") ~ ' %feature) + self.formula
        md = smf.ols(formula, data=self.data)
        model_fit = md.fit()

        if self.neg:
            b = model_fit.params.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos) ]
            p = model_fit.pvalues.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos) ]
            bse = model_fit.bse.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos) ]

        else:
            b = model_fit.params.loc[ "%s" %self.predictor ]
            p = model_fit.pvalues.loc[ "%s" %self.predictor ]
            bse = model_fit.bse.loc[ "%s" %self.predictor ]

        return b, bse, p

    def correlation_of_study(self, feature):
        self.data[feature] = self.data[feature].astype(float)
        formula = ('Q("%s") ~ ' %feature) + self.formula
        md = smf.ols(formula, data=self.data)
        model_fit = md.fit()

        ##print(model_fit.summary())
        ##predictor = self.predictor if (not "Treatment" in self.predictor) else (self.predictor.split(",")[0].replace("C(", "")) + ("[T.%s]" %self.pos)

        t = model_fit.tvalues.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos)]
        n1 = float(len(self.data.loc[(self.data[ self.predictor ]==self.neg)]))
        n2 = float(len(self.data.loc[(self.data[ self.predictor ]==self.pos)]))
        d = (t*(n1+n2))/float(np.sqrt(n1*n2)*np.sqrt(n1+n2-2))
        SEd = np.sqrt(((n1+n2-1)/float(n1+n2-3)) * ((4./float(n1+n2))*(1+((d**2.)/8.))))
        d_lw = d-(1.96*SEd)
        d_up = d+(1.96*SEd)
        return d, SEd, model_fit.pvalues.loc[ "C(%s, Treatment('%s'))[T.%s]" %(self.predictor, self.neg, self.pos)  ], n1, n2
