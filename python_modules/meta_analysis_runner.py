#!/usr/bin/env python

import sys, os
import pandas as pd
import numpy as np
import itertools as it
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
        fm_covs = formula.split("+")
        self.predictor = fm_covs[0].strip()
        self.pos, self.neg = pos, neg
        self.covariates = [ c.strip() for c in fm_covs[1:]]

        for cv in self.covariates:
            if not cv.startswith("Q"):
                self.data[cv] = self.data[cv].astype(float)
        
        if self.neg:
            self.formula = formula.replace(self.predictor, "C(%s), Treatment('%s')" %(self.predictor, self.neg))

        self.outfile = outfile
        self.het = heterogeneity
        self.studies = list(self.data[self.studyid].unique())
        self.N_studies = dict([(study, self.data.loc[self.data[self.studyid].isin([study]), :].shape[0]) for study in self.studies])
        self.cls_or_reg = cls_or_reg


    def correlation_of_study(self, study, feature):

        #print(self.data, " quesrui sono i miei dati")
        #print(study in self.data["dataset_name"].tolist())
        #print(feature, self.predictor, self.covariates, self.studyid)

        data_here = self.data.loc[self.data[self.studyid].isin([study]), [feature, self.predictor] + self.covariates]
        data_here[feature] = data_here[feature].astype(float)
        formula = ('Q("%s") ~ ' %feature) + self.formula

        md = smf.ols(formula, data=data_here)
        model_fit = md.fit()
        #print(model_fit.summary())
        #print(self.predictor , self.pos, self.neg, "self ghye")

        if self.cls_or_reg == "CLS":
            #predictor = self.predictor if (not "Treatment" in self.predictor) else (self.predictor.split(",")[0].replace("C(", "")) + ("[T.%s]" %self.pos)
            t = model_fit.tvalues.loc[self.predictor + ("[T.%s]" %self.pos)]
            n1 = float(len(data_here.loc[(data_here[self.predictor]==self.neg)]))
            n2 = float(len(data_here.loc[(data_here[self.predictor]==self.pos)]))
            d = (t*(n1+n2))/float(np.sqrt(n1*n2)*np.sqrt(n1+n2-2))
            SEd = np.sqrt(((n1+n2-1)/float(n1+n2-3)) * ((4./float(n1+n2))*(1+((d**2.)/8.))))
            d_lw = d-(1.96*SEd)
            d_up = d+(1.96*SEd)
            return d, model_fit.pvalues.loc[self.predictor + ("[T.%s]" %self.pos)], (n1,n2), False

        elif self.cls_or_reg == "REG":
            t = model_fit.tvalues.loc[self.predictor]
            n = float(len(data_here))
            r = float(t) / np.sqrt(np.float((t**2.) + (n - 1.)))
            Zr = np.arctanh(r)
            SEr = 1/np.sqrt(n - 3)
            r_lw = Zr - (1.96*SEr)
            r_up = Zr + (1.96*SEr)
            return np.tanh(r), model_fit.pvalues.loc[self.predictor]


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
 
        _, FDR = fdrcorrection(result["RE_Pvalue"].values.astype(float), alpha=.05)
        result.insert((len(self.studies)*2)+1, "RE_%s_Qvalue" %("Effect" if (self.cls_or_reg == "CLS") else "Correlation"), FDR)
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
        result = result[["Feature"] + \
            list(it.chain.from_iterable([[S+"_Correlation",  S+"_Pvalue", S+"_Qvalue"] for S in self.studies])) + \
            ["RE_Correlation", "RE_Pvalue", "RE_Correlation_Qvalue", "RE_stdErr", "RE_conf_int", "Zscore", "Tau2_DL", "Tau2_PM", "I2"] \
        if (self.cls_or_reg=="REG") else \
            list(it.chain.from_iterable([[S+"_Effect",  S+"_Pvalue", S+"_Qvalue"] for S in self.studies])) + \
            ["RE_Effect", "RE_Var", "RE_Pvalue", "RE_Effect_Qvalue", "RE_stdErr", "RE_conf_int", "Zscore", "Tau2_DL", "Tau2_PM", "I2"]]
        result.to_csv(self.outfile, sep="\t", header=True, index=True)


class analysis_with_linear_model(object):
    def __init__(self, dataTable, formula, feats, outfile, pos, neg):

        self.data = dataTable.T
        self.feats = feats
        fm_covs = formula.split("+")
        self.predictor = fm_covs[0].strip()
        self.pos, self.neg = pos, neg
        self.covariates = [c.strip() for c in fm_covs[1:]]

        for cv in self.covariates:
            if not cv.startswith("Q"):
                self.data[cv] = self.data[cv].astype(float)

        self.formula = formula #+ " + 1 "
        self.outfile = outfile
        
        ## main
        effects = []
        pvals = []
        Vars = []
        for ft in self.feats:
            effect, StdErr, p_value, n1, n2 = self.correlation_of_study(ft)
            effects += [effect]
            pvals += [p_value]
            Vars += [StdErr**2.]
  
        self.frame = pd.DataFrame({"RE_Effect": effects, "RE_Pvalue": pvals, "RE_Var": Vars}, index=self.feats)
        _, FDR = fdrcorrection(np.array(pvals, dtype=np.float64), alpha=0.05)

        self.frame["Qvalue"] = FDR
        self.frame.index.name = "Feature"
        for ft in self.frame.index.tolist():
            print(ft, self.frame.loc[ft, "RE_Effect"], self.frame.loc[ft,"Qvalue"])

    def write_out(self):
        self.frame.to_csv(self.outfile, sep="\t", header=True, index=True)

    def correlation_of_study(self, feature):
        self.data[feature] = self.data[feature].astype(float)
        formula = ('Q("%s") ~ ' %feature) + self.formula
        md = smf.ols(formula, data=self.data)
        model_fit = md.fit()
        print(model_fit.summary())

        predictor = self.predictor if (not "Treatment" in self.predictor) else (self.predictor.split(",")[0].replace("C(", "")) + ("[T.%s]" %self.pos)
        t = model_fit.tvalues.loc[predictor]
        n1 = float(len(data_here.loc[(data_here[predictor.replace(("[T.%s]" %self.pos),"")]==self.neg)]))
        n2 = float(len(data_here.loc[(data_here[predictor.replace(("[T.%s]" %self.pos),"")]==self.pos)]))
        d = (t*(n1+n2))/float(np.sqrt(n1*n2)*np.sqrt(n1+n2-2))
        SEd = np.sqrt(((n1+n2-1)/float(n1+n2-3)) * ((4./float(n1+n2))*(1+((d**2.)/8.))))
        d_lw = d-(1.96*SEd)
        d_up = d+(1.96*SEd)
        return d, SEd, model_fit.pvalues.loc[predictor], n1, n2
