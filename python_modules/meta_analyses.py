#!/usr/bin/env python

import sys, os
import pandas as pd
import numpy as np
from scipy import stats as sts
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection

def paule_mandel_tau(eff, var_eff, tau2_start=0, atol=1e-5, maxiter=50):
    tau2 = tau2_start
    k = eff.shape[0]
    converged = False
    for i in range(maxiter):
        w = 1 / (var_eff + tau2)
        m = w.dot(eff) / w.sum(0)
        resid_sq = (eff - m)**2
        q_w = w.dot(resid_sq)
        # estimating equation
        ee = q_w - (k - 1)
        if ee < 0:
            tau2 = 0
            converged = 0
            break
        if np.allclose(ee, 0, atol=atol):
            converged = True
            break
        # update tau2
        delta = ee / (w**2).dot(resid_sq)
        tau2 += delta
    return tau2, converged

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

## class for meta-analysis (with cohen d)
class RE_meta_binary(object):
    def __init__(self, effects, Pvalues, studies, n_cases, n_controls, \
	    responseName, EFF="D", variances_from_outside=False, CI=False, HET="PM"):

        self.responseName = responseName
        self.studies = studies
        self.n_cases = n_cases
        self.n_controls = n_controls
        self.singleStudyPvalues = Pvalues
        self.effects = np.array(effects, dtype=np.float64)
        self.n = float(len(studies))
        self.HET = HET

        if (EFF.lower() != "d") and (EFF.lower() != "precomputed"):
            raise NotImplementedError("Sorry: this script works just with EFF=D or precomputed")

        if EFF.lower() != "precomputed":
            self.n_studies = [(a+b) for a,b in zip(n_cases,n_controls)]
        else:
            self.n_studies = "NULL"

        if EFF.lower() == "precomputed":
            self.vi = np.array(variances_from_outside, dtype=np.float64)
            if not CI:
                self.devs = np.sqrt( self.vi )
                self.CI_of_d = [ (d-( 1.96*dv ), d+( 1.96*dv ))  for d,dv in zip(self.effects, self.devs)]
            else:
                self.CI_of_d = CI
        else:
            self.n_studies = "NULL"
            self.vi = np.array([(((nc+nt-1)/float(nc+nt-3)) * ((4./float(nc+nt))*(1+((eff**2.)/8.)))) for nt,nc,eff in \
                zip(n_cases,n_controls,effects)], dtype=np.float64)
            if not CI:
                self.devs = np.sqrt(self.vi)
                self.CI_of_d = []
                for d,n1,n2 in zip(self.effects, n_controls, n_cases):
                    SEd = np.sqrt(((n1+n2-1)/float(n1+n2-3)) * ((4./float(n1+n2))*(1+((d**2.)/8.))))
                    d_lw = d-(1.96*SEd)
                    d_up = d+(1.96*SEd)
                    self.CI_of_d += [[d_lw, d_up]]
            else:
                self.CI_of_d = CI
        
        self.w = np.array([(1./float(v)) for v in self.vi], dtype=np.float64)
        mu_bar = np.sum(a*b for a,b in zip(self.w, self.effects))/np.sum(self.w)
        self.Q = np.sum(a*b for a,b in zip(self.w, [(x - mu_bar)**2 for x in self.effects]))
        H = np.sqrt(self.Q/(self.n - 1))
        self.I2 = np.max([0., (self.Q-(len(self.vi)-1))/float(self.Q)])
        self.t2_PM, self.t2PM_conv = paule_mandel_tau(self.effects, self.vi)
        self.t2_DL = ((self.Q - self.n + 1) / self.scaling( self.w )) if (self.Q > (self.n-1)) else 0.

        if self.HET == "PM":
            self.W = [(1./float(v+self.t2_PM)) for v in self.vi]
        elif self.HET.startswith("FIX"):
            self.W = [(1./float(v)) for v in self.vi]
        else:
            self.W = [(1./float(v+self.t2_DL)) for v in self.vi]

        self.RE = self.CombinedEffect()
        self.stdErr = self.StdErrCombinedEffect(self.CombinedEffectVar())
        self.Zscore = self.CombinedEffectZScore(self.RE, self.stdErr)
        self.REvar = self.CombinedEffectVar()
        self.Pval = self.Pvalue(self.Zscore)
        self.conf_int = self.CombinedEffectConfInt(self.RE, self.stdErr)
        self.result = self.nice_shape()

    def tot_var(self, Effects, Weights):
        Q = np.sum(Weights * [x**2 for x in Effects]) - ((np.sum(Weights*Effects)**2)/np.sum(Weights))
        return Q

    def scaling(self, W):
        C = np.sum(W) - (np.sum([w**2 for w in W])/float(np.sum(W)))
        return C

    def tau_squared_DL(self, Q, df, C):
        return (Q-df)/float(C) if (Q>df) else 0.

    def CombinedEffect(self):
        return np.sum(self.W*self.effects)/float(np.sum(self.W))

    def CombinedEffectVar(self):
        return 1/float(np.sum(self.W))

    def StdErrCombinedEffect(self, CVar):
        return np.sqrt(CVar)

    def CombinedEffectConfInt(self, CE, SE):
        low = CE - 1.96*SE
        upp = CE + 1.96*SE
        return low, upp

    def CombinedEffectZScore(self, CE, SE):
        return CE / float(SE)

    def Pvalue(self, Z):
        return 2.*(1 - sts.norm.cdf(np.abs(Z)))

    def nice_shape(self):
        NS = {}
        for eff,P,study in zip(self.effects, self.singleStudyPvalues, self.studies):
            NS[str(study) + "_Effect"] = eff
            NS[str(study) + "_Pvalue"] = P
        NS["RE_Effect"] = self.RE
        NS["RE_Pvalue"] = self.Pval
        NS["RE_stdErr"] = self.stdErr
        NS["RE_conf_int"] = ";".join(list(map(str,self.conf_int)))
        NS["RE_Var"] = self.REvar
        NS["Zscore"] = self.Zscore
        NS["Tau2_DL"] = self.t2_DL
        NS["Tau2_PM"] = self.t2_PM
        NS["I2"] = self.I2
        NS = pd.DataFrame(NS, index=[self.responseName])
        return NS

## Class for meta-regression (with Fisher-Z included)
class RE_meta(object):
    def __init__(self, effects, Pvalues, studies, n_studies, responseName, het="PM", REG=True):
        self.HET = het

        self.responseName = responseName
        self.studies = studies
        self.singleStudyPvalues = Pvalues

        self.effects = np.arctanh(np.array(effects, dtype=np.float64)) #if REG else effects

        self.n_studies = n_studies
        self.n = float(len(studies))

        self.vi = np.array([(1./float(n-3)) for n in self.n_studies], dtype=np.float64)
        self.devs = np.sqrt( self.vi )
        self.CI_of_z = [ (z-( 1.96*dv ), z+( 1.96*dv )) for z,dv in zip(self.effects, self.devs)]
        
        self.w = np.array([(1./float(v)) for v in self.vi], dtype=np.float64)
        mu_bar = np.sum(a*b for a,b in zip(self.w, self.effects))/np.sum(self.w)
        self.Q = np.sum(a*b for a,b in zip(self.w, [(x - mu_bar)**2 for x in self.effects]))
        H = np.sqrt(self.Q/(self.n - 1))
        self.I2 = np.max([0., (self.Q-(len(self.vi)-1))/float(self.Q)])
        self.t2_PM, self.t2PM_conv = paule_mandel_tau(self.effects, self.vi)
        self.t2_DL = ((self.Q - self.n + 1) / self.scaling( self.w )) if (self.Q > (self.n-1)) else 0.

        if self.HET == "PM":
            self.W = [(1./float(v+self.t2_PM)) for v in self.vi]
        elif self.HET.startswith("FIX"):
            self.W = [(1./float(v)) for v in self.vi]
        else:
            self.W = [(1./float(v+self.t2_DL)) for v in self.vi]

        self.RE = self.CombinedEffect()
        self.stdErr = self.StdErrCombinedEffect(self.CombinedEffectVar())
        self.Zscore = self.CombinedEffectZScore(self.RE, self.stdErr)
        self.REvar = self.CombinedEffectVar()
        self.Pval = self.Pvalue(self.Zscore)
        self.conf_int = self.CombinedEffectConfInt(self.RE, self.stdErr)

        self.result = self.nice_shape(True)

    def tot_var(self, Effects, Weights):
        Q = np.sum(Weights * [x**2 for x in Effects]) - ((np.sum(Weights*Effects)**2)/np.sum(Weights))
        return Q

    def scaling(self, W):
        C = np.sum(W) - (np.sum([w**2 for w in W])/float(np.sum(W)))
        return C

    def tau_squared_DL(self, Q, df, C):
        return (Q-df)/float(C) if (Q>df) else 0.

    def CombinedEffect(self):
        return np.sum(self.W*self.effects)/float(np.sum(self.W))

    def CombinedEffectVar(self):
        return 1/float(np.sum(self.W))

    def StdErrCombinedEffect(self, CVar):
        return np.sqrt(CVar)

    def CombinedEffectConfInt(self, CE, SE):
        low = CE - 1.96*SE
        upp = CE + 1.96*SE
        return low, upp

    def CombinedEffectZScore(self, CE, SE):
        return CE / float(SE)

    def Pvalue(self, Z):
        return 2.*(1 - sts.norm.cdf(np.abs(Z)))

    def nice_shape(self, REG):
        NS = {}
        for rho,ci,P,study in zip(self.effects, self.CI_of_z, self.singleStudyPvalues, self.studies):
            eff = rho if (not REG) else np.tanh(rho)
            NS[study + "_Correlation"] = eff
            NS[study + "_Pvalue"] = P
            NS[study + "_conf_int"] = ";".join(list(map(str, [np.tanh(c) for c in ci])))

        NS["RE_Correlation"] = np.tanh(self.RE)
        NS["RE_Pvalue"] = self.Pval
        NS["RE_stdErr"] = np.tanh(self.stdErr)
        NS["RE_conf_int"] = ";".join(list(map(str, [np.tanh(c) for c in self.conf_int])))
        NS["RE_Var"] = self.REvar if (not REG) else np.tanh(self.REvar)
        NS["Zscore"] = self.Zscore if (not REG) else np.tanh(self.Zscore)
        NS["Tau2_DL"] = self.t2_DL
        NS["Tau2_PM"] = self.t2_PM
        NS["I2"] = self.I2
        NS = pd.DataFrame(NS, index=[self.responseName])
        return NS
