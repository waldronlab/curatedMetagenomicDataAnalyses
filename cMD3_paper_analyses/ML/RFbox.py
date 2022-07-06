#!/usr/bin/env python

import subprocess as sb
import os

class RFreg(object):
    BASE_ARC = "../../../metaml/"
    #"/shares/CIBIO-Storage/CM/news/users/paolo.manghi/git/metaml/"
    metamlREG_exec = "regression_dev.py"

    def __init__(self, response, input_, output_, target=None, nt=10000, nsl=5, mf=0.01, cc="40,80,120", ncores=10, runs=3, featid="s__"):
        self.string = ["python " + os.path.join(self.BASE_ARC, self.metamlREG_exec)]
        self.string += [" %s %s " %(input_, output_)]
        self.string += [" -l rf -nsl %i -nt %i -mf %s -cc %s -nc %i -z %s" %(nsl, nt, mf, cc, ncores, featid)]
        self.string += [" -d %s -r 10" %response]
        if target:
            self.string += [" -t %s " %target]

    def disable_feats_ranking(self):
        self.string += [" -df "]

    def run(self, input_, outfile, target, outstream=None):
        #if not filename:
        self.string[1] = " %s %s " %(input_, outfile)
        self.string += [" -t %s " %(target)]
        if not outstream:
            sb.call("".join(self.string), shell=True)
        else:
            outstream.write("".join(self.string) + "\n")
 
    def runCV(self):
        sb.call("".join(self.string), shell=True)

    #def run_normalizing_target(self, outfile, target):
    #    self.string[1] = " %s %s " %(input_, outfile)
    #    self.string[4] = " -d %s " %(target)
    #    sb.call(self.string, shell=True)


class RFcls(object):
    BASE_ARC = "../../../metaml/"
    metamlCLS_exec = "classification_thomas-manghi.py"

    def __init__(self, response, input_, output_, target=None, nt=10000, nsl=5, mf=0.01, cc="40,80,120", ncores=10, runs=3, featid="s__"):
        self.string = ["python " + os.path.join(self.BASE_ARC, self.metamlCLS_exec)]
        self.string += [" %s %s " %(input_, output_)]
        self.string += [" -l rf -nsl %i -nt %i -mf %s -cc %s -nc %i -z %s" %(nsl, nt, mf, cc, ncores, featid)]
        self.string += [" -d %s -r %i --no_norm" %(response, runs)]
        if target:
            self.string += [" -t %s " %target]

    def disable_feats_ranking(self):
        self.string += [" -df "]

    def run(self, input_, outfile, target, outstream=None):
        #if not filename:
        self.string[1] = " %s %s " %(input_, outfile)
        self.string += [" -t %s " %(target)]
        print("".join(self.string), "Mi piacerebbe chiamare... ")
        if not outstream:
            sb.call("".join(self.string), shell=True)
        else:
            outstream.write("".join(self.string) + "\n")

    def runCV(self):
        sb.call("".join(self.string), shell=True)

if __name__ == "__main__":
    print("trovati un altro pollo da spennare")
