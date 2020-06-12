import math
import sys,getopt
import os
import numpy as np
from scipy.stats import rankdata


inputfile = sys.argv[1]
outputfile = sys.argv[2]


min_rep_coverage = 20

class X2(object):
    def __init__(self):
        pass

    def approx_gamma(self,Z):
        RECIP_E = 0.36787944117144232159552377016147
        TWOPI = 6.283185307179586476925286766559
        D = 1.0 / (10.0 * Z)
        D = 1.0 / ((12 * Z) - D)
        D = (D + Z) * RECIP_E;
        D = math.pow(D, Z);
        D *= math.sqrt(TWOPI / Z)
        return D

    def igf(self,S,Z):
        if Z < 0.0:
            return 0.0
        Sc = (1.0 / S)
        Sc *= math.pow(Z, S)
        Sc *= math.exp(-Z)
        Sum = Nom = Denom = 1.0
        for I in range(200):
            Nom *= Z
            S += 1
            Denom *= S
            Sum += (Nom / Denom)
        return Sum * Sc

    def chisqr(self,Dof,Cv):
        if Cv < 0 or Dof < 1:
            return 0.0
        K = Dof * 0.5
        X = Cv * 0.5
        if Dof == 2:
            return math.exp(-1.0 * X)
        PValue = self.igf(K,X)
        if np.isnan(PValue) or math.isinf(PValue) or PValue <= 1e-8:
            return 1e-14
        PValue /= self.approx_gamma(K)
        return 1.0 - PValue

    def get_pvalue(self,FirstOriginal,SecondOriginal):
        dispots = len(FirstOriginal)
        SumFirst, SumSecond, SumTotal = 0,0,0
        SumBoth = [0] * dispots
        ExpFirst = [0] * dispots
        ExpSecond  = [0] * dispots
        for i in range(dispots):
            SumBoth[i] = FirstOriginal[i] + SecondOriginal[i]
            SumFirst += FirstOriginal[i]
            SumSecond += SecondOriginal[i]
        SumTotal = SumFirst + SumSecond

        for i in range(dispots):
            ExpFirst[i] = SumBoth[i] * SumFirst / SumTotal
            ExpSecond[i] = SumBoth[i] * SumSecond / SumTotal

        result = 0
        Dregree = 0
        for i in range(dispots):
            if FirstOriginal[i] + SecondOriginal[i] > 0.0:
                Dregree += 1
                if ExpFirst[i]:
                    result += (FirstOriginal[i] - ExpFirst[i]) * (FirstOriginal[i] - ExpFirst[i]) / ExpFirst[i]
                if ExpSecond[i]:
                    result += (SecondOriginal[i] - ExpSecond[i]) * (SecondOriginal[i] - ExpSecond[i]) / ExpSecond[i]

        if Dregree == 1 or result ==0:
            PValue = 1.0
        else:
            Dregree -= 1
            PValue = self.chisqr(Dregree,result)

        if PValue < 0:
            PValue = PValue * (-1)

        return PValue


    def load_dis(self):

        f_somatic_dis = open(outputfile,"w")
        pvalue_list = []
        site_list = []
        in_file_feature = open(inputfile, "r")
        l = in_file_feature.readline()
        while l:
            tmp_site = l
            data = l.strip().split(" ")
            site = data[0] + "_" +data[1]

            l = in_file_feature.readline()
            if l[0] != "N":
                print("distribution file error...")
                exit()
            tmp_N = l
            FirstOriginal = list(map(float, l.strip().split(" ")[1::]))
            if sum(FirstOriginal) < min_rep_coverage:
                l = in_file_feature.readline()
                l = in_file_feature.readline()
                continue

            l = in_file_feature.readline()
            if l[0] != "T":
                print("distribution file error...")
                exit()
            tmp_T = l
            SecondOriginal = list(map(float, l.strip().split(" ")[1::]))

            if sum(SecondOriginal) < min_rep_coverage:
                l = in_file_feature.readline()
                continue

            l = in_file_feature.readline()
            PValue = self.get_pvalue(FirstOriginal,SecondOriginal)

            if PValue < 0.05:
                f_somatic_dis.write(tmp_site)
                f_somatic_dis.write(tmp_N)
                f_somatic_dis.write(tmp_T)
        f_somatic_dis.close()


if __name__ == '__main__':

    x = X2()
    x.load_dis()
