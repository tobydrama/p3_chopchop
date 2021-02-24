#!/usr/bin/env python3
import math

from Guide import Guide
from Bio.Seq import Seq
from Functions import gccontent

SCORE = {"INPAIR_OFFTARGET_0": 5000,
         "INPAIR_OFFTARGET_1": 3000,
         "INPAIR_OFFTARGET_2": 2000,
         "INPAIR_OFFTARGET_3": 1000,
         "OFFTARGET_PAIR_SAME_STRAND": 10000,
         "OFFTARGET_PAIR_DIFF_STRAND": 5000,
         "PAM_IN_PENALTY": 1000,
         "MAX_OFFTARGETS": 20000,  # FIX: SPECIFIC FOR TALEN AND CRISPR
         "COEFFICIENTS": 100,  # also used for RNA folding in ISOFORM mode
         "CRISPR_BAD_GC": 300,
         "FOLDING": 1}

STEM_LEN = 4

GC_LOW = 40
GC_HIGH = 70


class Cas9(Guide):

    def scoregRNA(seq, PAM, tail, lookup):
        """ Calculate score from model coefficients. score is 0-1, higher is better """
        score = 0
        if lookup in "Intercept":
            score = lookup["Intercept"]

        seq = seq[::-1]  # we calculate from PAM in a way: 321PAM123

        if lookup in "gc_low":
            gc = seq[:20].count('G') + seq[:20].count('C')
            if gc < 10:
                score = score + (abs(gc - 10) * lookup["gc_low"])
            elif gc > 10:
                score = score + ((gc - 10) * lookup["gc_high"])

        for i in range(len(seq)):
            key = seq[i] + str(i + 1)
            if lookup in key:
                score += lookup[key]

            if i + 1 < len(seq):
                double_key = seq[i] + seq[i + 1] + str(i + 1)
                if lookup in double_key:
                    score += lookup[double_key]

            if i == 0:
                double_key = PAM[0] + seq[0] + str(0)
                if lookup in double_key:
                    score += lookup[double_key]

        for i in range(len(PAM)):
            key = 'PAM' + PAM[i] + str(i + 1)
            if lookup in key:
                score += lookup[key]

            if i + 1 < len(PAM):
                double_key = 'PAM' + PAM[i] + PAM[i + 1] + str(i + 1)
                if lookup in double_key:
                    score += lookup[double_key]

        for i in range(len(tail)):
            key = str(i + 1) + tail[i]
            if lookup in key:
                score += lookup[key]

            if i + 1 < len(tail):
                double_key = str(i + 1) + tail[i] + tail[i + 1]
                if lookup in double_key:
                    score += lookup[double_key]

        score = 1 / (1 + math.e ** -score)
        return score

    def __init__(self, *args, **kwargs):
        super(Cas9, self).__init__(*args, **kwargs)
        self.CoefficientsScore = {"XU_2015": 0,
                                  "DOENCH_2014": 0,
                                  "DOENCH_2016": 0,
                                  "MORENO_MATEOS_2015": 0,
                                  "CHARI_2015": 0,
                                  "G_20": 0,
                                  "ALKAN_2018": 0,
                                  "ZHANG_2019": 0}
        self.repProfile = None  # Shen et al 2018 prediction of repair profile
        self.repStats = None

        if self.scoringMethod not in ["CHARI_2015", "DOENCH_2016", "ALKAN_2018", "ZHANG_2019", "ALL"]:
            self.CoefficientsScore[self.scoringMethod] = self.scoregRNA(
                self.downstream5prim + self.strandedGuideSeq[:-len(self.PAM)],
                self.strandedGuideSeq[-len(self.PAM):], self.downstream3prim, globals()[self.scoringMethod])
            self.score -= self.CoefficientsScore[self.scoringMethod] * SCORE['COEFFICIENTS']

        if self.scoringMethod == "ALKAN_2018" or self.scoringMethod == "ALL":
            from CRISPRoff.CRISPRoff_specificity import CRISPRoff_score
            self.CoefficientsScore[self.scoringMethod] = CRISPRoff_score(self.strandedGuideSeq)
            self.score -= self.CoefficientsScore[self.scoringMethod] * SCORE['COEFFICIENTS']
            #testing shit
            #CRISPR_CALC.CRISPR_CALCULATIONS_COEFFICIENTSSCORE(self)
            #CRISPR_CALC.CRISPR_CALCULATIONS_SCORE(self)


        if self.scoringMethod == "ALL":
            for met in ["XU_2015", "DOENCH_2014", "MORENO_MATEOS_2015", "G_20"]:
                self.CoefficientsScore[met] = self.scoregRNA(
                    self.downstream5prim + self.strandedGuideSeq[:-len(self.PAM)],
                    self.strandedGuideSeq[-len(self.PAM):], self.downstream3prim, globals()[met])

    def __str__(self):
        self.sort_offTargets()
        if self.scoringMethod == "ALL":
            return "%s\t%s:%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (
                self.strandedGuideSeq,
                self.chrom,
                self.start,
                self.strand,
                self.GCcontent,
                self.folding,
                self.offTargetsMM[0],
                self.offTargetsMM[1],
                self.offTargetsMM[2],
                ">=" + str(self.offTargetsMM[3]) if self.isKmaxed else self.offTargetsMM[3],
                self.CoefficientsScore["XU_2015"],
                self.CoefficientsScore["DOENCH_2014"],
                self.CoefficientsScore["DOENCH_2016"],
                self.CoefficientsScore["MORENO_MATEOS_2015"],
                self.CoefficientsScore["CHARI_2015"],
                self.CoefficientsScore["G_20"],
                self.CoefficientsScore["ALKAN_2018"],
                self.CoefficientsScore["ZHANG_2019"])
        else:
            return "%s\t%s:%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s\t%.2f" % (self.strandedGuideSeq,
                                                                      self.chrom,
                                                                      self.start,
                                                                      self.strand,
                                                                      self.GCcontent,
                                                                      self.folding,
                                                                      self.offTargetsMM[0],
                                                                      self.offTargetsMM[1],
                                                                      self.offTargetsMM[2],
                                                                      ">=" + str(
                                                                          self.offTargetsMM[3]) if self.isKmaxed else
                                                                      self.offTargetsMM[3],
                                                                      self.CoefficientsScore[self.scoringMethod])

    def calcSelfComplementarity(self, scoreSelfComp, backbone_regions, PAM, replace5prime=None):
        if replace5prime:
            fwd = replace5prime + self.strandedGuideSeq[len(replace5prime):(
                None if PAM == "" else -len(PAM))]  # Replace the 2 first bases with e.g. "GG"
        else:
            fwd = self.guideSeq[
                  0:(None if PAM == "" else -len(PAM))]  # Do not include PAM motif in folding calculations

        rvs = str(Seq(fwd).reverse_complement())
        L = len(fwd) - STEM_LEN - 1

        self.folding = 0

        for i in range(0, len(fwd) - STEM_LEN):
            if gccontent(fwd[i:i + STEM_LEN]) >= 0.5:
                if fwd[i:i + STEM_LEN] in rvs[0:(L - i)] or any(
                        [fwd[i:i + STEM_LEN] in item for item in backbone_regions]):
                    # sys.stderr.write("%s\t%s\n" % (fwd, fwd[i:i+STEM_LEN]))
                    self.folding += 1

        self.score += self.folding * SCORE['FOLDING']

    def calcGCContent(self, scoreGC):
        """ Calculate the GC content of the guide """
        if self.PAM is not None and self.strandedGuideSeq is not None:
            gSeq = self.strandedGuideSeq[0:(None if self.PAM == "" else -len(self.PAM))]
            Gcount = gSeq.count('G')
            Ccount = gSeq.count('C')
            self.GCcontent = (100 * (float(Gcount + Ccount) / int(len(gSeq))))
        else:
            self.GCcontent = 0

        if scoreGC:
            if self.GCcontent > GC_HIGH or self.GCcontent < GC_LOW:
                self.score += SCORE['CRISPR_BAD_GC']
