#!/usr/bin/env python3
import math

from Bio.Seq import Seq

import config
from classes.Guide import Guide
from constants import *
from functions.evaluate import gc_content


class Cas9(Guide):

    def __init__(self, *args, **kwargs):
        super(Cas9, self).__init__(*args, **kwargs)
        self.GCcontent = 0
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

    def __str__(self):
        self.sort_off_targets()
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

    def calc_self_complementarity(self, score_self_comp, backbone_regions, PAM, replace_5prime=None):
        if replace_5prime:
            fwd = replace_5prime + self.strandedGuideSeq[len(replace_5prime):(
                None if PAM == "" else -len(PAM))]  # Replace the 2 first bases with e.g. "GG"
        else:
            fwd = self.guideSeq[
                  0:(None if PAM == "" else -len(PAM))]  # Do not include PAM motif in folding calculations

        rvs = str(Seq(fwd).reverse_complement())
        L = len(fwd) - STEM_LEN - 1

        self.folding = 0

        for i in range(0, len(fwd) - STEM_LEN):
            if gc_content(fwd[i:i + STEM_LEN]) >= 0.5:
                if fwd[i:i + STEM_LEN] in rvs[0:(L - i)] or any(
                        [fwd[i:i + STEM_LEN] in item for item in backbone_regions]):
                    # sys.stderr.write("%s\t%s\n" % (fwd, fwd[i:i+STEM_LEN]))
                    self.folding += 1

        self.score += self.folding * config.score('FOLDING')

    def calc_GC_content(self, score_GC):
        """ Calculate the GC content of the guide """
        if self.PAM is not None and self.strandedGuideSeq is not None:
            g_seq = self.strandedGuideSeq[0:(None if self.PAM == "" else -len(self.PAM))]
            G_count = g_seq.count('G')
            C_count = g_seq.count('C')
            self.GCcontent = (100 * (float(G_count + C_count) / int(len(g_seq))))

        if score_GC:
            if self.GCcontent > GC_HIGH or self.GCcontent < GC_LOW:
                self.score += config.score('CRISPR_BAD_GC')
