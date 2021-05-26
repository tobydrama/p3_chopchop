#!/usr/bin/env python3

from Bio.Seq import Seq

import config
from classes.Guide import Guide
from constants import *
from functions.evaluate import gc_content


class Cas9(Guide):

    def __init__(self, *args, **kwargs):
        super(Cas9, self).__init__(*args, **kwargs)
        self.coefficients_score = {"XU_2015": 0,
                                   "DOENCH_2014": 0,
                                   "DOENCH_2016": 0,
                                   "MORENO_MATEOS_2015": 0,
                                   "CHARI_2015": 0,
                                   "G_20": 0,
                                   "ALKAN_2018": 0,
                                   "ZHANG_2019": 0}
        self.repair_profile = None  # Shen et al 2018 prediction of repair profile
        self.repair_stats = None

    def __str__(self):
        self.sort_off_targets()
        if self.scoring_method == "ALL":
            return "%s\t%s:%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (
                self.stranded_guide_seq, self.chrom, self.start, self.strand, self.gc_content, self.folding,
                self.off_targets_mm[0], self.off_targets_mm[1], self.off_targets_mm[2],
                ">=" + str(self.off_targets_mm[3]) if self.is_kmaxed else self.off_targets_mm[3],
                self.coefficients_score["XU_2015"], self.coefficients_score["DOENCH_2014"],
                self.coefficients_score["DOENCH_2016"], self.coefficients_score["MORENO_MATEOS_2015"],
                self.coefficients_score["CHARI_2015"], self.coefficients_score["G_20"],
                self.coefficients_score["ALKAN_2018"], self.coefficients_score["ZHANG_2019"]
            )
        else:
            return "%s\t%s:%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s\t%.2f" % (
                self.stranded_guide_seq, self.chrom, self.start, self.strand, self.gc_content, self.folding,
                self.off_targets_mm[0], self.off_targets_mm[1], self.off_targets_mm[2],
                ">=" + str(self.off_targets_mm[3]) if self.is_kmaxed else self.off_targets_mm[3],
                self.coefficients_score[self.scoring_method]
            )

    def calc_self_complementarity(self, score_self_comp, backbone_regions, pam, replace_5prime=None):
        if replace_5prime:
            fwd = replace_5prime + self.stranded_guide_seq[len(replace_5prime):(
                None if pam == "" else -len(pam))]  # Replace the 2 first bases with e.g. "GG"
        else:
            fwd = self.guide_seq[
                  0:(None if pam == "" else -len(pam))]  # Do not include PAM motif in folding calculations

        rvs = str(Seq(fwd).reverse_complement())
        l = len(fwd) - STEM_LEN - 1

        self.folding = 0

        for i in range(0, len(fwd) - STEM_LEN):
            if gc_content(fwd[i:i + STEM_LEN]) >= 0.5:
                if fwd[i:i + STEM_LEN] in rvs[0:(l - i)] or any(
                        [fwd[i:i + STEM_LEN] in item for item in backbone_regions]):
                    # sys.stderr.write("%s\t%s\n" % (fwd, fwd[i:i+STEM_LEN]))
                    self.folding += 1

        self.score += self.folding * config.score('FOLDING')

    def calc_gc_content(self, score_gc):
        """ Calculate the GC content of the guide """
        if self.pam is not None and self.stranded_guide_seq is not None:
            g_seq = self.stranded_guide_seq[0:(None if self.pam == "" else -len(self.pam))]
            g_count = g_seq.count('G')
            c_count = g_seq.count('C')
            self.gc_content = (100 * (float(g_count + c_count) / int(len(g_seq))))

        if score_gc:
            if self.gc_content > GC_HIGH or self.gc_content < GC_LOW:
                self.score += config.score('CRISPR_BAD_GC')
