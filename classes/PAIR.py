from Bio.Seq import Seq
from functions.make_primers import find_restriction_sites
import config

class Pair:
    """ Pair class for 2 TALEs that are the correct distance apart """
    def __init__(self, tale1, tale2, spacerSeq, spacerSize, offTargetPairs, enzymeCo, maxOffTargets, g_RVD, minResSiteLen):
        self.tale1 = tale1
        self.tale2 = tale2
        self.chrom = tale1.chrom
        self.strand = tale1.strand
        self.ID = ""
        self.tale1.rvd = ""
        self.tale2.rvd = ""
        self.restrictionSites = ""

        # Start of region covered by tale pair
        self.start = tale1.start

        # End of region covered by tale pair
        self.end = tale2.end  # + tale2.guideSize
        self.spacerSeq = spacerSeq
        self.targetSize = spacerSize
        self.spacerSize = spacerSize
        self.offTargetPairs = offTargetPairs
        self.diffStrandOffTarget = 0
        self.sameStrandOffTarget = 0

        # Start cluster as -1, but will increment from 1
        self.cluster = -1
        self.spacerStart = tale1.start + tale1.guideSize
        self.spacerEnd = tale2.start - 1

        self.enzymeCo = enzymeCo
        self.strandedGuideSeq = str(self.tale1.guideSeq) + "\n" + self.spacerSeq + "\n" + str(self.tale2.guideSeq)

        # Calculate RVD for TALEs; FIX: use mapping
        for base in tale1.guideSeq:
            if base == "A":
                tale1.rvd += "NI "
            elif base == "T":
                tale1.rvd += "NG "
            elif base == "C":
                tale1.rvd += "HD "
            elif base == "G":
                tale1.rvd += g_RVD

        for base in Seq(tale2.guideSeq).reverse_complement():
            if base == "A":
                tale2.rvd += "NI "
            elif base == "T":
                tale2.rvd += "NG "
            elif base == "C":
                tale2.rvd += "HD "
            elif base == "G":
                tale2.rvd += g_RVD

        self.offTargetPairCount = 0

        # Use bitwise operator to compare flag sum to see whether off-target TALEs are on different strands (bad = good cutting ability)
        # or on the same strand (not so bad = FokI domains probably too far apart to cut)
        indivScore = 0

        for (hit1, hit2) in offTargetPairs:
            # Using boolean, count number of offtarget pairs on different strands
            if hit2.flagSum & hit1.flagSum == 0:
                self.diffStrandOffTarget += 1

            # Using boolean, count number of offtarget pairs on same strand
            elif hit2.flagSum & hit1.flagSum == 16:
                self.sameStrandOffTarget += 1

            for opt in [hit1.opts, hit2.opts]:
                if opt == "NM:i:0":
                    indivScore += config.score('INPAIR_OFFTARGET_0')
                if opt == "NM:i:1":
                    indivScore += config.score('INPAIR_OFFTARGET_1')
                if opt == "NM:i:2":
                    indivScore += config.score('INPAIR_OFFTARGET_2')
                if opt == "NM:i:3":
                    indivScore += SCORE['INPAIR_OFFTARGET_3']

        # Compute penalties (scores) for off-target hits. Worst = off-target pair, Not so bad = off-target single tale
        self.score = (self.sameStrandOffTarget * SCORE['OFFTARGET_PAIR_SAME_STRAND']) + (self.diffStrandOffTarget * SCORE['OFFTARGET_PAIR_DIFF_STRAND']) + tale1.score + tale2.score + indivScore
        resSites = find_restriction_sites(self.spacerSeq, enzymeCo, minResSiteLen)
        self.restrictionSites = ";".join(map(lambda x: "%s:%s" % (str(x), ",".join(map(str, resSites[x]))), resSites))

    def __str__(self):
        # This creates a tab delimited list of output, with the final column as a semicolon-separated list of REs that cut in the spacer
        sequence = str(self.tale1.guideSeq) + "*" + self.spacerSeq + "*" + str(self.tale2.guideSeq)

        return "%s\t%s:%s\t%s\t%s\t%s\t%s\t%s/%s\t%s/%s\t%s/%s\t%s/%s\t%s" % (
                sequence, self.chrom, self.start, self.tale1.rvd,
                self.tale2.rvd, self.cluster, len(self.offTargetPairs), self.tale1.offTargetsMM[0],
                self.tale2.offTargetsMM[0], self.tale1.offTargetsMM[1], self.tale2.offTargetsMM[1],
                self.tale1.offTargetsMM[2], self.tale2.offTargetsMM[2],
                ">=" + str(self.tale1.offTargetsMM[3]) if self.tale1.isKmaxed else self.tale1.offTargetsMM[3],
                ">=" + str(self.tale2.offTargetsMM[3]) if self.tale2.isKmaxed else self.tale2.offTargetsMM[3],
                self.restrictionSites)

    def as_off_target_string(self, label, maxOffTargets):
        pairs = []

        # Add any off-target pairs
        if self.offTargetPairs:
            for offTargetPair in self.offTargetPairs:
                pairs.append("%s,%s" % (offTargetPair[0].as_off_target_string(label, maxOffTargets), offTargetPair[1].as_off_target_string(label, maxOffTargets)))
        else:
            pairs.append("")

        pairs = ";".join(pairs)

        return "\n".join([pairs, self.tale1.as_off_target_string("TALE 1", maxOffTargets), self.tale2.as_off_target_string("TALE 2", maxOffTargets)])
