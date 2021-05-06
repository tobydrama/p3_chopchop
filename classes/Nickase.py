import config
from constants import SINGLE_OFFTARGET_SCORE
from functions.make_primers import find_restriction_sites


class Nickase:
    """ Pair class for 2 Cas9 that are the correct distance apart """

    def __init__(self, tale1, tale2, spacer_seq, spacer_size, off_target_pairs, enzyme_co, max_off_targets,
                 min_res_site_len):
        self.tale1 = tale1
        self.tale2 = tale2
        self.chrom = tale1.chrom
        self.strand = tale1.strand
        self.ID = ""
        self.restrictionSites = ""

        # Start of region covered by tale pair
        self.start = tale1.start

        # End of region covered by tale pair
        self.end = tale2.end
        self.spacerSeq = spacer_seq
        self.targetSize = spacer_size
        self.spacerSize = spacer_size
        self.offTargetPairs = off_target_pairs
        self.diffStrandOffTarget = 0
        self.sameStrandOffTarget = 0

        # Start cluster as -1, but will increment from 1
        self.cluster = -1
        self.spacerStart = tale1.start + tale1.guideSize
        self.spacerEnd = tale2.start - 1

        self.enzymeCo = enzyme_co
        self.strandedGuideSeq = str(self.tale1.guideSeq) + "\n" + self.spacerSeq + "\n" + str(self.tale2.guideSeq)
        self.offTargetPairCount = 0

        # Use bitwise operator to compare flag sum to see whether off-target TALEs are on different strands
        # (bad = good cutting ability) or on the same strand (not so bad = FokI domains probably too far apart to cut)
        indiv_score = 0

        for (hit1, hit2) in off_target_pairs:
            # Using boolean, count number of offtarget pairs on different strands
            if hit2.flagSum & hit1.flagSum == 0:
                self.diffStrandOffTarget += 1

            for opt in [hit1.opts, hit2.opts]:
                if opt == "NM:i:0":
                    indiv_score += SINGLE_OFFTARGET_SCORE[0]
                if opt == "NM:i:1":
                    indiv_score += SINGLE_OFFTARGET_SCORE[1]
                if opt == "NM:i:2":
                    indiv_score += SINGLE_OFFTARGET_SCORE[2]
                if opt == "NM:i:3":
                    indiv_score += SINGLE_OFFTARGET_SCORE[3]

        # Compute penalties (scores) for off-target hits. Worst = off-target pair, Not so bad = off-target single tale
        self.score = (self.diffStrandOffTarget * config.score(
            'OFFTARGET_PAIR_DIFF_STRAND')) + tale1.score + tale2.score - indiv_score + (
                                 tale1.strand == "+") * config.score('PAM_IN_PENALTY')
        res_sites = find_restriction_sites(self.spacerSeq, enzyme_co, min_res_site_len)
        self.restrictionSites = ";".join(map(lambda x: "%s:%s" % (str(x), ",".join(map(str, res_sites[x]))), res_sites))

    def __str__(self):
        # This creates a tab delimited list of output, with the final column as a semicolon-separated list of REs that
        # cut in the spacer
        sequence = str(self.tale1.guideSeq) + "*" + self.spacerSeq + "*" + str(self.tale2.guideSeq)

        return "%s\t%s:%s\t%s\t%s\t%s/%s\t%s/%s\t%s/%s\t%s/%s\t%s" % (
            sequence, self.chrom, self.start, self.cluster,
            len(self.offTargetPairs), self.tale1.offTargetsMM[0], self.tale2.offTargetsMM[0],
            self.tale1.offTargetsMM[1], self.tale2.offTargetsMM[1], self.tale1.offTargetsMM[2],
            self.tale2.offTargetsMM[2],
            ">=" + str(self.tale1.offTargetsMM[3]) if self.tale1.isKmaxed else self.tale1.offTargetsMM[3],
            ">=" + str(self.tale2.offTargetsMM[3]) if self.tale2.isKmaxed else self.tale2.offTargetsMM[3],
            self.restrictionSites)

    def as_off_target_string(self, label, max_off_targets):
        pairs = []

        # Add any off-target pairs
        if self.offTargetPairs:
            for offTargetPair in self.offTargetPairs:
                pairs.append("%s,%s" % (offTargetPair[0].as_off_target_string(label, max_off_targets),
                                        offTargetPair[1].as_off_target_string(label, max_off_targets)))
        else:
            pairs.append("")

        pairs = ";".join(pairs)

        return "\n".join([pairs, self.tale1.as_off_target_string("TALE 1", max_off_targets),
                          self.tale2.as_off_target_string("TALE 2", max_off_targets)])
