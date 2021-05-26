from Bio.Seq import Seq

import config
from functions.make_primers import find_restriction_sites


class Pair:
    """ Pair class for 2 TALEs that are the correct distance apart """

    def __init__(self, tale1, tale2, spacer_seq, spacer_size, off_target_pairs, enzyme_co, max_off_targets, g_rvd,
                 min_res_site_len):
        self.tale1 = tale1
        self.tale2 = tale2
        self.chrom = tale1.chrom
        self.strand = tale1.strand
        self.id = ""
        self.tale1.rvd = ""
        self.tale2.rvd = ""
        self.restriction_sites = ""

        # Start of region covered by tale pair
        self.start = tale1.start

        # End of region covered by tale pair
        self.end = tale2.end  # + tale2.guideSize
        self.spacer_seq = spacer_seq
        self.target_size = spacer_size
        self.spacer_size = spacer_size
        self.off_target_pairs = off_target_pairs
        self.diff_strand_off_target = 0
        self.same_strand_off_target = 0

        # Start cluster as -1, but will increment from 1
        self.cluster = -1
        self.spacer_start = tale1.start + tale1.guide_size
        self.spacer_end = tale2.start - 1

        self.enzyme_co = enzyme_co
        self.stranded_guide_seq = str(self.tale1.guide_seq) + "\n" + self.spacer_seq + "\n" + str(self.tale2.guide_seq)

        # Calculate RVD for TALEs; FIX: use mapping
        for base in tale1.guide_seq:
            if base == "A":
                tale1.rvd += "NI "
            elif base == "T":
                tale1.rvd += "NG "
            elif base == "C":
                tale1.rvd += "HD "
            elif base == "G":
                tale1.rvd += g_rvd

        for base in Seq(tale2.guide_seq).reverse_complement():
            if base == "A":
                tale2.rvd += "NI "
            elif base == "T":
                tale2.rvd += "NG "
            elif base == "C":
                tale2.rvd += "HD "
            elif base == "G":
                tale2.rvd += g_rvd

        self.offTargetPairCount = 0

        # Use bitwise operator to compare flag sum to see whether off-target TALEs are on different strands
        # (bad = good cutting ability) or on the same strand (not so bad = FokI domains probably too far apart to cut)
        indiv_score = 0

        for (hit1, hit2) in off_target_pairs:
            # Using boolean, count number of offtarget pairs on different strands
            if hit2.flag_sum & hit1.flag_sum == 0:
                self.diff_strand_off_target += 1

            # Using boolean, count number of offtarget pairs on same strand
            elif hit2.flag_sum & hit1.flag_sum == 16:
                self.same_strand_off_target += 1

            for opt in [hit1.opts, hit2.opts]:
                if opt == "NM:i:0":
                    indiv_score += config.score('INPAIR_OFFTARGET_0')
                if opt == "NM:i:1":
                    indiv_score += config.score('INPAIR_OFFTARGET_1')
                if opt == "NM:i:2":
                    indiv_score += config.score('INPAIR_OFFTARGET_2')
                if opt == "NM:i:3":
                    indiv_score += config.score('INPAIR_OFFTARGET_3')

        # Compute penalties (scores) for off-target hits. Worst = off-target pair, Not so bad = off-target single tale
        self.score = (self.same_strand_off_target * config.score('OFFTARGET_PAIR_SAME_STRAND')) + \
                     (self.diff_strand_off_target * config.score('OFFTARGET_PAIR_DIFF_STRAND')) + \
                     (tale1.score + tale2.score + indiv_score)
        res_sites = find_restriction_sites(self.spacer_seq, enzyme_co, min_res_site_len)
        self.restriction_sites = ";".join(
            map(lambda x: "%s:%s" % (str(x), ",".join(map(str, res_sites[x]))), res_sites))

    def __str__(self):
        # This creates a tab delimited list of output, with the final column as a semicolon-separated list of REs that
        # cut in the spacer
        sequence = str(self.tale1.guide_seq) + "*" + self.spacer_seq + "*" + str(self.tale2.guide_seq)

        return "%s\t%s:%s\t%s\t%s\t%s\t%s\t%s/%s\t%s/%s\t%s/%s\t%s/%s\t%s" % (
            sequence, self.chrom, self.start, self.tale1.rvd,
            self.tale2.rvd, self.cluster, len(self.off_target_pairs), self.tale1.off_targets_mm[0],
            self.tale2.off_targets_mm[0], self.tale1.off_targets_mm[1], self.tale2.off_targets_mm[1],
            self.tale1.off_targets_mm[2], self.tale2.off_targets_mm[2],
            ">=" + str(self.tale1.off_targets_mm[3]) if self.tale1.is_kmaxed else self.tale1.off_targets_mm[3],
            ">=" + str(self.tale2.off_targets_mm[3]) if self.tale2.is_kmaxed else self.tale2.off_targets_mm[3],
            self.restriction_sites)

    def as_off_target_string(self, label, max_off_targets):
        pairs = []

        # Add any off-target pairs
        if self.off_target_pairs:
            for off_target_pair in self.off_target_pairs:
                pairs.append("%s,%s" % (off_target_pair[0].as_off_target_string(label, max_off_targets),
                                        off_target_pair[1].as_off_target_string(label, max_off_targets)))
        else:
            pairs.append("")

        pairs = ";".join(pairs)

        return "\n".join([pairs, self.tale1.as_off_target_string("TALE 1", max_off_targets),
                          self.tale2.as_off_target_string("TALE 2", max_off_targets)])
