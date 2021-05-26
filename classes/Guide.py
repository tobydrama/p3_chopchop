#!/usr/bin/env python3
import re
from operator import attrgetter

from Bio.Seq import Seq

import config
from constants import STEM_LEN, GC_HIGH, GC_LOW, SINGLE_OFFTARGET_SCORE


# So it doesn't crash
def gc_content(seq):
    gc = 0
    for i in seq:
        if i == 'G' or i == 'g' or i == 'C' or i == 'c':
            gc += 1
    return float(gc) / float(len(seq))


# Used in Guide
def get_mismatch_pos(mismatch_string):
    mismatches = []

    if mismatch_string.isdigit():
        return []

    current = 0
    for c in range(0, len(mismatch_string) - 1):

        # If the character is a digit, check if the next character is a digit (>9) and add number to current
        if mismatch_string[c].isdigit():
            if mismatch_string[c + 1].isdigit():
                current += (int(mismatch_string[c]) * 10)
            else:
                current += int(mismatch_string[c])

            # if character is a letter, it's a mismatch => add to results
        else:
            mismatches.append(current)
            current += 1

    # Include letters at the end
    if mismatch_string[-1].isalpha():
        mismatches.append(current)

    return mismatches


class Guide(object):
    """ This defines a class for each guide. The (off-target) hits for each guide form a separate class. The functions
    "add_off_target" and "sort_off_targets" applies to just the Tale class """

    def __init__(self, name, flag_sum, guide_size, guide_seq, score_gc, score_self_comp,
                 backbone_regions, pam, replace_5prime=None, scoring_method=None,
                 genome=None, gene=None, isoform=None, gene_isoforms=None, is_kmaxed=False):

        self.gc_content = 0
        self.is_kmaxed = is_kmaxed  # to print possibility of more mismatches
        self.scoring_method = scoring_method
        self.genome = genome
        self.gene = gene
        self.isoform = isoform
        self.gene_isoforms = gene_isoforms
        self.off_targets_iso = {0: set(), 1: set(), 2: set(), 3: set()}
        self.constitutive = False  # conservation of guide across all isoforms
        self.pam = pam
        # From the guide's name we can get the chromosome
        self.flag_sum = str(flag_sum)
        elements = name.split(":")
        self.id = elements[0]
        self.chrom = elements[1]
        coord = elements[2]

        self.name = ":".join(elements[0:3])

        if len(elements) > 3:
            self.downstream_5_prim = elements[3]
            self.downstream_3_prim = elements[4]
            self.strand = elements[5]
        else:
            self.downstream_5_prim = ''
            self.downstream_3_prim = ''
            self.strand = None

        self.guide_size = guide_size
        self.target_size = guide_size
        self.cluster = -1
        self.score = 0
        self.all_scores = [0, 0, 0, 0, 0, 0]
        self.mean_bpp = 0  # in ISOFORM mode median of base pair probabilities

        # Off target count
        self.off_targets_mm = [0] * 4

        # The location of the last digit of the exon start in the name string
        mid = coord.find('-')
        # The location of the first digit of the guide position in the exon
        end = (coord.find('_')) + 1

        # The full position of the guide in the exon
        region = coord[end:]

        # The location of the last digit of the guide position in the exon
        location = region.find('-')

        # The start of the exon containing the guide
        self.exon_start = int(coord[0:mid])

        # The number of bases after the exon start
        guide_pos = int(region[:location]) + 1

        # guide start coordinate
        self.start = self.exon_start + guide_pos
        self.end = self.start + guide_size
        self.guide_seq = guide_seq

        # Record which strand the guide is on
        if self.flag_sum == "16" or config.isoforms:  # due to reverse complementing before alignments
            self.stranded_guide_seq = guide_seq
            if self.strand is None:
                self.strand = '+'
        else:
            self.stranded_guide_seq = str(Seq(guide_seq).reverse_complement())
            if self.strand is None:
                self.strand = '-'

        # Initiate offTargets list
        self.off_targets = []
        self.off_target_hash = {}
        self.off_targets_sorted = False

        if score_self_comp:
            self.calc_self_complementarity(score_self_comp, backbone_regions, pam, replace_5prime)
        else:
            self.folding = "N/A"

        # Scoring
        self.calc_gc_content(score_gc)

    def calc_self_complementarity(self, score_self_comp, backbone_regions, pam, replace_5_prime=None):
        if replace_5_prime:
            # Replace the 2 first bases with e.g. "GG"
            fwd = self.stranded_guide_seq[len(pam):-len(replace_5_prime)] + replace_5_prime
        else:
            fwd = self.guide_seq[len(pam):]  # Do not include PAM motif in folding calculations

        rvs = str(Seq(fwd).reverse_complement())
        l = len(fwd) - STEM_LEN - 1

        self.folding = 0

        for i in range(0, len(fwd) - STEM_LEN):
            if gc_content(fwd[i:i + STEM_LEN]) >= 0.5:
                if fwd[i:i + STEM_LEN] in rvs[0:(l - i)] or any(
                        [fwd[i:i + STEM_LEN] in item for item in backbone_regions]):
                    self.folding += 1

        self.score += self.folding * config.score('FOLDING')

    def calc_gc_content(self, score_gc):
        """ Calculate the GC content of the guide """
        if self.pam is not None and self.stranded_guide_seq is not None:
            g_seq = self.stranded_guide_seq[len(self.pam):]
            g_count = g_seq.count('G')
            c_count = g_seq.count('C')
            self.gc_content = (100 * (float(g_count + c_count) / int(len(g_seq))))

        if score_gc:
            if self.gc_content > GC_HIGH or self.gc_content < GC_LOW:
                self.score += config.score('CRISPR_BAD_GC')

    def add_off_target(self, hit, check_mismatch, max_off_targets, count_mm_pos):
        """ Add off target hits (and not original hit) to list for each guide RNA """

        hit_id = "%s:%s" % (hit.chrom, hit.start)
        n_miss = 0
        mm_pattern = re.compile(r'NM:i:(\d+)')

        # If the hit is identical to the guide coord it is the original correct hit
        if self.chrom == hit.chrom and self.start == hit.start:  # never true for isoforms
            # This is the original/main hit
            self.correct_hit = hit
            return

        if config.isoforms and self.isoform == hit.chrom and self.stranded_guide_seq == hit.matchSeq:
            # This is the original/main hit
            self.correct_hit = hit
            return

        # Do not count off targets twice, e.g. for TALENs valid on both strands.
        if hit_id in self.off_target_hash:
            return

        # Reverse count+allowed arrays if on the reverse strand
        if check_mismatch and hit.flag_sum == 0 and not config.isoforms:
            count_mm_pos = count_mm_pos[::-1]

        self.off_target_hash[hit_id] = hit
        if check_mismatch:
            mms = get_mismatch_pos(hit.mismatch_pos[5:])
            for mm in mms:
                if not count_mm_pos[mm]:
                    del (self.off_target_hash[hit_id])
                    return

                elif not count_mm_pos[mm]:
                    n_miss += 1

        # Calculate score
        for opt in hit.opts:
            m = mm_pattern.match(opt)
            if m:
                mm = int(m.group(1)) - n_miss

                # ugly repeat to save time from iterating all isoforms
                if config.isoforms and check_mismatch:
                    if hit.chrom in self.gene_isoforms:  # and hit.chrom not in self.offTargetsIso[mm]:
                        self.off_targets_iso[mm].add(hit.chrom)
                        # don't count/score isoform mismatches but display which isoforms have them
                    else:
                        self.off_targets_mm[mm] += 1
                        self.score += SINGLE_OFFTARGET_SCORE[mm]
                else:
                    self.off_targets_mm[mm] += 1
                    self.score += SINGLE_OFFTARGET_SCORE[mm]

            if opt == "XM:i:" + str(max_off_targets):
                self.score += config.score('MAX_OFFTARGETS')
                self.off_targets_mm[0] += max_off_targets
                self.off_targets_mm[1] += max_off_targets
                self.off_targets_mm[2] += max_off_targets
                self.off_targets_mm[3] += max_off_targets

        self.off_targets_sorted = False

    def num_off_targets(self):
        """ Returns the number of off-target hits for each guide """
        self.sort_off_targets()
        return len(self.off_targets)

    def sort_off_targets(self):
        """ Sort off-target hits according to chromosome and genomic coordinate """

        if not self.off_targets_sorted:
            self.off_targets = self.off_target_hash.values()
            self.off_targets = sorted(self.off_targets, key=attrgetter('chrom', 'start'))
            self.off_targets_sorted = True

    def __str__(self) -> str:
        self.sort_off_targets()
        if config.isoforms:
            return "%s\t%s:%s\t%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
                self.stranded_guide_seq,
                self.chrom, self.start,
                self.gene, self.isoform,
                self.gc_content, self.folding, self.mean_bpp,
                self.off_targets_mm[0], self.off_targets_mm[1],
                self.off_targets_mm[2], self.off_targets_mm[3],
                self.constitutive, ",".join(set(self.off_targets_iso[0])),
                ",".join(set(self.off_targets_iso[1])),
                ",".join(set(self.off_targets_iso[2])),
                ",".join(set(self.off_targets_iso[3])))
        else:
            return "%s\t%s:%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s" % (
                self.stranded_guide_seq, self.chrom, self.start,
                self.strand, self.gc_content, self.folding,
                self.off_targets_mm[0], self.off_targets_mm[1],
                self.off_targets_mm[2],
                ">=" + str(self.off_targets_mm[3]) if self.is_kmaxed else self.off_targets_mm[3])

    def as_off_target_string(self, label, max_off_targets):
        self.sort_off_targets()
        off_targets = map(lambda x: x.as_off_target_string(label, max_off_targets), self.off_targets)

        return ";".join(off_targets)
