#!/usr/bin/env python3
import re
import uuid
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

    def __init__(self, name, flag_sum, guide_size, guide_seq, score_GC, score_self_comp,
                 backbone_regions, PAM, replace_5prime=None, scoring_method=None,
                 genome=None, gene=None, isoform=None, gene_isoforms=None, is_kmaxed=False, id=None):

        self.GCcontent = 0
        self.isKmaxed = is_kmaxed  # to print possibility of more mismatches
        self.scoringMethod = scoring_method
        self.genome = genome
        self.gene = gene
        self.isoform = isoform
        self.gene_isoforms = gene_isoforms
        self.offTargetsIso = {0: set(), 1: set(), 2: set(), 3: set()}
        self.constitutive = False  # conservation of guide across all isoforms
        self.PAM = PAM
        # From the guide's name we can get the chromosome
        self.flagSum = str(flag_sum)
        elements = name.split(":")
        self.ID = elements[0]
        self.chrom = elements[1]
        coord = elements[2]

        self.name = ":".join(elements[0:3])

        if len(elements) > 3:
            self.downstream5prim = elements[3]
            self.downstream3prim = elements[4]
            self.strand = elements[5]
        else:
            self.downstream5prim = ''
            self.downstream3prim = ''
            self.strand = None

        self.guideSize = guide_size
        self.targetSize = guide_size
        self.cluster = -1
        self.score = 0
        self.ALL_scores = [0, 0, 0, 0, 0, 0]
        self.meanBPP = 0  # in ISOFORM mode median of base pair probabilities

        # Unique identifier
        if id is None:
            self.uuid = str(uuid.uuid4())
        else:
            self.uuid = str(id)

        # Off target count
        self.offTargetsMM = [0] * 4

        # The location of the last digit of the exon start in the name string
        mid = coord.find('-')
        # The location of the first digit of the guide position in the exon
        end = (coord.find('_')) + 1

        # The full position of the guide in the exon
        region = coord[end:]

        # The location of the last digit of the guide position in the exon
        location = region.find('-')

        # The start of the exon containing the guide
        self.exonStart = int(coord[0:mid])

        # The number of bases after the exon start
        guide_pos = int(region[:location]) + 1

        # guide start coordinate
        self.start = self.exonStart + guide_pos
        self.end = self.start + guide_size
        self.guideSeq = guide_seq

        # Record which strand the guide is on
        if self.flagSum == "16" or config.isoforms:  # due to reverse complementing before alignments
            self.strandedGuideSeq = guide_seq
            if self.strand is None:
                self.strand = '+'
        else:
            self.strandedGuideSeq = str(Seq(guide_seq).reverse_complement())
            if self.strand is None:
                self.strand = '-'

        # Initiate offTargets list
        self.offTargets = []
        self.offTarget_hash = {}
        self.offTargets_sorted = False

        if score_self_comp:
            self.calc_self_complementarity(score_self_comp, backbone_regions, PAM, replace_5prime)
        else:
            self.folding = "N/A"

        # Scoring
        self.calc_GC_content(score_GC)

    def reinitialize_flag_sum(self, flag_sum):
        self.flagSum = flag_sum 
        
        if self.flagSum == "16" or config.isoforms:  # due to reverse complementing before alignments
            self.strandedGuideSeq = self.guideSeq
            if self.strand is None:
                self.strand = '+'
        else:
            self.strandedGuideSeq = str(Seq(self.guideSeq).reverse_complement())
            if self.strand is None:
                self.strand = '-'

    def calc_self_complementarity(self, score_self_comp, backbone_regions, PAM, replace5prime=None):
        if replace5prime:
            # Replace the 2 first bases with e.g. "GG"
            fwd = self.strandedGuideSeq[len(PAM):-len(replace5prime)] + replace5prime
        else:
            fwd = self.guideSeq[len(PAM):]  # Do not include PAM motif in folding calculations

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
            g_seq = self.strandedGuideSeq[len(self.PAM):]
            G_count = g_seq.count('G')
            C_count = g_seq.count('C')
            self.GCcontent = (100 * (float(G_count + C_count) / int(len(g_seq))))

        if score_GC:
            if self.GCcontent > GC_HIGH or self.GCcontent < GC_LOW:
                self.score += config.score('CRISPR_BAD_GC')

    def add_off_target(self, hit, check_mismatch, max_off_targets, count_MM_pos):
        """ Add off target hits (and not original hit) to list for each guide RNA """

        hit_id = "%s:%s" % (hit.chrom, hit.start)
        n_miss = 0
        mm_pattern = re.compile(r'NM:i:(\d+)')

        # If the hit is identical to the guide coord it is the original correct hit
        if self.chrom == hit.chrom and self.start == hit.start:  # never true for isoforms
            # This is the original/main hit
            self.correct_hit = hit
            return

        if config.isoforms and self.isoform == hit.chrom and self.strandedGuideSeq == hit.matchSeq:
            # This is the original/main hit
            self.correct_hit = hit
            return

        # Do not count off targets twice, e.g. for TALENs valid on both strands.
        if hit_id in self.offTarget_hash:
            return

        # Reverse count+allowed arrays if on the reverse strand
        if check_mismatch and hit.flagSum == 0 and not config.isoforms:
            count_MM_pos = count_MM_pos[::-1]

        self.offTarget_hash[hit_id] = hit
        if check_mismatch:
            MMs = get_mismatch_pos(hit.mismatchPos[5:])
            for mm in MMs:
                if not count_MM_pos[mm]:
                    del (self.offTarget_hash[hit_id])
                    return

                elif not count_MM_pos[mm]:
                    n_miss += 1

        # Calculate score
        for opt in hit.opts:
            m = mm_pattern.match(opt)
            if m:
                mm = int(m.group(1)) - n_miss

                # ugly repeat to save time from iterating all isoforms
                if config.isoforms and check_mismatch:
                    if hit.chrom in self.gene_isoforms:  # and hit.chrom not in self.offTargetsIso[mm]:
                        self.offTargetsIso[mm].add(hit.chrom)
                        # don't count/score isoform mismatches but display which isoforms have them
                    else:
                        self.offTargetsMM[mm] += 1
                        self.score += SINGLE_OFFTARGET_SCORE[mm]
                else:
                    self.offTargetsMM[mm] += 1
                    self.score += SINGLE_OFFTARGET_SCORE[mm]

            if opt == "XM:i:" + str(max_off_targets):
                self.score += config.score('MAX_OFFTARGETS')
                self.offTargetsMM[0] += max_off_targets
                self.offTargetsMM[1] += max_off_targets
                self.offTargetsMM[2] += max_off_targets
                self.offTargetsMM[3] += max_off_targets

        self.offTargets_sorted = False

    def num_off_targets(self):
        """ Returns the number of off-target hits for each guide """
        self.sort_off_targets()
        return len(self.offTargets)

    def sort_off_targets(self):
        """ Sort off-target hits according to chromosome and genomic coordinate """

        if self.offTargets_sorted:
            return

        self.offTargets = self.offTarget_hash.values()
        self.offTargets = sorted(self.offTargets, key=attrgetter('chrom', 'start'))
        self.offTargets_sorted = True

    def __str__(self):
        self.sort_off_targets()
        if config.isoforms:
            return "%s\t%s:%s\t%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
                self.strandedGuideSeq,
                self.chrom, self.start,
                self.gene, self.isoform,
                self.GCcontent, self.folding, self.meanBPP,
                self.offTargetsMM[0], self.offTargetsMM[1],
                self.offTargetsMM[2], self.offTargetsMM[3],
                self.constitutive, ",".join(set(self.offTargetsIso[0])),
                ",".join(set(self.offTargetsIso[1])),
                ",".join(set(self.offTargetsIso[2])),
                ",".join(set(self.offTargetsIso[3])))
        return "%s\t%s:%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s" % (
            self.strandedGuideSeq, self.chrom, self.start,
            self.strand, self.GCcontent, self.folding,
            self.offTargetsMM[0], self.offTargetsMM[1],
            self.offTargetsMM[2],
            ">=" + str(self.offTargetsMM[3]) if self.isKmaxed else self.offTargetsMM[3])

    def as_off_target_string(self, label, max_off_targets):
        self.sort_off_targets()
        off_targets = map(lambda x: x.as_off_target_string(label, max_off_targets), self.offTargets)

        return ";".join(off_targets)
