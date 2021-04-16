#!/usr/bin/env python3

from Vars import *
import re
from Bio.Seq import Seq
from operator import attrgetter

#So it doesn't crash
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
    for c in range(0, len(mismatch_string)-1):

        # If the character is a digit, check if the next character is a digit (>9) and add number to current
        if mismatch_string[c].isdigit():
            if mismatch_string[c+1].isdigit():
                current += (int(mismatch_string[c])*10)
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

    def __init__(self, name, flagSum, guideSize, guideSeq, scoreGC, scoreSelfComp,
                 backbone_regions, PAM, replace5prime=None, scoringMethod=None,
                 genome=None, gene=None, isoform=None, gene_isoforms=None, isKmaxed=False):

        self.isKmaxed = isKmaxed  # to print possibility of more mismatches
        self.scoringMethod = scoringMethod
        self.genome = genome
        self.gene = gene
        self.isoform = isoform
        self.gene_isoforms = gene_isoforms
        self.offTargetsIso = {0: set(), 1: set(), 2: set(), 3: set()}
        self.constitutive = False  # conservation of guide across all isoforms
        self.PAM = PAM
        # From the guide's name we can get the chromosome
        self.flagSum = str(flagSum)
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

        self.guideSize = guideSize
        self.targetSize = guideSize
        self.cluster = -1
        self.score = 0
        self.ALL_scores = [0, 0, 0, 0, 0, 0]
        self.meanBPP = 0  # in ISOFORM mode median of base pair probabilities

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
        guidePos = int(region[:location]) + 1

        # guide start coordinate
        self.start = self.exonStart + guidePos
        self.end = self.start + guideSize
        self.guideSeq = guideSeq

        # Record which strand the guide is on
        if self.flagSum == "16" or ISOFORMS:  # due to reverse complementing before alignments
            self.strandedGuideSeq = guideSeq
            if self.strand is None:
                self.strand = '+'
        else:
            self.strandedGuideSeq = str(Seq(guideSeq).reverse_complement())
            if self.strand is None:
                self.strand = '-'

        # Initiate offTargets list
        self.offTargets = []
        self.offTarget_hash = {}
        self.offTargets_sorted = False

        if scoreSelfComp:
            self.calc_self_complementarity(scoreSelfComp, backbone_regions, PAM, replace5prime)
        else:
            self.folding = "N/A"

        # Scoring
        self.calc_GC_content(scoreGC)

    def calc_self_complementarity(self, scoreSelfComp, backbone_regions, PAM, replace5prime = None):
        if replace5prime:
            fwd = self.strandedGuideSeq[len(PAM):-len(replace5prime)] + replace5prime #Replace the 2 first bases with e.g. "GG"
        else:
            fwd = self.guideSeq[len(PAM):] # Do not include PAM motif in folding calculations

        rvs = str(Seq(fwd).reverse_complement())
        L = len(fwd)-STEM_LEN-1

        self.folding = 0

        for i in range(0,len(fwd)-STEM_LEN):
            if gc_content(fwd[i:i + STEM_LEN]) >= 0.5:
                if fwd[i:i+STEM_LEN] in rvs[0:(L-i)] or any([fwd[i:i+STEM_LEN] in item for item in backbone_regions]):
                    #sys.stderr.write("%s\t%s\n" % (fwd, fwd[i:i+STEM_LEN]))
                    self.folding += 1

        self.score += self.folding * SCORE['FOLDING']


    def calc_GC_content(self, scoreGC):
        """ Calculate the GC content of the guide """
        if self.PAM is not None and self.strandedGuideSeq is not None:
            g_seq = self.strandedGuideSeq[len(self.PAM):]
            G_count = g_seq.count('G')
            C_count = g_seq.count('C')
            self.GCcontent = (100*(float(G_count+C_count)/int(len(g_seq))))
        else:
            self.GCcontent = 0

        if scoreGC:
            if self.GCcontent > GC_HIGH or self.GCcontent < GC_LOW:
                self.score += SCORE['CRISPR_BAD_GC']


    def add_off_target(self, hit, checkMismatch, maxOffTargets, countMMPos):
        """ Add off target hits (and not original hit) to list for each guide RNA """

        hit_id = "%s:%s" % (hit.chrom, hit.start)
        nmiss = 0
        mm_pattern = re.compile('NM:i:(\d+)')

        # If the hit is identical to the guide coord it is the original correct hit
        if self.chrom == hit.chrom and self.start == hit.start: # never true for isoforms
            # This is the original/main hit
            self.correct_hit = hit
            return

        if ISOFORMS and self.isoform == hit.chrom and self.strandedGuideSeq == hit.matchSeq:
            # This is the original/main hit
            self.correct_hit = hit
            return

        # Do not count off targets twice, e.g. for TALENs valid on both strands.
        if hit_id in self.offTarget_hash:
            return

        # Reverse count+allowed arrays if on the reverse strand
        if checkMismatch and hit.flagSum == 0 and not ISOFORMS:
            countMMPos = countMMPos[::-1]

        self.offTarget_hash[hit_id] = hit
        if checkMismatch:
            MMs = get_mismatch_pos(hit.mismatchPos[5:])
            for mm in MMs:
                if not countMMPos[mm]:
                    del(self.offTarget_hash[hit_id])
                    return

                elif not countMMPos[mm]:
                    nmiss += 1

        # Calculate score
        for opt in hit.opts:
            m = mm_pattern.match(opt)
            if m:
                mm = int(m.group(1)) - nmiss

                # ugly repeat to save time from iterating all isoforms
                if ISOFORMS and checkMismatch:
                    if self.gene_isoforms in hit.chrom: # and hit.chrom not in self.offTargetsIso[mm]:
                        self.offTargetsIso[mm].add(hit.chrom)
                        # don't count/score isoform mismatches but display which isoforms have them
                    else:
                        self.offTargetsMM[mm] += 1
                        self.score += SINGLE_OFFTARGET_SCORE[mm]
                else:
                    self.offTargetsMM[mm] += 1
                    self.score += SINGLE_OFFTARGET_SCORE[mm]

            if opt == "XM:i:" + str(maxOffTargets):
                self.score += SCORE['MAX_OFFTARGETS']
                self.offTargetsMM[0] += maxOffTargets
                self.offTargetsMM[1] += maxOffTargets
                self.offTargetsMM[2] += maxOffTargets
                self.offTargetsMM[3] += maxOffTargets

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
        if ISOFORMS:
            return "%s\t%s:%s\t%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.strandedGuideSeq,
                                                                                            self.chrom, self.start,
                                                                                            self.gene, self.isoform,
                                                                                            self.GCcontent, self.folding, self.meanBPP,
                                                                                            self.offTargetsMM[0], self.offTargetsMM[1],
                                                                                            self.offTargetsMM[2], self.offTargetsMM[3],
                                                                                            self.constitutive, (",").join(set(self.offTargetsIso[0])),
                                                                                            (",").join(set(self.offTargetsIso[1])),
                                                                                            (",").join(set(self.offTargetsIso[2])),
                                                                                            (",").join(set(self.offTargetsIso[3])))
        return "%s\t%s:%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s" % (self.strandedGuideSeq, self.chrom, self.start,
                                                                self.strand, self.GCcontent, self.folding,
                                                                self.offTargetsMM[0], self.offTargetsMM[1],
                                                                self.offTargetsMM[2],
                                                                ">=" + str(self.offTargetsMM[3]) if self.isKmaxed else self.offTargetsMM[3])

    def as_off_target_string(self, label, max_off_targets):
        self.sort_off_targets()
        off_targets = map(lambda x: x.as_off_target_string(label, max_off_targets), self.offTargets)

        return ";".join(off_targets)