import re
import sys
import json
import os
import logging
import pandas

from Vars import PRIMER3_CONFIG, CONFIG, EXIT, PRIMER_OFF_TARGET_MIN
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from subprocess import Popen, PIPE
from classes.Guide import Guide
from classes.Hit import Hit
from operator import itemgetter
from Bio.Restriction import Analysis, RestrictionBatch


# Used in makePrimersFasta and makePrimersGenome
def parse_primer3_output(target, region, primer3output, primerFastaFile):
    posPattern = re.compile('PRIMER_(\w+)_(\d+)')
    attPattern = re.compile('PRIMER_(\w+)_(\d+)_(\w+)')
    primers = {}
    primerPos = {}

    for line in primer3output.split("\n"):
        if line[0] == "=":
            break

        label, value = line.split("=")
        m = attPattern.match(label)
        if m:
            primers[(m.group(2), m.group(1), m.group(3))] = value
        else:
            m = posPattern.match(label)

            if m:
                position, length = value.split(",")

                s, e = int(position), int(position) + int(length)
                if m.group(1) == "RIGHT":
                    s, e = int(position) - int(length) + 1, int(position) + 1
                primerPos[label] = [s, e]

                primerFastaFile.write(">%s_%s_%s:%s_%s-%s\n%s\n" % (
                    target.ID, m.group(2), m.group(1), region, s, e, primers[(m.group(2), m.group(1), "SEQUENCE")]))

    return primers, primerPos


# Used in makePrimersFasta and makePrimersGenome
def get_primer_options(options):
    # Parse primer3 options. Update config if known option, otherwise append to primer3 input file
    primerOpt = ""

    if options:
        for opt in options.split(","):
            key, value = opt.split("=")
            if key in PRIMER3_CONFIG:
                PRIMER3_CONFIG[key] = value
            else:
                primerOpt += opt + "\n"

    return primerOpt


# Used in makePrimersFasta
def get_primer_query_sequence_fasta(target, outputDir, flank, fastaSequence):
    s = target.start - flank
    e = target.end + flank
    seqLenBeforeTarget = flank

    if s < 0:
        seqLenBeforeTarget -= abs(s)
        s = 0

    if e > len(fastaSequence):
        e = len(fastaSequence)

    return fastaSequence[s:e], seqLenBeforeTarget


# Used in makePrimerGnome
def get_primer_query_sequence_2bit(target, outputDir, flank, genome, twoBitToFaIndexDir, strand):
    s = target.start - flank
    seqLenBeforeTarget = flank

    if s < 0:
        seqLenBeforeTarget -= abs(s)
        s = 0

    prog = Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2>> %s/twoBitToFa.err" % (
        CONFIG["PATH"]["TWOBITTOFA"], target.chrom, s, target.end + flank, twoBitToFaIndexDir, genome, outputDir),
                 stdout=PIPE, shell=True)
    output = prog.communicate()

    if prog.returncode != 0:
        sys.stderr.write("Running twoBitToFa failed\n")
        sys.exit(EXIT['TWOBITTOFA_ERROR'])

    output = output[0].decode().split("\n")
    del (output[0])
    seq = "".join(output)
    return seq, seqLenBeforeTarget


# Used in makePrimersFasta and makePrimersGenome
def run_bowtie_primers(primerFastaFileName, outputDir, genome, bowtieIndexDir, maxOffTargets):
    command = "%s -v 0 --best --sam-nohead -k 10 %s/%s -f %s -S %s/primer_results.sam 2> %s/bowtie_primers.err" % (
        CONFIG["PATH"]["BOWTIE"], bowtieIndexDir, genome, primerFastaFileName, outputDir, outputDir)
    prog = Popen(command, shell=True)
    prog.wait()

    if prog.returncode != 0:
        sys.stderr.write("Running bowtie on primers failed\n")
        sys.exit(EXIT['BOWTIE_PRIMER_ERROR'])

    return parse_bowtie(Guide, "%s/primer_results.sam" % outputDir, False, False, False, None, None,
                        maxOffTargets, None, None, False, None, None)


# Used in Nickase, Pair, dump_restriction_sites
def find_restriction_sites(sequence, enzymeCompany, minSize=1):
    # Take spacer_seq as DNA input for restriction site search
    mySeq = Seq(sequence)

    # Restricts enzyme possibilities to NEB enzymes. Can ultimately change to any supplier.
    rb = RestrictionBatch(first=[], suppliers=[enzymeCompany])

    # Filter binding sites shorter than given length
    rb = filter(lambda x: len(x) > minSize, rb)

    # Determine which restriction enzymes cut in the sequence provided
    analyze = Analysis(rb, mySeq)
    return analyze.with_sites()


# Used in makePrimersFasta and makePrimersGenome
def dump_restriction_sites(target, seq, flanks, enzymeCo, outputDir, minResSiteLen):
    sites = find_restriction_sites(seq, enzymeCo, minResSiteLen)
    out = [map(lambda x: [str(enzyme), x + target.start - flanks, enzyme.size], sites[enzyme]) for enzyme in sites]
    out = [item for sublist in out for item in sublist]
    out = sorted(out, key=itemgetter(1))

    # Assign tier to avoid overlaps
    siteCount = {}
    tiers = [0] * 23
    for site in out:
        tier = 0

        # count number of sites for each enzyme
        if not site[0] in siteCount:
            siteCount[site[0]] = 0
        siteCount[site[0]] += 1

        for j in range(len(tiers)):
            if site[1] > tiers[j]:
                tier = j
                tiers[j] = site[1] + site[2]
                break
        site.append(tier)

    # Assign colors depending on uniqueness
    for site in out:
        if siteCount[site[0]] == 1:
            site.append("green")
        else:
            site.append("red")

    outputFile = open("%s/restriction_%s.json" % (outputDir, target.ID), 'w')
    json.dump(out, outputFile)
    outputFile.close()

    return sites


# Used in make_primers_fasta and make_primers_Genome
def dump_locus_sequence(target, outputDir, seq, seqLenBeforeTarget, strand):
    if strand == "-":
        seq = str(Seq(seq).complement())
    out = [[target.start - seqLenBeforeTarget, target.end, seq]]
    outputFile = open("%s/locusSeq_%s.json" % (outputDir, target.ID), 'w')
    json.dump(out, outputFile)
    outputFile.close()


# Used in makePrimersFasta and makePrimersGenome
def dump_genbank_file(seq, target, restSites, primers, outputDir, geneID, lociStart, strand):
    name = "%s, locus %s" % (geneID, target.ID)
    desc = "CHOPCHOP prediction for gene %s, target %s" % (geneID, target.ID)
    annotation = {"organism": "Danio rerio", "Target location": "chrT:1-20"}

    # Genbank file
    genbankFile = open('%s/%s_%s.gb' % (outputDir, geneID, target.ID), 'w')
    record = SeqRecord(Seq(seq), description=desc, name="CHOPCHOP", id=name)
    record.annotation = annotation

    if target.strand == "+":
        ts = 1
    else:
        ts = -1

    record.features.append(
        SeqFeature(FeatureLocation(target.start - lociStart - 1, target.end - lociStart - 1, strand=ts),
                   type="Target"))

    for primer in primers:
        record.features.append(SeqFeature(FeatureLocation(primers[primer][0], primers[primer][1]), type=primer))

    if strand == "-":
        record = record.reverse_complement()

    # TODO: Have an actual bioinformatics person look at this.
    record.annotations["molecule_type"] = "DNA"

    SeqIO.write(record, genbankFile, "genbank")
    genbankFile.close()

    pass


def has_Off_targets(tale1, tale2, offTargetMin, offTargetMax):
    """ Returns the number of off-targets for a pair of TALENs (10-24bp apart) """

    offTargetPairs = []

    # Calls sort function to sort off-targets by chromosome and chromosome position.
    # Bowtie ranks them according to quality of hit
    tale1.sort_off_targets()
    tale2.sort_off_targets()

    # TODO: Eivind to write this code properly. Include a way to step backwards, so as not to miss any hits.
    # Need to make a queue..?
    for i in range(len(tale1.offTargets)):
        hit1 = tale1.offTargets[i]

        for j in range(len(tale2.offTargets)):
            hit2 = tale2.offTargets[j]

            # Determines whether 2 tales are on the same chromosome and 10-24 bp apart.
            if hit2.chrom == hit1.chrom and offTargetMin <= abs(hit2.start-hit1.start) <= offTargetMax:
                offTargetPairs.append([hit1, hit2])

    return offTargetPairs


# Used in makePrimersFasta and makePrimersGenome
def pair_primers(primerAttributes, primerList, outputDir):
    primers = {}

    for primer in primerList:
        guide, primerPairID, side = primer.ID.split("_")

        s = 0
        if side == "RIGHT":
            s = 1

        if guide not in primers:
            primers[guide] = {}

        if primerPairID not in primers[guide]:
            primers[guide][primerPairID] = [None, None]

        primers[guide][primerPairID][s] = primer

    for guideID in primers:
        guide = primers[guideID]

        att = primerAttributes[int(guideID)]

        outputFile = open("%s/primer_%s.json" % (outputDir, guideID), 'w')
        output = []
        i = 0

        for pairID in guide:
            pair = guide[pairID]

            size = att[(pairID, "PAIR", "PRODUCT_SIZE")]
            ltm = "%.1f" % float(att[(pairID, "LEFT", "TM")])
            rtm = "%.1f" % float(att[(pairID, "RIGHT", "TM")])

            lsq = Seq(att[(pairID, "LEFT", "SEQUENCE")])
            rsq = Seq(att[(pairID, "RIGHT", "SEQUENCE")])

            offTargetPairs = has_Off_targets(pair[0], pair[1], PRIMER_OFF_TARGET_MIN, PRIMER_OFF_TARGET_MIN)
            output.append([pair[0].chrom, pair[0].start, pair[0].end, pair[1].start, pair[1].end, i, pair[0].strand,
                           "%s" % lsq, "%s" % rsq, len(pair[0].offTargets), len(pair[1].offTargets),
                           len(offTargetPairs), ltm, rtm, size])

            i += 1

        json.dump(output, outputFile)
        outputFile.close()


# Used in makePrimersFasta and makePrimersGenome
def make_primer_for_target(guide, outputDir, sequence, seqLenBeforeTarget, primer3options, padding):
    template = """PRIMER_SEQUENCE_ID={PRIMER_SEQUENCE_ID:s}
SEQUENCE_TEMPLATE={SEQUENCE_TEMPLATE:s}
SEQUENCE_TARGET={SEQUENCE_TARGET_START:s},{SEQUENCE_TARGET_LEN:s}
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE={PRIMER_OPT_SIZE:s}
PRIMER_MIN_SIZE={PRIMER_MIN_SIZE:s}
PRIMER_MAX_SIZE={PRIMER_MAX_SIZE:s}
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_PRODUCT_SIZE_RANGE={PRODUCT_SIZE_MIN:s}-{PRODUCT_SIZE_MAX:s}
P3_FILE_FLAG=0
PRIMER_EXPLAIN_FLAG=1
"""

    primConfig = PRIMER3_CONFIG.copy()
    primConfig['PRIMER_SEQUENCE_ID'] = str(guide.ID)
    primConfig['SEQUENCE_TEMPLATE'] = sequence
    primConfig['SEQUENCE_TARGET_START'] = str(seqLenBeforeTarget - padding)
    primConfig['SEQUENCE_TARGET_LEN'] = str(guide.targetSize + (2 * padding))

    primer3InputFile = '%s/%s.primer3Input' % (outputDir, guide.ID)
    f = open(primer3InputFile, 'w')
    f.write(template.format(**primConfig))
    f.write(primer3options)
    f.write("=\n")
    f.close()

    command = "%s < %s 2>> %s/primer3.error" % (CONFIG["PATH"]["PRIMER3"], primer3InputFile, outputDir)
    # sys.stderr.write("%s\n" % command)
    prog = Popen(command, stdout=PIPE, shell=True)
    output = prog.communicate()

    if prog.returncode != 0:
        sys.stderr.write("Running Primer3 failed\n")
        sys.exit(EXIT['PRIMER3_ERROR'])

    return output[0].decode()


# Used in main, He lives once more
def make_primers_fasta(targets, outputDir, flanks, displayFlanks, genome, limitPrintResults, bowtieIndexDir,
                       fastaSequence, primer3options, guidePadding, enzymeCo, minResSiteLen, geneID, maxOffTargets):
    primers = {}
    primerOpt = get_primer_options(primer3options)

    primerFastaFileName = '%s/primers.fa' % outputDir
    primerFastaFile = open(primerFastaFileName, 'w')
    for i in range(min(limitPrintResults-1, len(targets))):
        target = targets[i]
        seq, seqLenBeforeTarget = get_primer_query_sequence_fasta(target, outputDir, flanks, fastaSequence)
        primer3_output = make_primer_for_target(target, outputDir, seq, seqLenBeforeTarget, primerOpt, guidePadding)
        region = "%s:%s-%s" % (target.chrom, max(0, target.start-flanks), min(len(fastaSequence), target.end+flanks))
        target_primers, primerPos = parse_primer3_output(target, region, primer3_output, primerFastaFile)
        primers[target.ID] = target_primers

        # Restriction sites
        restSites = dump_restriction_sites(target, seq, flanks, enzymeCo, outputDir, minResSiteLen)
        # Sequence for visualization of locus
        seq2, seqLenBeforeTarget2 = get_primer_query_sequence_fasta(target, outputDir, displayFlanks, fastaSequence)
        dump_locus_sequence(target, outputDir, seq2, seqLenBeforeTarget2, "+")
        # Genbank file for download
        dump_genbank_file(seq, target, restSites, primerPos, outputDir, geneID, target.start-seqLenBeforeTarget, "+")

    primerFastaFile.close()

    primerResults = run_bowtie_primers(primerFastaFileName, outputDir, genome, bowtieIndexDir, maxOffTargets)
    pair_primers(primers, primerResults, outputDir)


# Used in main, Zombie funky
def make_primers_genome(targets, outputDir, flanks, display_seq_len, genome, limitPrintResults, bowtieIndexDir, twoBitToFaIndexDir,
                        primer3options, guidePadding, enzymeCo, minResSiteLen, strand, geneID, maxOffTargets):
    primers = {}

    primerOpt = get_primer_options(primer3options)

    # RUN PRIMER3 ON TARGET SITES AND CREATE FASTA FILE OF PRIMERS FOR BOWTIE
    primerFastaFileName = '%s/primers.fa' % outputDir
    primerFastaFile = open(primerFastaFileName, 'w')
    for i in range(min(limitPrintResults-1, len(targets))):
        target = targets[i]
        seq, seqLenBeforeTarget = get_primer_query_sequence_2bit(
            target, outputDir, flanks, genome, twoBitToFaIndexDir, strand)
        primer3_output = make_primer_for_target(target, outputDir, seq, seqLenBeforeTarget, primerOpt, guidePadding)
        region = "%s:%s-%s" % (target.chrom, max(0, target.start-flanks), target.end+flanks)
        target_primers, primerPos = parse_primer3_output(target, region, primer3_output, primerFastaFile)
        primers[target.ID] = target_primers

        # Restriction sites
        restSites = dump_restriction_sites(target, seq, flanks, enzymeCo, outputDir, minResSiteLen)
        # Sequence for visualization of locus
        seq2, seqLenBeforeTarget2 = get_primer_query_sequence_2bit(
            target, outputDir, display_seq_len, genome, twoBitToFaIndexDir, strand)
        dump_locus_sequence(target, outputDir, seq2, seqLenBeforeTarget2, strand)
        # Genbank file for download
        dump_genbank_file(seq, target, restSites, primerPos, outputDir, geneID, target.start-seqLenBeforeTarget, strand)

    primerFastaFile.close()

    primerResults = run_bowtie_primers(primerFastaFileName, outputDir, genome, bowtieIndexDir, maxOffTargets)
    pair_primers(primers, primerResults, outputDir)


# Used in main and runbowtiePrimer
def parse_bowtie(guideClass, bowtieResultsFile, checkMismatch, scoreGC, scoreSelfComp,
                 backbone, replace5prime, maxOffTargets, countMM, PAM, mode, scoringMethod=None,
                 genome=None, gene=None, isoform=None, gene_isoforms=None):
    """ Parses bowtie hits and build list of guides"""
    logging.info("Parsing bowtie file '%s'." % bowtieResultsFile)

    curr_guide = None
    guide_list = []

    if os.stat(bowtieResultsFile).st_size == 0:  # file is empty
        return guide_list

    sam = pandas.read_csv(bowtieResultsFile, sep='\t', names=list(range(14)),
                          header=None, index_col=False,
                          dtype={0: str, 1: int, 2: str, 3: int, 4: int, 5: str, 6: str, 7: int,
                                 8: int, 9: str, 10: str, 11: str, 12: str, 13: str, 14: str})
    sam_name = sam.iloc[:, 0].value_counts()
    sam_name = sam_name >= maxOffTargets
    if mode:  # Cas9, Cpf1, Nickase and not TALEN
        sam[14] = sam[0].str[-(len(PAM) + 1):]
        sam[0] = sam[0].str[:-(len(PAM) + 1)]
        sam_name = sam.groupby(0).apply(lambda x, m=maxOffTargets: any(x.iloc[:, 14].value_counts() >= m))
        sam = sam.drop([14], axis=1)

        sam = sam.groupby([0, 1, 2, 3]).apply(  # remove duplicates
            lambda x: x.sort_values(by=11).iloc[0])
        sam.rename(columns={0: "name", 11: "mm", 1: "str", 2: "chr", 3: "loc"}, inplace=True)
        sam = sam.sort_values(by=["name", "mm", "str", "chr", "loc"])
        sam = sam.reset_index(drop=True)

    for idx, row in sam.iterrows():
        line = list(row)
        if line[12] != line[12]:
            line = line[:-2]

        #  Encountered a new guide RNA (not a new hit for the same guide)
        elements = line[0].split(":")  # removes from name 5' and 3' tails
        name = ":".join(elements[0:3])
        is_kmaxed = sam_name[line[0]]
        line[0] = ":".join(elements[0:6])
        if len(elements) == 7 and line[1] == 16:
            elements[6] = str(Seq(elements[6]).reverse_complement())
        if curr_guide is None or name != curr_guide.name:
            curr_guide = guideClass(line[0], line[1], len(line[9]),
                                    elements[6] if len(elements) == 7 else line[9], scoreGC, scoreSelfComp,
                                    backbone, PAM, replace5prime, scoringMethod,
                                    genome, gene, isoform, gene_isoforms,
                                    isKmaxed=is_kmaxed)
            guide_list.append(curr_guide)

        # Adds hit to off-target list of current guide.
        curr_guide.addOffTarget(Hit(line), checkMismatch, maxOffTargets, countMM)

    logging.debug("Parsed %d guides from bowtie file '%s'." % (len(guide_list), bowtieResultsFile))

    return guide_list


__all__ = ["make_primers_genome", "make_primers_fasta", "getAllowedFivePrimes"]
