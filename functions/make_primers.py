import json
import logging
import os
import re
import sys
from operator import itemgetter
from subprocess import Popen, PIPE

import pandas
from Bio import SeqIO
from Bio.Restriction import Analysis, RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

import config
from classes.Guide import Guide
from classes.Hit import Hit
from constants import PRIMER3_CONFIG, EXIT, PRIMER_OFF_TARGET_MIN


# Used in makePrimersFasta and makePrimersGenome
def parse_primer3_output(target, region, primer3output, primer_fasta_file):
    pos_pattern = re.compile(r'PRIMER_(\w+)_(\d+)')
    att_pattern = re.compile(r'PRIMER_(\w+)_(\d+)_(\w+)')
    primers = {}
    primer_pos = {}

    for line in primer3output.split("\n"):
        if line[0] == "=":
            break

        label, value = line.split("=")
        m = att_pattern.match(label)
        if m:
            primers[(m.group(2), m.group(1), m.group(3))] = value
        else:
            m = pos_pattern.match(label)

            if m:
                position, length = value.split(",")

                s, e = int(position), int(position) + int(length)
                if m.group(1) == "RIGHT":
                    s, e = int(position) - int(length) + 1, int(position) + 1
                primer_pos[label] = [s, e]

                primer_fasta_file.write(">%s_%s_%s:%s_%s-%s\n%s\n" % (
                    target.id, m.group(2), m.group(1), region, s, e, primers[(m.group(2), m.group(1), "SEQUENCE")]))

    return primers, primer_pos


# Used in makePrimersFasta and makePrimersGenome
def get_primer_options(options):
    # Parse primer3 options. Update config if known option, otherwise append to primer3 input file
    primer_opt = ""

    if options:
        for opt in options.split(","):
            key, value = opt.split("=")
            if key in PRIMER3_CONFIG:
                PRIMER3_CONFIG[key] = value
            else:
                primer_opt += opt + "\n"

    return primer_opt


# Used in makePrimersFasta
def get_primer_query_sequence_fasta(target, output_dir, flank, fasta_sequence):
    s = target.start - flank
    e = target.end + flank
    seq_len_before_target = flank

    if s < 0:
        seq_len_before_target -= abs(s)
        s = 0

    if e > len(fasta_sequence):
        e = len(fasta_sequence)

    return fasta_sequence[s:e], seq_len_before_target


# Used in makePrimerGnome
def get_primer_query_sequence_2bit(target, output_dir, flank, genome, twobittofa_index_dir, strand):
    s = target.start - flank
    seq_len_before_target = flank

    if s < 0:
        seq_len_before_target -= abs(s)
        s = 0

    prog = Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2>> %s/twoBitToFa.err" % (
        config.path("TWOBITTOFA"), target.chrom, s, target.end + flank, twobittofa_index_dir, genome, output_dir),
                 stdout=PIPE, shell=True)
    output = prog.communicate()

    if prog.returncode != 0:
        sys.stderr.write("Running twoBitToFa failed\n")
        sys.exit(EXIT['TWOBITTOFA_ERROR'])

    output = output[0].decode().split("\n")
    del (output[0])
    seq = "".join(output)
    return seq, seq_len_before_target


# Used in makePrimersFasta and makePrimersGenome
def run_bowtie_primers(primer_fasta_file_name, output_dir, genome, bowtie_index_dir, max_off_targets):
    command = "%s -v 0 --best --sam-nohead -k 10 %s/%s -f %s -S %s/primer_results.sam 2> %s/bowtie_primers.err" % (
        config.path("BOWTIE"), bowtie_index_dir, genome, primer_fasta_file_name, output_dir, output_dir)
    prog = Popen(command, shell=True)
    prog.wait()

    if prog.returncode != 0:
        sys.stderr.write("Running bowtie on primers failed\n")
        sys.exit(EXIT['BOWTIE_PRIMER_ERROR'])

    return parse_bowtie(Guide, "%s/primer_results.sam" % output_dir, False, False, False, None, None,
                        max_off_targets, None, None, False, None, None)


# Used in Nickase, Pair, dump_restriction_sites
def find_restriction_sites(sequence, enzyme_company, min_size=1):
    # Take spacer_seq as DNA input for restriction site search
    my_seq = Seq(sequence)

    # Restricts enzyme possibilities to NEB enzymes. Can ultimately change to any supplier.
    rb = RestrictionBatch(first=[], suppliers=[enzyme_company])

    # Filter binding sites shorter than given length
    rb = filter(lambda x: len(x) > min_size, rb)

    # Determine which restriction enzymes cut in the sequence provided
    analyze = Analysis(rb, my_seq)
    return analyze.with_sites()


# Used in makePrimersFasta and makePrimersGenome
def dump_restriction_sites(target, seq, flanks, enzyme_co, output_dir, min_res_site_len):
    sites = find_restriction_sites(seq, enzyme_co, min_res_site_len)
    out = [map(lambda x: [str(enzyme), x + target.start - flanks, enzyme.size], sites[enzyme]) for enzyme in sites]
    out = [item for sublist in out for item in sublist]
    out = sorted(out, key=itemgetter(1))

    # Assign tier to avoid overlaps
    site_count = {}
    tiers = [0] * 23
    for site in out:
        tier = 0

        # count number of sites for each enzyme
        if not site[0] in site_count:
            site_count[site[0]] = 0
        site_count[site[0]] += 1

        for j in range(len(tiers)):
            if site[1] > tiers[j]:
                tier = j
                tiers[j] = site[1] + site[2]
                break
        site.append(tier)

    # Assign colors depending on uniqueness
    for site in out:
        if site_count[site[0]] == 1:
            site.append("green")
        else:
            site.append("red")

    output_file = open("%s/restriction_%s.json" % (output_dir, target.id), 'w')
    json.dump(out, output_file)
    output_file.close()

    return sites


# Used in make_primers_fasta and make_primers_Genome
def dump_locus_sequence(target, output_dir, seq, seq_len_before_target, strand):
    if strand == "-":
        seq = str(Seq(seq).complement())
    out = [[target.start - seq_len_before_target, target.end, seq]]
    output_file = open("%s/locusSeq_%s.json" % (output_dir, target.id), 'w')
    json.dump(out, output_file)
    output_file.close()


# Used in makePrimersFasta and makePrimersGenome
def dump_genbank_file(seq, target, rest_sites, primers, output_dir, gene_id, loci_start, strand):
    name = "%s, locus %s" % (gene_id, target.id)
    desc = "CHOPCHOP prediction for gene %s, target %s" % (gene_id, target.id)
    annotation = {"organism": "Danio rerio", "Target location": "chrT:1-20"}

    # Genbank file
    genbank_file = open('%s/%s_%s.gb' % (output_dir, gene_id, target.id), 'w')
    record = SeqRecord(Seq(seq), description=desc, name="CHOPCHOP", id=name)
    record.annotation = annotation

    if target.strand == "+":
        ts = 1
    else:
        ts = -1

    record.features.append(
        SeqFeature(FeatureLocation(target.start - loci_start - 1, target.end - loci_start - 1, strand=ts),
                   type="Target"))

    for primer in primers:
        record.features.append(SeqFeature(FeatureLocation(primers[primer][0], primers[primer][1]), type=primer))

    if strand == "-":
        record = record.reverse_complement()

    # TODO: Have an actual bioinformatics person look at this.
    record.annotations["molecule_type"] = "DNA"

    SeqIO.write(record, genbank_file, "genbank")
    genbank_file.close()

    pass


def has_off_targets(tale1, tale2, off_target_min, off_target_max):
    """ Returns the number of off-targets for a pair of TALENs (10-24bp apart) """

    off_target_pairs = []

    # Calls sort function to sort off-targets by chromosome and chromosome position.
    # Bowtie ranks them according to quality of hit
    tale1.sort_off_targets()
    tale2.sort_off_targets()

    # TODO: Eivind to write this code properly. Include a way to step backwards, so as not to miss any hits.
    # Need to make a queue..?
    for i in range(len(tale1.off_targets)):
        hit1 = tale1.off_targets[i]

        for j in range(len(tale2.off_targets)):
            hit2 = tale2.off_targets[j]

            # Determines whether 2 tales are on the same chromosome and 10-24 bp apart.
            if hit2.chrom == hit1.chrom and off_target_min <= abs(hit2.start - hit1.start) <= off_target_max:
                off_target_pairs.append([hit1, hit2])

    return off_target_pairs


# Used in makePrimersFasta and makePrimersGenome
def pair_primers(primer_attributes, primer_list, output_dir):
    primers = {}

    for primer in primer_list:
        guide, primer_pair_id, side = primer.id.split("_")

        s = 0
        if side == "RIGHT":
            s = 1

        if guide not in primers:
            primers[guide] = {}

        if primer_pair_id not in primers[guide]:
            primers[guide][primer_pair_id] = [None, None]

        primers[guide][primer_pair_id][s] = primer

    for guideID in primers:
        guide = primers[guideID]

        att = primer_attributes[int(guideID)]

        output_file = open("%s/primer_%s.json" % (output_dir, guideID), 'w')
        output = []
        i = 0

        for pairID in guide:
            pair = guide[pairID]

            size = att[(pairID, "PAIR", "PRODUCT_SIZE")]
            ltm = "%.1f" % float(att[(pairID, "LEFT", "TM")])
            rtm = "%.1f" % float(att[(pairID, "RIGHT", "TM")])

            lsq = Seq(att[(pairID, "LEFT", "SEQUENCE")])
            rsq = Seq(att[(pairID, "RIGHT", "SEQUENCE")])

            off_target_pairs = has_off_targets(pair[0], pair[1], PRIMER_OFF_TARGET_MIN, PRIMER_OFF_TARGET_MIN)
            output.append([pair[0].chrom, pair[0].start, pair[0].end, pair[1].start, pair[1].end, i, pair[0].strand,
                           "%s" % lsq, "%s" % rsq, len(pair[0].off_targets), len(pair[1].off_targets),
                           len(off_target_pairs), ltm, rtm, size])

            i += 1

        json.dump(output, output_file)
        output_file.close()


# Used in makePrimersFasta and makePrimersGenome
def make_primer_for_target(guide, output_dir, sequence, seq_len_before_target, primer3_options, padding):
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

    prim_config = PRIMER3_CONFIG.copy()
    prim_config['PRIMER_SEQUENCE_ID'] = str(guide.id)
    prim_config['SEQUENCE_TEMPLATE'] = sequence
    prim_config['SEQUENCE_TARGET_START'] = str(seq_len_before_target - padding)
    prim_config['SEQUENCE_TARGET_LEN'] = str(guide.target_size + (2 * padding))

    primer_3_input_file = '%s/%s.primer3Input' % (output_dir, guide.id)
    f = open(primer_3_input_file, 'w')
    f.write(template.format(**prim_config))
    f.write(primer3_options)
    f.write("=\n")
    f.close()

    command = "%s < %s 2>> %s/primer3.error" % (config.path("PRIMER3"), primer_3_input_file, output_dir)
    # sys.stderr.write("%s\n" % command)
    prog = Popen(command, stdout=PIPE, shell=True)
    output = prog.communicate()

    if prog.returncode != 0:
        sys.stderr.write("Running Primer3 failed\n")
        sys.exit(EXIT['PRIMER3_ERROR'])

    return output[0].decode()


# Used in main, He lives once more
def make_primers_fasta(targets, output_dir, flanks, display_flanks, genome, limit_print_results, bowtie_index_dir,
                       fasta_sequence, primer3_options, guide_padding, enzyme_co, min_res_site_len, gene_id,
                       max_off_targets):
    primers = {}
    primer_opt = get_primer_options(primer3_options)

    primer_fasta_file_name = '%s/primers.fa' % output_dir
    primer_fasta_file = open(primer_fasta_file_name, 'w')
    for i in range(min(limit_print_results - 1, len(targets))):
        target = targets[i]
        seq, seq_len_before_target = get_primer_query_sequence_fasta(target, output_dir, flanks, fasta_sequence)
        primer3_output = make_primer_for_target(target, output_dir, seq, seq_len_before_target, primer_opt,
                                                guide_padding)
        region = "%s:%s-%s" % (target.chrom, max(0, target.start - flanks),
                               min(len(fasta_sequence), target.end + flanks))
        target_primers, primer_pos = parse_primer3_output(target, region, primer3_output, primer_fasta_file)
        primers[target.id] = target_primers

        # Restriction sites
        rest_sites = dump_restriction_sites(target, seq, flanks, enzyme_co, output_dir, min_res_site_len)
        # Sequence for visualization of locus
        seq2, seq_len_before_target2 = get_primer_query_sequence_fasta(target, output_dir, display_flanks,
                                                                       fasta_sequence)
        dump_locus_sequence(target, output_dir, seq2, seq_len_before_target2, "+")
        # Genbank file for download
        dump_genbank_file(seq, target, rest_sites, primer_pos, output_dir, gene_id,
                          target.start - seq_len_before_target, "+")

    primer_fasta_file.close()

    primer_results = run_bowtie_primers(primer_fasta_file_name, output_dir, genome, bowtie_index_dir, max_off_targets)
    pair_primers(primers, primer_results, output_dir)


# Used in main, Zombie funky
def make_primers_genome(targets, output_dir, flanks, display_seq_len, genome, limit_print_results, bowtie_index_dir,
                        twobittofa_index_dir, primer3options, guide_padding, enzyme_co, min_res_site_len, strand,
                        gene_id, max_off_targets):
    primers = {}

    primer_opt = get_primer_options(primer3options)

    # RUN PRIMER3 ON TARGET SITES AND CREATE FASTA FILE OF PRIMERS FOR BOWTIE
    primer_fasta_file_name = '%s/primers.fa' % output_dir
    primer_fasta_file = open(primer_fasta_file_name, 'w')
    for i in range(min(limit_print_results - 1, len(targets))):
        target = targets[i]
        seq, seq_len_before_target = get_primer_query_sequence_2bit(
            target, output_dir, flanks, genome, twobittofa_index_dir, strand)
        primer3_output = make_primer_for_target(target, output_dir, seq, seq_len_before_target, primer_opt,
                                                guide_padding)
        region = "%s:%s-%s" % (target.chrom, max(0, target.start - flanks), target.end + flanks)
        target_primers, primer_pos = parse_primer3_output(target, region, primer3_output, primer_fasta_file)
        primers[target.id] = target_primers

        # Restriction sites
        rest_sites = dump_restriction_sites(target, seq, flanks, enzyme_co, output_dir, min_res_site_len)
        # Sequence for visualization of locus
        seq2, seq_len_before_target2 = get_primer_query_sequence_2bit(
            target, output_dir, display_seq_len, genome, twobittofa_index_dir, strand)
        dump_locus_sequence(target, output_dir, seq2, seq_len_before_target2, strand)
        # Genbank file for download
        dump_genbank_file(seq, target, rest_sites, primer_pos, output_dir, gene_id,
                          target.start - seq_len_before_target, strand)

    primer_fasta_file.close()

    primer_results = run_bowtie_primers(primer_fasta_file_name, output_dir, genome, bowtie_index_dir, max_off_targets)
    pair_primers(primers, primer_results, output_dir)


# Used in main and runbowtiePrimer
def parse_bowtie(guide_class, bowtie_results_file, check_mismatch, score_gc, score_self_comp,
                 backbone, replace_5prime, max_off_targets, count_mm, pam, mode, scoring_method=None,
                 genome=None, gene=None, isoform=None, gene_isoforms=None):
    """ Parses bowtie hits and build list of guides"""
    logging.info("Parsing bowtie file '%s'." % bowtie_results_file)

    curr_guide = None
    guide_list = []

    if os.stat(bowtie_results_file).st_size == 0:  # file is empty
        return guide_list

    sam = pandas.read_csv(bowtie_results_file, sep='\t', names=list(range(14)),
                          header=None, index_col=False,
                          dtype={0: str, 1: int, 2: str, 3: int, 4: int, 5: str, 6: str, 7: int,
                                 8: int, 9: str, 10: str, 11: str, 12: str, 13: str, 14: str})
    sam_name = sam.iloc[:, 0].value_counts()
    sam_name = sam_name >= max_off_targets
    if mode:  # Cas9, Cpf1, Nickase and not TALEN
        sam[14] = sam[0].str[-(len(pam) + 1):]
        sam[0] = sam[0].str[:-(len(pam) + 1)]
        sam_name = sam.groupby(0).apply(lambda x, m=max_off_targets: any(x.iloc[:, 14].value_counts() >= m))
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
            curr_guide = guide_class(line[0], line[1], len(line[9]),
                                     elements[6] if len(elements) == 7 else line[9], score_gc, score_self_comp,
                                     backbone, pam, replace_5prime, scoring_method,
                                     genome, gene, isoform, gene_isoforms,
                                     is_kmaxed=is_kmaxed)
            guide_list.append(curr_guide)

        # Adds hit to off-target list of current guide.
        curr_guide.add_off_target(Hit(line), check_mismatch, max_off_targets, count_mm)

    logging.debug("Parsed %d guides from bowtie file '%s'." % (len(guide_list), bowtie_results_file))

    return guide_list


__all__ = ["make_primers_genome", "make_primers_fasta", "parse_bowtie", "find_restriction_sites", "has_off_targets"]
