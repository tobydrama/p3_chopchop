import csv
import logging
import re
import sys
from subprocess import Popen, PIPE

import numpy
from Bio.Seq import Seq

import config
from constants import EXIT, TARGET_MAX


# Used in ParseTargets
def bins(x):  # from ranges to bins
    x = list(x)
    x.sort()
    x = numpy.array(x)
    if x.size == 1:
        return x, x
    dx, = numpy.nonzero(numpy.diff(x) > 1)
    starts = numpy.append(x[0], x[dx + 1])
    ends = numpy.append(x[dx], x[-1])
    return starts, ends


# Used in ParseTargets
def get_isoforms(gene, table_file):
    gene_isoforms = set()
    table_r = open(table_file, 'r')
    table_reader = csv.DictReader(table_r, delimiter='\t', quoting=csv.QUOTE_NONE)
    for row in table_reader:
        if row['name2'] == gene:
            gene_isoforms.add(row['name'])
    table_r.close()
    return gene_isoforms


# Used in parse_targets
def filter_repeating_names(tx_info, filter_names=["fix", "random", "alt"]):
    # if more isoforms have exact same name filter the ones
    # with "alt", "fix", "random" in chr names
    # then take the first one
    seen = []
    same_name_tx = []
    is_special = []
    for x in tx_info:
        if str(x[3]) not in seen:
            seen.append(str(x[3]))
            same_name_tx.append([x])
            is_special.append([any(fn in str(x[0]) for fn in filter_names)])
        else:
            idx = seen.index(str(x[3]))
            same_name_tx[idx].append(x)
            is_special[idx].append(any(fn in str(x[0]) for fn in filter_names))

    tx_info_ = []
    for i, tx in enumerate(same_name_tx):
        if any(is_special[i]) and sum(is_special[i]) < len(is_special[i]):
            idx = [i for i, x in enumerate(is_special[i]) if not x]
            tx_info_.append(tx[idx[0]])
        else:
            tx_info_.append(tx[0])

    return tx_info_


# Used in subsetExon
def hyphen_range(s):
    """ Takes a range in form of "a-b" and generate a list of numbers between a and b inclusive.
    Also accepts comma separated ranges like "a-b,c-d,f" will build a list which will include
    Numbers from a to b, a to d and f"""

    s = "".join(s.split())  # removes white space
    r = set()

    for x in s.split(','):
        t = x.split('-')
        if len(t) not in [1, 2]:
            raise SyntaxError("Range is not properly formatted: " + s)
        if len(t) == 1:
            r.add(int(t[0]))
        else:
            r.update(set(range(int(t[0]), int(t[1]) + 1)))

    l = list(r)
    l.sort()

    return l


# Used in parse_targets
def subset_exons(exons, targets):
    if exons:
        indices = hyphen_range(exons)
        for index in indices:
            if int(index) > len(targets):
                sys.stderr.write("That exon does not exist\n")
                sys.exit(EXIT['PYTHON_ERROR'])
        targets = [targets[int(i) - 1] for i in indices]  # indices is a list of exon numbers -1 e.g. exon 2 is [1]
    return targets


# Used in parse_targets
def truncate_to_utr5(cds_start, exons):
    """ Truncates the gene to only target 5' UTR """

    end_exon = 0
    for exon in range(len(exons)):
        if (cds_start > exons[exon][1]) and (cds_start < exons[exon][2]):
            exons[exon][2] = cds_start
            end_exon = exon
            break

    return exons[:end_exon + 1]


# Used in parse_targets
def truncate_to_promoter(strand, exons, ups_bp, down_bp):
    """ Truncates the gene to only target promoter +-bp TSS """

    if strand == "+":
        first_exon = exons[0]
        first_exon[2] = first_exon[1] + down_bp
        first_exon[1] = first_exon[1] - ups_bp
        return [first_exon]
    else:
        first_exon = exons[-1]
        first_exon[1] = first_exon[2] - down_bp
        first_exon[2] = first_exon[2] + ups_bp
        return [first_exon]


# Used in parse_targets
def truncate_to_utr3(cds_end, exons):
    """ Truncates the gene to only target 3' UTR """

    start_exon = 0
    for exon in range(len(exons)):
        if (cds_end > exons[exon][1]) and (cds_end < exons[exon][2]):
            exons[exon][1] = cds_end
            start_exon = exon

    return exons[start_exon:]


# Used in parse_targets
def truncate_to_splice(exons):
    """ Truncates the gene to only target splice sites """

    splice_sites = []
    for ind in range(0, len(exons)):
        splice_sites.append([exons[ind][0], exons[ind][1] - 1, exons[ind][1] + 1])
        splice_sites.append([exons[ind][0], exons[ind][2] - 1, exons[ind][2] + 1])
    # Remove first and last (i.e. transcription start and termination site)
    return splice_sites[1:len(splice_sites) - 1]


# Used in parse_targets
def truncate_to_coding(cds_start, cds_end, exons):
    """ Truncates the gene to only consider the coding region """

    start_exon, end_exon = 0, len(exons) - 1
    # Shortens the coding region to the exons and coordinates between the cds start and cds end
    for exon in range(len(exons)):
        if (cds_start >= exons[exon][1]) and (cds_start <= exons[exon][2]):
            exons[exon][1] = cds_start
            start_exon = exon

        if (cds_end >= exons[exon][1]) and (cds_end <= exons[exon][2]):
            # replace the end with the cds end
            exons[exon][2] = cds_end
            end_exon = exon

    if start_exon > end_exon:
        start_exon, end_exon = end_exon, start_exon

    # Shorten list to include exons from cds start to end
    return exons[start_exon:(end_exon + 1)]


# Used in parse_targets
def gene_to_coord_db(gene, organism, db):
    """ Gets genomic coordinates for a gene from a database """

    # Try refseq first
    lines = db.execute(
        "SELECT chrom, exonStarts, exonEnds, r.name, cdsStart, cdsEnd, strand, txStart, txEnd FROM organism o, refGene"
        " r WHERE o.assembly='%s' AND o.organism_id=r.organism_id AND (r.name='%s' OR r.name2='%s')" % (
            organism, gene, gene))

    # Then Ensembl
    if lines == 0:
        lines = db.execute(
            "SELECT chrom, exonStarts, exonEnds, r.name, cdsStart, cdsEnd, strand, txStart, txEnd FROM organism o,"
            " ensGene r LEFT OUTER JOIN ensemblToGeneName g ON r.name=g.name WHERE o.assembly='%s'"
            " AND o.organism_id=r.organism_id AND  (r.name='%s' OR r.name2='%s' OR g.value='%s')" % (
                organism, gene, gene, gene))

    # Then the general genePred table
    if lines == 0:
        lines = db.execute(
            "SELECT chrom, exonStarts, exonEnds, r.name, cdsStart, cdsEnd, strand, txStart, txEnd FROM organism o, "
            "gpGene r WHERE o.assembly='%s' AND o.organism_id=r.organism_id AND (r.name='%s' OR r.name2='%s')" % (
                organism, gene, gene))

    # Then wormbase. FIX: NO HARDCODED ASSEMBLY!!!
    if organism == "ce6" and lines == 0:
        lines = db.execute(
            "SELECT chrom, exonStarts, exonEnds, name, cdsStart, cdsEnd, strand, txStart, txEnd FROM sangerGene WHERE"
            " (name='%s' OR proteinID='%s')" % (
                gene, gene))

    # Then flybase. FIX: NO HARDCODED ASSEMBLY!!!
    if organism == "dm3" and lines == 0:
        lines = db.execute(
            "SELECT chrom, exonStarts, exonEnds, name, cdsStart, cdsEnd, strand, txStart, txEnd FROM flyBaseGene WHERE"
            " name='%s'" % (
                gene))

    if lines == 0:
        sys.stderr.write(
            "The gene name %s was not found in the gene sets for assembly %s. Consider trying an alternative ID"
            " (see the instruction page for supported gene identifiers) or using genomic coordinates. If you believe"
            " this type of ID should be supported for your organism contact us and we will do our best to support it."
            " \n" % (
                gene, organism))
        sys.exit(EXIT['GENE_ERROR'])

    tx_info = []
    for i in range(lines):
        tx_info.append(db.fetchone())

    return tx_info


# Used in parse_targets
def gene_to_coord_file(gene_in, table_file):
    """ Extracts coordinates of genomic regions to parse for suitable guide binding sites """
    table_r = open(table_file, 'r')

    table_reader = csv.DictReader(table_r, delimiter='\t', quoting=csv.QUOTE_NONE)
    tx_info = []
    gene = None
    # Look in genome table for gene of question
    for row in table_reader:
        if row['name'] == gene_in or row['name2'] == gene_in or row['name'] == gene_in.upper() \
                or row['name2'] == gene_in.upper():
            tx_info.append([row['chrom'], row['exonStarts'], row['exonEnds'], row['name'],
                            row['cdsStart'], row['cdsEnd'], row['strand'],
                            row['txStart'], row['txEnd']])
            gene = row['name2']
    table_r.close()

    if len(tx_info) == 0:
        sys.stderr.write("The gene name %s does not exist in file %s. Please try again.\n" % (gene_in, table_file))
        sys.exit(EXIT['GENE_ERROR'])

    return gene, tx_info


def coordinate_search(is_coordinate, target_string, pattern, target_size, vis_coords, targets, pad_size, make_vis):
    if config.isoforms:
        sys.stderr.write("--isoforms is not working with coordinate search.\n")
        sys.exit(EXIT['ISOFORMS_ERROR'])

    chrom = is_coordinate.group(2)
    vis_coords.append({"exons": [], "ATG": [], "name": chrom})

    for target in target_string.split(";"):
        m = pattern.match(target)
        if m:
            if m.group(2) is not None and chrom != m.group(2):
                sys.stderr.write(
                    "Can't target regions on separate chromosomes (%s != %s).\n" % (chrom, m.group(2)))
                sys.exit(EXIT['GENE_ERROR'])

            start_pos = m.group(3)
            end_pos = m.group(4)
            start_pos = int(start_pos.replace(",", "").replace(".", ""))
            end_pos = int(end_pos.replace(",", "").replace(".", ""))
            target_size += end_pos - start_pos + 1

            if start_pos >= end_pos:
                sys.stderr.write(
                    "Start position (%s) must be smaller than end position (%s)\n" % (start_pos, end_pos))
                sys.exit(EXIT['GENE_ERROR'])

            targets.append("%s:%s-%s" % (chrom, max(0, start_pos - pad_size), end_pos + pad_size))
            if make_vis:
                vis_coords[0]["exons"].append([chrom, start_pos, end_pos, 0, True, "+"])
        else:
            sys.stderr.write("Unknown format: %s\n" % target)
            sys.exit(EXIT['GENE_ERROR'])
    return target_size, vis_coords


def make_vis_coords(starts_v, ends_v, tx, tx_vis, index_dir, genome, output_dir, vis_coords):
    intron_size = [int(starts_v[x + 1]) - int(ends_v[x]) for x in range(len(starts_v) - 1)]
    intron_size.append(0)
    # tx_vis exons are [chr, start, end, intron_size, isIntron, strand]
    for e in range(len(starts_v)):
        if ends_v[e] <= tx[4] or starts_v[e] >= tx[5]:
            tx_vis["exons"].append([tx[0], starts_v[e], ends_v[e], intron_size[e], True, tx[6]])
        else:
            if starts_v[e] < tx[4] < ends_v[e]:
                tx_vis["exons"].append([tx[0], starts_v[e], tx[4], 0, True, tx[6]])
                starts_v[e] = tx[4]

            if starts_v[e] < tx[5] < ends_v[e]:
                tx_vis["exons"].append([tx[0], tx[5], ends_v[e], intron_size[e], True, tx[6]])
                ends_v[e] = tx[5]
                intron_size[e] = 0

            tx_vis["exons"].append([tx[0], starts_v[e], ends_v[e], intron_size[e], False, tx[6]])

    tx_vis["exons"].sort(key=lambda x: x[1])  # sort on starts
    # ATG locations
    logging.debug(f"Running twoBitToFa on file '{index_dir + '/' + genome}.2bit'")

    prog = Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2> %s/twoBitToFa.err" % (
        config.path("TWOBITTOFA"), tx[0], int(tx[4]) + 1, int(tx[5]) + 1, index_dir,
        genome, output_dir), stdout=PIPE, shell=True)
    iso_seq = prog.communicate()
    if prog.returncode != 0:
        sys.stderr.write("Running twoBitToFa when searching isoform sequence failed\n")
        sys.exit(EXIT['TWOBITTOFA_ERROR'])

    iso_seq = iso_seq[0].decode()
    iso_seq = iso_seq.split("\n")
    iso_seq = Seq(''.join(iso_seq[1:]).upper())
    # splicing
    iso_seq_spl = ""
    for e in tx_vis["exons"]:
        if not e[4]:
            iso_seq_spl += iso_seq[(e[1] - tx[4]):(e[2] - tx[4])]
    atg = "ATG" if tx[6] != "-" else "CAT"
    tx_atg = [m.start() for m in re.finditer(atg, str(iso_seq_spl)) if m.start() % 3 == 0]
    tx_atg.sort()
    for atg1 in tx_atg:  # every ATG as 3 x 1bp as they can span across two exons...
        atg2 = atg1 + 1
        atg3 = atg1 + 2
        shift_atg1, shift_atg2, shift_atg3, exon_len = 0, 0, 0, 0
        for e in tx_vis["exons"]:  # exons are sorted
            if not e[4]:
                exon_len += (e[2] - e[1])
                if atg1 > exon_len:
                    shift_atg1 += e[3]
                if atg2 > exon_len:
                    shift_atg2 += e[3]
                if atg3 > exon_len:
                    shift_atg3 += e[3]
        tx_vis["ATG"].extend([atg1 + shift_atg1 + tx[4], atg2 + shift_atg2 + tx[4],
                              atg3 + shift_atg3 + tx[4]])

    vis_coords.append(tx_vis)
    return vis_coords


def truncate_to_region(target_region, tx, coords, ups_bp, down_bp):
    if target_region == "CODING":
        coords = truncate_to_coding(tx[4], tx[5], coords)
    elif target_region == "UTR5":
        if tx[6] == "+":
            coords = truncate_to_utr5(tx[4], coords)
        else:
            coords = truncate_to_utr3(tx[5], coords)
    elif target_region == "PROMOTER":
        coords = truncate_to_promoter(tx[6], coords, ups_bp, down_bp)
    elif target_region == "UTR3":
        if tx[6] == "+":
            coords = truncate_to_utr3(tx[5], coords)
        else:
            coords = truncate_to_utr5(tx[4], coords)
    elif target_region == "SPLICE":
        coords = truncate_to_splice(coords)
    elif target_region != "WHOLE":
        sys.stderr.write("Unknown region: %s\n" % target_region)
        sys.exit(EXIT['PYTHON_ERROR'])
    return coords


def compute_intersection_union_all_exions(tx_info, tx, coords, targets, use_union, guide_len, isoforms):
    if tx_info[0][3] == tx[3]:  # if this is first of the isoforms
        for x in coords:
            targets.extend(range(x[1], x[2] + 1))
        targets = set(targets)
    else:
        if not use_union:
            targets_ = []
            for x in coords:
                targets_.extend(range(x[1], x[2] + 1))

            if len(targets_) >= guide_len:  # cover cases where some transcripts provide short or none bp
                targets &= set(targets_)

            if len(targets) < guide_len:
                sys.stderr.write(
                    "Computing intersection over specified isoforms resulted in lack of targets." +
                    " Consider either using specific isoform as input: " + ', '.join(isoforms) +
                    " or using --consensusUnion to compute union instead of intersection " +
                    "of your isoforms (on the website you can find it in " +
                    "Options -> General -> Isoform consensus determined by -> Union.")
                sys.exit(EXIT['GENE_ERROR'])
        else:
            targets_ = []
            for x in coords:
                targets_.extend(range(x[1], x[2] + 1))
            targets |= set(targets_)
    return targets


# Used in main
def parse_targets(target_string, genome, use_db, data, pad_size, target_region, exon_subset, ups_bp, down_bp,
                  index_dir, output_dir, use_union, make_vis, guide_len):
    targets = []
    vis_coords = []
    target_strand = "+"
    target_size = 0
    gene, isoform, gene_isoforms = (None, None, set())

    pattern = re.compile(r"(([.\w]+):)?([.,\d]+)-([.,\d]+)")
    is_coordinate = pattern.match(str(target_string))

    if is_coordinate:
        target_size, vis_coords = coordinate_search(is_coordinate, target_string, pattern, target_size, vis_coords,
                                                    targets, pad_size, make_vis)

    else:
        if use_db:
            if config.isoforms:
                sys.stderr.write("--isoforms is not working with database search.\n")
                sys.exit(EXIT['ISOFORMS_ERROR'])
            tx_info = gene_to_coord_db(target_string, genome, data)
            tx_info = filter_repeating_names(tx_info)
        else:
            gene, tx_info = gene_to_coord_file(target_string, data)
            tx_info = filter_repeating_names(tx_info)
            isoform = "union" if use_union else "intersection"
            gene_isoforms = set([str(x[3]) for x in tx_info])
            if target_string in gene_isoforms:
                isoform = target_string
                gene_isoforms = get_isoforms(gene, data)

        target_chr = set([x[0] for x in tx_info])
        target_strand = set([x[6] for x in tx_info])
        isoforms = [str(x[3]) for x in tx_info]
        if len(target_strand) > 1 or len(target_chr) > 1:
            sys.stderr.write(
                "Specify which isoform you want to target as your query " + str(target_string) +
                " returns many isoforms: " + ', '.join(isoforms) +
                " which are from either inconsistent strands or chromosomes.\n")
            sys.exit(EXIT['GENE_ERROR'])
        else:
            target_strand = list(target_strand)[0]
            target_chr = list(target_chr)[0]

        for tx in tx_info:
            tx = list(tx)
            tx[4] = int(tx[4])
            tx[5] = int(tx[5])
            starts = tx[1].split(",")
            ends = tx[2].split(",")
            del starts[-1]
            del ends[-1]
            starts = list(map(int, starts))
            ends = list(map(int, ends))
            starts_v = starts[:]
            ends_v = ends[:]
            tx_vis = {"exons": [], "ATG": [], "name": tx[3]}

            if make_vis:
                vis_coords = make_vis_coords(starts_v, ends_v, tx, tx_vis, index_dir, genome, output_dir, vis_coords)

            # restrict isoforms
            coords = list(map(lambda x: [tx[0], x[0], x[1]], zip(starts, ends)))
            if tx[6] == "-":
                coords.reverse()
            coords = subset_exons(exon_subset, coords)
            if tx[6] == "-":
                coords.reverse()

            # Truncate to region
            coords = truncate_to_region(target_region, tx, coords, ups_bp, down_bp)

            # filter exons that are too truncated
            coords = [x for x in coords if x[1] < x[2]]
            if not coords:
                if gene_isoforms:
                    gene_isoforms.remove(tx[3])
                if vis_coords:
                    del vis_coords[-1]

            # compute intersection/union on all exons
            targets = compute_intersection_union_all_exions(tx_info, tx, coords, targets, use_union, guide_len,
                                                            isoforms)

        target_size = len(targets)
        if target_size < guide_len:
            sys.stderr.write("Search region is too small. You probably want to specify -t option as WHOLE")
            sys.exit(EXIT['GENE_ERROR'])

        starts, ends = bins(targets)
        if config.isoforms:
            targets = list(map(lambda x: "%s:%s-%s" % (target_chr, x[0], x[1]), zip(starts, ends)))
        else:
            targets = list(
                map(lambda x: "%s:%s-%s" % (target_chr, x[0] - pad_size, x[1] + pad_size), zip(starts, ends)))

    if target_size > TARGET_MAX:
        sys.stderr.write("Search region is too large (%s nt). Maximum search region is %s nt.\n" % (
            target_size, TARGET_MAX))
        sys.exit(EXIT['GENE_ERROR'])

    return targets, vis_coords, target_strand, gene, isoform, gene_isoforms


__all__ = ["parse_targets"]
