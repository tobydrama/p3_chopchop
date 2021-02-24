from Bio.Restriction import Analysis, RestrictionBatch
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from subprocess import Popen, PIPE
from Vars import *
from Bio.SeqRecord import SeqRecord

import scipy.stats as ss
import sys
import numpy
import pandas
import csv
import re
import Hit



#Used in Cas9, Guide, and selfComp
def gccontent(seq):
    gc = 0
    for i in seq:
        if i == 'G' or i == 'g' or i == 'C' or i == 'c':
            gc += 1
    return float(gc) / float(len(seq))


#Used in Nickase, Pair, dump_restriction_sites
def findRestrictionSites(sequence, enzymeCompany, minSize=1):
    # Take spacerSeq as DNA input for restriction site search
    mySeq = Seq(sequence, IUPACAmbiguousDNA())

    # Restricts enzyme possibilities to NEB enzymes. Can ultimately change to any supplier.
    rb = RestrictionBatch(first=[], suppliers=[enzymeCompany])

    # Filter binding sites shorter than given length
    rb = filter(lambda x: len(x) > minSize, rb)

    # Determine which restriction enzymes cut in the sequence provided
    analyze = Analysis(rb, mySeq)
    return analyze.with_sites()


#Used in parseTargets
def truncateToUTR5(cds_start, exons):
    """ Truncates the gene to only target 5' UTR """

    end_exon = 0
    for exon in range(len(exons)):
        if (cds_start > exons[exon][1]) and (cds_start < exons[exon][2]):
            exons[exon][2] = cds_start
            end_exon = exon
            break

    return exons[:end_exon + 1]

#Used in parseTargets
def truncateToPROMOTER(strand, exons, ups_bp, down_bp):
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

    return exons

#Used in parseTargets
def truncateToUTR3(cds_end, exons):
    """ Truncates the gene to only target 3' UTR """

    start_exon = 0
    for exon in range(len(exons)):
        if (cds_end > exons[exon][1]) and (cds_end < exons[exon][2]):
            exons[exon][1] = cds_end
            start_exon = exon

    return exons[start_exon:]

#Used in parseTargets
def truncateToSplice(exons):
    """ Truncates the gene to only target splice sites """

    splice_sites = []
    for ind in range(0, len(exons)):
        splice_sites.append([exons[ind][0], exons[ind][1]-1, exons[ind][1]+1])
        splice_sites.append([exons[ind][0], exons[ind][2]-1, exons[ind][2]+1])
    # Remove first and last (i.e. transcription start and termination site)
    return splice_sites[1:len(splice_sites)-1]

#Used in parseTargets
def truncateToCoding(cds_start, cds_end, exons):
    """ Truncates the gene to only consider the coding region """

    start_exon, end_exon = 0, len(exons)-1
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
    return exons[start_exon:(end_exon+1)]

#Used in parseTargets
def geneToCoord_db(gene, organism, db):
    """ Gets genomic coordinates for a gene from a database """

    # Try refseq first
    lines = db.execute("SELECT chrom, exonStarts, exonEnds, r.name, cdsStart, cdsEnd, strand, txStart, txEnd FROM organism o, refGene r WHERE o.assembly='%s' AND o.organism_id=r.organism_id AND (r.name='%s' OR r.name2='%s')" % (organism, gene, gene))

    # Then Ensembl
    if lines == 0:
        lines = db.execute("SELECT chrom, exonStarts, exonEnds, r.name, cdsStart, cdsEnd, strand, txStart, txEnd FROM organism o, ensGene r LEFT OUTER JOIN ensemblToGeneName g ON r.name=g.name WHERE o.assembly='%s' AND o.organism_id=r.organism_id AND  (r.name='%s' OR r.name2='%s' OR g.value='%s')" % (organism, gene, gene, gene))

    # Then the general genePred table
    if lines == 0:
        lines = db.execute("SELECT chrom, exonStarts, exonEnds, r.name, cdsStart, cdsEnd, strand, txStart, txEnd FROM organism o, gpGene r WHERE o.assembly='%s' AND o.organism_id=r.organism_id AND (r.name='%s' OR r.name2='%s')" % (organism, gene, gene))

    # Then wormbase. FIX: NO HARDCODED ASSEMBLY!!!
    if organism == "ce6" and lines == 0:
        lines = db.execute("SELECT chrom, exonStarts, exonEnds, name, cdsStart, cdsEnd, strand, txStart, txEnd FROM sangerGene WHERE (name='%s' OR proteinID='%s')" % (gene, gene))

    # Then flybase. FIX: NO HARDCODED ASSEMBLY!!!
    if organism == "dm3" and lines == 0:
        lines = db.execute("SELECT chrom, exonStarts, exonEnds, name, cdsStart, cdsEnd, strand, txStart, txEnd FROM flyBaseGene WHERE name='%s'" % (gene))

    if lines == 0:
        sys.stderr.write("The gene name %s was not found in the gene sets for assembly %s. Consider trying an alternative ID (see the instruction page for supported gene identifiers) or using genomic coordinates. If you believe this type of ID should be supported for your organism contact us and we will do our best to support it. \n" % (gene, organism))
        sys.exit(EXIT['GENE_ERROR'])

    txInfo = []
    for i in range(lines):
        txInfo.append(db.fetchone())

    return txInfo

#Used in parseTargets
def geneToCoord_file(gene_in, table_file):
    """ Extracts coordinates of genomic regions to parse for suitable guide binding sites """

    table_r = open(table_file, 'rb')
    tablereader = csv.DictReader(table_r, delimiter='\t', quoting=csv.QUOTE_NONE)

    tx_info = []
    gene = None
    # Look in genome table for gene of question
    for row in tablereader:
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

#Used in main and bowtiePrimer
def parseBowtie(guideClass, bowtieResultsFile, checkMismatch, scoreGC, scoreSelfComp,
                backbone, replace5prime, maxOffTargets, countMM, PAM, mode, scoringMethod=None,
                genome=None, gene=None, isoform=None, gene_isoforms=None):
    """ Parses bowtie hits and build list of guides"""

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
    if mode: # Cas9, Cpf1, Nickase and not TALEN
        sam[14] = sam[0].str[-(len(PAM) + 1):]
        sam[0] = sam[0].str[:-(len(PAM) + 1)]
        sam_name = sam.groupby(0).apply(lambda x, m=maxOffTargets: any(x.iloc[:, 14].value_counts() >= m))
        sam = sam.drop([14], axis=1)

        sam = sam.groupby([0, 1, 2, 3]).apply(# remove duplicates
            lambda x: x.sort_values(by=11).iloc[0])
        sam.rename(columns={0: "name", 11: "mm", 1: "str", 2: "chr", 3: "loc"}, inplace=True)
        sam = sam.sort_values(by=["name", "mm", "str", "chr", "loc"])
        sam = sam.reset_index(drop=True)

    for idx, row in sam.iterrows():
        line = list(row)
        if line[12] != line[12]:
            line = line[:-2]

        #  Encountered a new guide RNA (not a new hit for the same guide)
        elements = line[0].split(":") #removes from name 5' and 3' tails
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

    return guide_list

#Used in makePrimersFasta and makePrimersGenome
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

                s, e = int(position), int(position)+int(length)
                if m.group(1) == "RIGHT":
                    s, e = int(position)-int(length)+1, int(position)+1
                primerPos[label] = [s,e]

                primerFastaFile.write(">%s_%s_%s:%s_%s-%s\n%s\n" % (
                    target.ID, m.group(2), m.group(1), region, s, e, primers[(m.group(2), m.group(1), "SEQUENCE")]))

    return primers, primerPos

#Used in makePrimersFasta and makePrimersGenome
def get_primer_options(options):
    # Parse primer3 options. Update config if known option, otherwise append to primer3 input file
    primerOpt = ""

    if options:
        for opt in options.split(","):
            key, value = opt.split("=")
            if PRIMER3_CONFIG.has_key(key):
                PRIMER3_CONFIG[key] = value
            else:
                primerOpt += opt + "\n"

    return primerOpt

#Used in makePrimersFasta
def get_primer_query_sequence_fasta(target, outputDir, flank, fastaSequence):
    s = target.start-flank
    e = target.end+flank
    seqLenBeforeTarget = flank

    if s < 0:
        seqLenBeforeTarget -= abs(s)
        s = 0

    if e > len(fastaSequence):
        e = len(fastaSequence)

    return fastaSequence[s:e], seqLenBeforeTarget

#Used in makePrimerGnome
def get_primer_query_sequence_2bit(target, outputDir, flank, genome, twoBitToFaIndexDir, strand):
    s = target.start-flank
    seqLenBeforeTarget = flank

    if s < 0:
        seqLenBeforeTarget -= abs(s)
        s = 0

    prog = Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2>> %s/twoBitToFa.err" % (
        CONFIG["PATH"]["TWOBITTOFA"], target.chrom, s, target.end+flank, twoBitToFaIndexDir, genome, outputDir),
                 stdout=PIPE, shell=True)
    output = prog.communicate()

    if prog.returncode != 0:
        sys.stderr.write("Running twoBitToFa failed\n")
        sys.exit(EXIT['TWOBITTOFA_ERROR'])

    output = output[0].split("\n")
    del(output[0])
    seq = "".join(output)
    return seq, seqLenBeforeTarget

#Used in makePrimersFasta and makePrimersGenome
def runBowtiePrimers(primerFastaFileName, outputDir, genome, bowtieIndexDir, maxOffTargets):
    command = "%s -v 0 --best --sam-nohead -k 10 %s/%s -f %s -S %s/primer_results.sam 2> %s/bowtie_primers.err" % (
        CONFIG["PATH"]["BOWTIE"], bowtieIndexDir, genome, primerFastaFileName, outputDir, outputDir)
    prog = Popen(command, shell = True)
    prog.wait()

    if prog.returncode != 0:
        sys.stderr.write("Running bowtie on primers failed\n")
        sys.exit(EXIT['BOWTIE_PRIMER_ERROR'])

    return parseBowtie(Guide, "%s/primer_results.sam" % outputDir, False, False, False, None, None,
                       maxOffTargets, None, None, False, None, None)


#Used in makePrimersFasta and makePrimersGenome
def dump_restriction_sites(target, seq, flanks, enzymeCo, outputDir, minResSiteLen):
    sites = findRestrictionSites(seq, enzymeCo, minResSiteLen)
    out = [map(lambda x : [str(enzyme), x + target.start-flanks, enzyme.size], sites[enzyme]) for enzyme in sites]
    out = [item for sublist in out for item in sublist]
    out = sorted(out, key=itemgetter(1))

    # Assign tier to avoid overlaps
    siteCount = {}
    tiers = [0] * 23
    for site in out:
        tier = 0

        # count number of sites for each enzyme
        if not siteCount.has_key(site[0]):
            siteCount[site[0]] = 0
        siteCount[site[0]] += 1

        for j in range(len(tiers)):
            if site[1] > tiers[j]:
                tier = j
                tiers[j] = site[1]+site[2]
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

#Used in makePrimersFasta and makePrimersGenome
def dump_locus_sequence(target, outputDir, seq, seqLenBeforeTarget, strand):
    if strand == "-":
        seq = str(Seq(seq).complement())
    out = [[target.start-seqLenBeforeTarget, target.end, seq]]
    outputFile = open("%s/locusSeq_%s.json" % (outputDir, target.ID), 'w')
    json.dump(out, outputFile)
    outputFile.close()

#Used in makePrimersFasta and makePrimersGenome
def dump_genbank_file(seq, target, restSites, primers, outputDir, geneID, lociStart, strand):
    name= "%s, locus %s" % (geneID, target.ID)
    desc = "CHOPCHOP prediction for gene %s, target %s" % (geneID, target.ID)
    annotation = {"organism" : "Danio rerio", "Target location" : "chrT:1-20"}

    # Genbank file
    genbankFile = open('%s/%s_%s.gb' % (outputDir, geneID, target.ID), 'w')
    record = SeqRecord(Seq(seq, IUPACAmbiguousDNA()), description=desc, name="CHOPCHOP", id=name)
    record.annotation = annotation

    if target.strand == "+":
        ts = 1
    else:
        ts = -1

    record.features.append(SeqFeature(FeatureLocation(target.start-lociStart-1, target.end-lociStart-1, strand=ts),
                                      type="Target"))

    for primer in primers:
        record.features.append(SeqFeature(FeatureLocation(primers[primer][0], primers[primer][1]), type=primer))

    if strand == "-":
        record = record.reverse_complement()

    SeqIO.write(record, genbankFile, "genbank")
    genbankFile.close()

    pass

#Used in makePrimersFasta and makePrimersGenome
def pairPrimers(primerAttributes, primerList, outputDir):
    primers = {}

    for primer in primerList:
        guide, primerPairID, side = primer.ID.split("_")

        s = 0
        if side == "RIGHT": s = 1
        if not primers.has_key(guide): primers[guide] = {}
        if not primers[guide].has_key(primerPairID): primers[guide][primerPairID] = [None, None]
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
            output.append([ pair[0].chrom, pair[0].start, pair[0].end, pair[1].start, pair[1].end, i, pair[0].strand,
                            "%s" % lsq, "%s" % rsq, len(pair[0].offTargets), len(pair[1].offTargets),
                            len(offTargetPairs), ltm, rtm, size ])

            i += 1

        json.dump(output, outputFile)
        outputFile.close()

#Used in makePrimersFasta and makePrimersGenome
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
    primConfig['SEQUENCE_TARGET_START'] = str(seqLenBeforeTarget-padding)
    primConfig['SEQUENCE_TARGET_LEN'] = str(guide.targetSize+(2*padding))


    primer3InputFile = '%s/%s.primer3Input' % (outputDir, guide.ID)
    f = open(primer3InputFile, 'w')
    f.write(template.format(**primConfig))
    f.write(primer3options)
    f.write("=\n")
    f.close()

    command = "%s < %s 2>> %s/primer3.error" % (CONFIG["PATH"]["PRIMER3"], primer3InputFile, outputDir)
    # sys.stderr.write("%s\n" % command)
    prog = Popen(command, stdout = PIPE, shell=True)
    output = prog.communicate()

    if (prog.returncode != 0):
        sys.stderr.write("Running Primer3 failed\n");
        sys.exit(EXIT['PRIMER3_ERROR']);

    return output[0]


#Used in eval_CPF1_sequence and eval_CRISPR_sequence
def comaprePAM(basePAM, baseDNA):
    if basePAM == "N":
        return True

    if basePAM == baseDNA:
        return True

    if basePAM == "W" and (baseDNA == "A" or baseDNA == "T"):
        return True

    if basePAM == "S" and (baseDNA == "C" or baseDNA == "G"):
        return True

    if basePAM == "M" and (baseDNA == "A" or baseDNA == "C"):
        return True

    if basePAM == "K" and (baseDNA == "G" or baseDNA == "T"):
        return True

    if basePAM == "R" and (baseDNA == "A" or baseDNA == "G"):
        return True

    if basePAM == "Y" and (baseDNA == "C" or baseDNA == "T"):
        return True

    if basePAM == "B" and baseDNA != "A":
        return True

    if basePAM == "D" and baseDNA != "C":
        return True

    if basePAM == "H" and baseDNA != "G":
        return True

    if basePAM == "V" and baseDNA != "T":
        return True

    return False

#Used in eval_CPF1_sequence and eval_CRISPR_sequence
def permPAM(PAM):
    PAM = PAM.upper()
    new_comb = [""] # in case no PAM
    if len(PAM) == 1:
        new_comb = codes[PAM]

    for i in range(len(PAM) - 1):
        if i == 0:
            comb = codes[PAM[0]]
            new_comb = []
        else:
            comb = new_comb
            new_comb = []

        for j in codes[PAM[i + 1]]:
            for c in comb:
                new_comb.append(c + j)

    return new_comb

#####################
##
## JASON visualization
##


def complement(sequence):
    return sequence.translate(string.maketrans("ACGT", "TGCA"))


def FastaToViscoords(sequences, strand):
    """ Makes the exons in 'sequences' array generated in coordToFasta json readable for visualization"""
    exonstart = []
    exonend = []
    exonsequence = []

    for exon in sequences:
        # sys.stderr.write("%s\n" % exon)
        exonlist = exon.split(':')
        exoncoord = exonlist[2].split('-')
        exonstart.append(exoncoord[0])
        exonend.append(exoncoord[1])
        seq = sequences[exon]
        if strand == "-":
            seq = complement(seq)

        exonsequence.append(seq)

    return zip(exonstart, exonend, exonsequence)

#####################
##
## MAIN
##


#Used in ParseTargets
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

#Used in ParseTargets
def get_isoforms(gene, table_file):
    gene_isoforms = set()
    tableR = open(table_file, 'rb')
    tablereader = csv.DictReader(tableR, delimiter='\t', quoting=csv.QUOTE_NONE)
    for row in tablereader:
        if row['name2'] == gene:
            gene_isoforms.add(row['name'])
    tableR.close()
    return gene_isoforms


#Used in parseTargets
def filterRepeatingNames(txInfo, filter_names=["fix", "random", "alt"]):
    # if more isoforms have exact same name filter the ones
    # with "alt", "fix", "random" in chr names
    # then take the first one
    seen = []
    same_name_tx = []
    is_special = []
    for x in txInfo:
        if str(x[3]) not in seen:
            seen.append(str(x[3]))
            same_name_tx.append([x])
            is_special.append([any(fn in str(x[0]) for fn in filter_names)])
        else:
            idx = seen.index(str(x[3]))
            same_name_tx[idx].append(x)
            is_special[idx].append(any(fn in str(x[0]) for fn in filter_names))

    txInfo_ = []
    for i, tx in enumerate(same_name_tx):
        if any(is_special[i]) and sum(is_special[i]) < len(is_special[i]):
            idx = [i for i, x in enumerate(is_special[i]) if not x]
            txInfo_.append(tx[idx[0]])
        else:
            txInfo_.append(tx[0])

    return txInfo_


#Used in subsetExon
def hyphen_range(s):
    """ Takes a range in form of "a-b" and generate a list of numbers between a and b inclusive.
    Also accepts comma separated ranges like "a-b,c-d,f" will build a list which will include
    Numbers from a to b, a to d and f"""

    s = "".join(s.split()) #removes white space
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

#Used in parseTargets
def subsetExons(exons, targets):
    if exons:
        indices = hyphen_range(exons)
        for index in indices:
            if int(index) > len(targets):
                sys.stderr.write("That exon does not exist\n")
                sys.exit(EXIT['PYTHON_ERROR'])
        targets = [targets[int(i)-1] for i in indices] # indices is a list of exon numbers -1 e.g. exon 2 is [1]
    return targets


