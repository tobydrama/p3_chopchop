#!/usr/bin/env python3

import argparse
import resource
import Guide
import Cas9
import MySQLdb
import re
import numpy
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from operator import attrgetter
from subprocess import Popen, PIPE
from Vars import *

soft, HARD_LIMIT = resource.getrlimit(resource.RLIMIT_NOFILE)
resource.setrlimit(resource.RLIMIT_NOFILE, (HARD_LIMIT, HARD_LIMIT))
def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-Target", "--targets", type=str, help="Target genes or regions", required=True)
    parser.add_argument("-r", "--gRVD", default="NH ", dest="g_RVD", action="store_const", const="NN ",
                        help="Use RVD 'NN' instead of 'NH' for guanine nucleotides. 'NH' appears to be more specific "
                             "than 'NN' but the choice depends on assembly kit.")
    parser.add_argument("-D", "--database",
                        help="Connect to a chopchop database to retrieve gene: user_name:passwd@host/database",
                        metavar="DATABASE", dest="database")
    parser.add_argument("-e", "--exon", help="Comma separated list of exon indices. Only find sites in this subset. ",
                        metavar="EXON_NUMBER", dest="exons")
    parser.add_argument("-TDP", "--targetDownstreamPromoter", default=200, type=int,
                        help="how many bp to target downstream of TSS")
    parser.add_argument("-TUP", "--targetUpstreamPromoter", default=200, type=int,
                        help="how many bp to target upstream of TSS")
    parser.add_argument("-G", "--genome", default="danRer7", metavar="GENOME",
                        help="The genome to search.")
    parser.add_argument("-g", "--guideSize", default=None, type=int, metavar="GUIDE_SIZE",
                        help="The size of the guide RNA.")
    parser.add_argument("-c", "--scoreGC", default=None, action="store_false",
                        help="Score GC content. True for CRISPR, False for TALENs.")
    parser.add_argument("-SC", "--noScoreSelfComp", default=None, action="store_false",
                        help="Do not penalize self-complementarity of CRISPR.")
    parser.add_argument("-BB", "--backbone", default=None, type=str,
                        help="Penalize self-complementarity versus backbone regions (comma-separated list, same strand "
                             "as guide). Requires -C.")
    parser.add_argument("-R5", "--replace5P", default=None, metavar="REPLACE_5P",
                        help="Replace bases from 5' end (with e.g. 'GG') ")  # FIX: AT THE MOMENT THIS IS ONLY APPLIES
    # TO FOLDING/SELF-COMPL
    parser.add_argument("-t", "--target", default="CODING", dest="targetRegion",
                        help="Target the whole gene CODING/WHOLE/UTR5/UTR3/SPLICE. Default is CODING.")
    parser.add_argument("-T", "--MODE", default=1, type=int, choices=[1, 2, 3, 4],
                        help="Set mode (int): default is Cas9 = 1, Talen = 2, Cpf1 = 3, Nickase = 4")
    parser.add_argument("-taleMin", "--taleMin", default=14, type=int,
                        help="Minimum distance between TALENs. Default is 14.")  # 14 + 18(length of TALE) = 32
    parser.add_argument("-taleMax", "--taleMax", default=20, type=int,
                        help="Maximum distance between TALENs. Default is 20.")  # 20 + 18(length of TALE) = 38
    parser.add_argument("-nickaseMin", "--nickaseMin", default=10, type=int,
                        help="Minimum distance between TALENs. Default is 10.")
    parser.add_argument("-nickaseMax", "--nickaseMax", default=31, type=int,
                        help="Maximum distance between TALENs. Default is 31.")
    parser.add_argument("-offtargetMaxDist", "--offtargetMaxDist", default=100, type=int,
                        help="Maximum distance between offtargets for Nickase. Default is 100.")
    parser.add_argument("-f", "--fivePrimeEnd", default="NN", type=str,
                        help="Specifies the requirement of the two nucleotides 5' end of the CRISPR guide: A/C/G/T/N. Default: NN.")
    parser.add_argument("-n", "--enzymeCo", default="N", metavar="ENZYME_CO",
                        help="The restriction enzyme company for TALEN spacer.")
    parser.add_argument("-R", "--minResSiteLen", type=int, default=4,
                        help="The minimum length of the restriction enzyme.")
    parser.add_argument("-v", "--maxMismatches", default=3, type=int, choices=[0, 1, 2, 3], metavar="MAX_MISMATCHES",
                        help="The number of mismatches to check across the sequence.")
    parser.add_argument("-m", "--maxOffTargets", metavar="MAX_HITS", help="The maximum number of off targets allowed.")
    parser.add_argument("-M", "--PAM", type=str, help="The PAM motif.")
    parser.add_argument("-o", "--outputDir", default="./", metavar="OUTPUT_DIR",
                        help="The output directory. Default is the current directory.")
    parser.add_argument("-F", "--fasta", default=False, action="store_true",
                        help="Use FASTA file as input rather than gene or genomic region.")
    parser.add_argument("-p", "--padSize", default=-1, type=int,
                        help="Extra bases searched outside the exon. Defaults to the size of the guide RNA for CRISPR and TALEN + maximum spacer for TALENS.")
    parser.add_argument("-P", "--makePrimers", default=False, action="store_true",
                        help="Designes primers using Primer3 to detect mutation.")
    parser.add_argument("-3", "--primer3options", default=None,
                        help="Options for Primer3. E.g. 'KEY1=VALUE1,KEY2=VALUE2'")
    parser.add_argument("-A", "--primerFlanks", default=300, type=int,
                        help="Size of flanking regions to search for primers.")
    parser.add_argument("-DF", "--displaySeqFlanks", default=300, type=int,
                        help="Size of flanking regions to output sequence into locusSeq_.")
    parser.add_argument("-a", "--guidePadding", default=20, type=int, help="Minimum distance of primer to target site.")
    parser.add_argument("-O", "--limitPrintResults", type=int, default=3000 if HARD_LIMIT > 3000 else HARD_LIMIT,
                        dest="limitPrintResults",
                        help="The number of results to print extended information for. Web server can handle 4k of these.")
    parser.add_argument("-w", "--uniqueMethod_Cong", default=False, dest="uniqueMethod_Cong", action="store_true",
                        help="A method to determine how unique the site is in the genome: allows 0 mismatches in last 15 bp.")
    parser.add_argument("-J", "--jsonVisualize", default=False, action="store_true",
                        help="Create files for visualization with json.")
    parser.add_argument("-nonO", "--nonOverlapping", default=False, action="store_true",
                        help="Will not produce overlapping guides, saves time, and recommended for permissive PAMs (e.g. Cas13d).")
    parser.add_argument("-scoringMethod", "--scoringMethod", default="G_20", type=str,
                        choices=["XU_2015", "DOENCH_2014", "DOENCH_2016", "MORENO_MATEOS_2015", "CHARI_2015", "G_20",
                                 "KIM_2018", "ALKAN_2018", "ZHANG_2019", "ALL"],
                        help="Scoring used for Cas9 and Nickase. Default is G_20. If a method fails to give scores, CHOPCHOP will output 0 instead of terminating.")
    parser.add_argument("-repairPredictions", "--repairPredictions", default=None, type=str,
                        choices=['mESC', 'U2OS', 'HEK293', 'HCT116', 'K562'],
                        help="Use inDelphi from Shen et al 2018 to predict repair profiles for every guideRNA, this will make .repProfile and .repStats files")
    parser.add_argument("-rm1perfOff", "--rm1perfOff", default=False, action="store_true",
                        help="For fasta input, don't score one off-target without mismatches.")
    parser.add_argument("-isoforms", "--isoforms", default=False, action="store_true",
                        help="Search for offtargets on the transcriptome.")
    parser.add_argument("-filterGCmin", "--filterGCmin", default=0, type=int,
                        help="Minimum required GC percentage. Default is 0.")
    parser.add_argument("-filterGCmax", "--filterGCmax", default=100, type=int,
                        help="Maximum allowed GC percentage. Default is 100.")
    parser.add_argument("-filterSelfCompMax", "--filterSelfCompMax", default=-1, type=int,
                        help="Maximum acceptable Self-complementarity score. Default is -1, no filter.")
    parser.add_argument("-consensusUnion", "--consensusUnion", default=False, action="store_true",
                        help="When calculating consensus sequence from multiple isoforms default uses intersection. This option specifies union of isoforms.")
    parser.add_argument("-BED", "--BED", default=False, action="store_true",
                        help="Create results as BED file, can be used for integration with UCSC.")
    parser.add_argument("-GenBank", "--GenBank", default=False, action="store_true",
                        help="Create results as GenBank file, sequence of targeted region with introns is included.")
    parser.add_argument("-offtargetsTable", "--offtargetsTable", default=False, action="store_true",
                        help="Create .tsv table with off-targets. Not all off-targets will be reported when early stopping will work on a guide! Limited also to CRISPR mode only and limited by --limitPrintResults option.")
    return parser.parse_args()


def mode_select(var, index, mode: ProgramMode):
    """ Selects a default depending on mode for options that have not been set """

    if var is not None:
        return var

    if mode == ProgramMode.CRISPR:
        return CRISPR_DEFAULT[index]

    elif mode == ProgramMode.TALENS:
        return TALEN_DEFAULT[index]

    elif mode == ProgramMode.CPF1:
        return CPF1_DEFAULT[index]

    elif mode == ProgramMode.NICKASE:
        return NICKASE_DEFAULT[index]

    sys.stderr.write("Unknown model %s\n" % mode)
    sys.exit(EXIT['PYTHON_ERROR'])


def getMismatchVectors(pam, gLength, cong):
    allowed = [True] * (gLength - len(pam))
    count = [True] * (gLength - len(pam))

    if cong:
        allowed = [True] * 9 + [False] * (gLength - len(pam) - 9)

    for char in pam:
        count.append(False)
        if char == "N":
            allowed.append(True)
        else:
            allowed.append(False)

    return allowed, count


def getCpf1MismatchVectors(pam, gLength):
    allowed = [True] * (gLength - len(pam))
    count = [True] * (gLength - len(pam))

    for char in pam[::-1]:
        count.insert(0, False)
        if char == "N":
            allowed.insert(0, True)
        else:
            allowed.insert(0, False)

    return allowed, count


def getAllowedFivePrime(allowed):
    new_allowed = []
    for el in allowed.split(","):
        if el[0] == 'N' and el[1] == 'N':
            return "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"
        elif el[0] == 'N':
            new_allowed.extend(["A" + el[1], "C" + el[1], "G" + el[1], "T" + el[1]])
        elif el[1] == 'N':
            new_allowed.extend([el[0] + "A", el[0] + "C", el[0] + "G", el[0] + "T"])
        else:
            new_allowed.append(el)
    return dict(zip(new_allowed, [True] * len(new_allowed)))


#####################
##
## CRISPR SPECIFIC FUNCTIONS
##
def gccontent(seq):
    gc = 0
    for i in seq:
        if i == 'G' or i == 'g' or i == 'C' or i == 'c':
            gc += 1
    return float(gc) / float(len(seq))


def selfComp(fwd, backbone):
    rvs = str(fwd.reverse_complement())
    fwd = str(fwd)
    L = len(fwd) - STEM_LEN - 1
    folding = 0
    for i in range(0, len(fwd) - STEM_LEN):
        if gccontent(fwd[i:i + STEM_LEN]) >= 0.5:
            if fwd[i:i + STEM_LEN] in rvs[0:(L - i)] or any(
                    [fwd[i:i + STEM_LEN] in item for item in backbone]):
                folding += 1

    return folding


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


def permPAM(PAM):
    PAM = PAM.upper()
    new_comb = [""]  # in case no PAM
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


def eval_CRISPR_sequence(name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, allowed, PAM,
                         filterGCmin, filterGCmax, filterSelfCompMax, replace5prime=None, backbone=None):
    """ Evaluates an k-mer as a potential CRISPR target site """

    gLen = guideSize - len(PAM)
    revCompPAM = str(Seq(PAM).reverse_complement())
    dna = Seq(dna)

    if str(dna[0:2]) in allowed:
        add = True
        for pos in range(len(PAM)):
            if comaprePAM(PAM[pos], dna[gLen + pos]):
                continue
            else:
                add = False
                break

        if add and (filterGCmin != 0 or filterGCmax != 100):
            gc = GC(
                dna[0:(None if PAM == "" else -len(PAM))])  # FIX EVERYWHERE GC content does not assumes 5' replacement
            if gc < filterGCmin or gc > filterGCmax:
                add = False

        if add and filterSelfCompMax != -1:
            if replace5prime:
                fwd = replace5prime + dna[len(replace5prime):(None if PAM == "" else -len(PAM))]
            else:
                fwd = dna[0:(None if PAM == "" else -len(PAM))]
            folding = selfComp(fwd, backbone)
            if folding > filterSelfCompMax:
                add = False

        # in order to control the number of mismatches to search in the last 8 or 3 bps,
        # need to reverse complement so the seed region can be at the start
        # rather than end of the sequence
        # not in isoforms case as we don't search reverse complement
        if add:
            if ISOFORMS:
                pam_comb = permPAM(PAM)
                for p in pam_comb:
                    fastaFile.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                        name, num, num + guideSize, downstream5prim, downstream3prim,
                        dna, p, dna[:gLen] + p))
                return True
            else:
                # all combinations of possible PAMs
                dna = dna.reverse_complement()
                pam_comb = permPAM(revCompPAM)
                for p in pam_comb:
                    fastaFile.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                        name, num, num + guideSize, downstream5prim, downstream3prim,
                        dna, p, p + dna[len(revCompPAM):]))
                return True

    if str(dna[-2:].reverse_complement()) in allowed and not ISOFORMS:
        add = True

        for pos in range(len(PAM)):
            if comaprePAM(revCompPAM[pos], dna[pos]):
                continue
            else:
                add = False
                break

        if add and (filterGCmin != 0 or filterGCmax != 100):
            gc = GC(dna[len(PAM):])
            if gc < filterGCmin or gc > filterGCmax:
                add = False

        if add and filterSelfCompMax != -1:
            if replace5prime:
                fwd = replace5prime + dna.reverse_complement()[len(PAM):-len(replace5prime)]
            else:
                fwd = dna.reverse_complement()[len(PAM):]
            folding = selfComp(fwd, backbone)
            if folding > filterSelfCompMax:
                add = False

        if add:
            pam_comb = permPAM(revCompPAM)
            for p in pam_comb:
                # on the reverse strand seq of 5' downstream becomes 3' downstream and vice versa
                fastaFile.write('>%s_%d-%d:%s:%s:-:%s:%s\n%s\n' % (
                    name, num, num + guideSize,
                    Seq(downstream3prim).reverse_complement(),
                    Seq(downstream5prim).reverse_complement(),
                    dna, p, p + dna[len(revCompPAM):]))
            return True

    return False


def sort_CRISPR_guides(guides):
    """ Sort pairs according to score  """
    return sorted(guides, key=attrgetter('score'))


#####################
##
## CPF1 SPECIFIC FUNCTIONS
##

def eval_CPF1_sequence(name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, PAM,
                       filterGCmin, filterGCmax, filterSelfCompMax, replace5prime=None, backbone=None):
    """ Evaluates an k-mer as a potential Cpf1 target site """

    gLen = guideSize - len(PAM)
    revCompPAM = str(Seq(PAM).reverse_complement())
    dna = Seq(dna)

    if replace5prime:
        fwd = dna[len(PAM):-len(replace5prime)] + replace5prime  # Replace the 2 first bases with e.g. "GG"
    else:
        fwd = dna[len(PAM):]  # Do not include PAM motif in folding calculations

    add = True
    for pos in range(len(PAM)):
        if comaprePAM(PAM[pos], dna[pos]):
            continue
        else:
            add = False
            break

    if add and (filterGCmin != 0 or filterGCmax != 100):
        gc = GC(dna[len(PAM):])
        if gc < filterGCmin or gc > filterGCmax:
            add = False

    if add and filterSelfCompMax != -1:
        if replace5prime:
            fwd = replace5prime + dna[len(PAM):-len(replace5prime)]
        else:
            fwd = dna[len(PAM):]
        folding = selfComp(fwd, backbone)
        if folding > filterSelfCompMax:
            add = False

    if add:
        if ISOFORMS:
            pam_comb = permPAM(PAM)
            for p in pam_comb:
                fastaFile.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                    name, num, num + guideSize, downstream5prim, downstream3prim,
                    dna, p, p + dna[len(PAM):]))
        else:
            dna = dna.reverse_complement()
            pam_comb = permPAM(revCompPAM)
            for p in pam_comb:
                fastaFile.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                    name, num, num + guideSize, downstream5prim, downstream3prim,
                    dna, p, dna[:gLen] + p))
        return True

    add = True and not ISOFORMS

    for pos in range(len(PAM)):
        if comaprePAM(revCompPAM[pos], dna[gLen + pos]):
            continue
        else:
            add = False
            break

    if add and (filterGCmin != 0 or filterGCmax != 100):
        gc = GC(dna.reverse_complement()[len(PAM):])
        if gc < filterGCmin or gc > filterGCmax:
            add = False

    if add and filterSelfCompMax != -1:
        if replace5prime:
            fwd = replace5prime + dna.reverse_complement()[len(PAM):-len(replace5prime)]
        else:
            fwd = dna.reverse_complement()[len(PAM):]
        folding = selfComp(fwd, backbone)
        if folding > filterSelfCompMax:
            add = False

    if add:
        pam_comb = permPAM(revCompPAM)
        for p in pam_comb:
            # on the reverse strand seq of 5' downstream becomes 3' downstream and vice versa
            fastaFile.write('>%s_%d-%d:%s:%s:-:%s:%s\n%s\n' % (
                name, num, num + guideSize,
                Seq(downstream3prim).reverse_complement(),
                Seq(downstream5prim).reverse_complement(),
                dna, p, dna[:gLen] + p))
        return True

    return False


def sort_TALEN_pairs(pairs):
    """ Sort pairs according to score and cluster """

    return sorted(pairs, key=attrgetter('score', 'cluster'))


def eval_TALENS_sequence(name, targetSize, dna, num, fastaFile, downstream5prim, downstream3prim):
    """ Evaluates an N-mer as a potential TALENs target site """
    del downstream5prim, downstream3prim
    found = False
    if dna[0] == "T":
        # dna = Seq(dna).reverse_complement()
        fastaFile.write('>%s_%d-%d\n%s\n' % (name, num, num + targetSize, dna))
        found = True
    elif dna[-1] == "A":
        fastaFile.write('>%s_%d-%d\n%s\n' % (name, num, num + targetSize, dna))
        found = True

    return found


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


def get_isoforms(gene, table_file):
    gene_isoforms = set()
    tableR = open(table_file, 'rb')
    tablereader = csv.DictReader(tableR, delimiter='\t', quoting=csv.QUOTE_NONE)
    for row in tablereader:
        if row['name2'] == gene:
            gene_isoforms.add(row['name'])
    tableR.close()
    return gene_isoforms


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


def subsetExons(exons, targets):
    if exons:
        indices = hyphen_range(exons)
        for index in indices:
            if int(index) > len(targets):
                sys.stderr.write("That exon does not exist\n")
                sys.exit(EXIT['PYTHON_ERROR'])
        targets = [targets[int(i) - 1] for i in indices]  # indices is a list of exon numbers -1 e.g. exon 2 is [1]
    return targets


def set_default_modes(args):
    if args.MODE == ProgramMode.CRISPR or ProgramMode.NICKASE:
        # Set mismatch checking policy
        (allowedMM, countMM) = getMismatchVectors(args.PAM, args.guideSize, args.uniqueMethod_Cong)
        allowed = getAllowedFivePrime(args.fivePrimeEnd)
        evalSequence = lambda name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim: eval_CRISPR_sequence(
            name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, allowed=allowed, PAM=args.PAM,
            filterGCmin=args.filterGCmin, filterGCmax=args.filterGCmax,
            filterSelfCompMax=args.filterSelfCompMax, replace5prime=args.replace5P, backbone=args.backbone)
        if args.MODE == ProgramMode.CRISPR:
            guideClass = Cas9 if not ISOFORMS else Guide
            sortOutput = sort_CRISPR_guides
        elif args.MODE == ProgramMode.NICKASE:
            guideClass = Cas9
            sortOutput = sort_TALEN_pairs

    elif args.MODE == ProgramMode.CPF1:
        (allowedMM, countMM) = getCpf1MismatchVectors(args.PAM, args.guideSize)
        evalSequence = lambda name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim: eval_CPF1_sequence(
            name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, PAM=args.PAM,
            filterGCmin=args.filterGCmin, filterGCmax=args.filterGCmax,
            filterSelfCompMax=args.filterSelfCompMax, replace5prime=args.replace5P, backbone=args.backbone)
        guideClass = ProgramMode.Cpf1 if not ISOFORMS else Guide
        sortOutput = sort_CRISPR_guides

    elif args.MODE == ProgramMode.TALENS:
        (allowedMM, countMM) = getMismatchVectors(args.PAM, args.guideSize, None)
        guideClass = Guide
        evalSequence = eval_TALENS_sequence
        sortOutput = sort_TALEN_pairs

    return evalSequence, guideClass, sortOutput


def connect_db(database_string):
    m = re.compile("(.+):(.+)@(.+)/(.+)").search(database_string)
    if not m:
        sys.stderr.write("Wrong syntax for connection string: username:pass@localhost/db_name")
        sys.exit(EXIT["DB_ERROR"])

    try:
        db = MySQLdb.connect(user=m.group(1), passwd=m.group(2), host=m.group(3), db=m.group(4))
    except:
        sys.stderr.write("Could not connect to database\n")
        sys.exit(EXIT['DB_ERROR'])

    return db


def geneToCoord_db(gene, organism, db):
    """ Gets genomic coordinates for a gene from a database """

    # Try refseq first
    lines = db.execute(
        "SELECT chrom, exonStarts, exonEnds, r.name, cdsStart, cdsEnd, strand, txStart, txEnd FROM organism o, refGene r WHERE o.assembly='%s' AND o.organism_id=r.organism_id AND (r.name='%s' OR r.name2='%s')" % (
            organism, gene, gene))

    # Then Ensembl
    if lines == 0:
        lines = db.execute(
            "SELECT chrom, exonStarts, exonEnds, r.name, cdsStart, cdsEnd, strand, txStart, txEnd FROM organism o, ensGene r LEFT OUTER JOIN ensemblToGeneName g ON r.name=g.name WHERE o.assembly='%s' AND o.organism_id=r.organism_id AND  (r.name='%s' OR r.name2='%s' OR g.value='%s')" % (
                organism, gene, gene, gene))

    # Then the general genePred table
    if lines == 0:
        lines = db.execute(
            "SELECT chrom, exonStarts, exonEnds, r.name, cdsStart, cdsEnd, strand, txStart, txEnd FROM organism o, gpGene r WHERE o.assembly='%s' AND o.organism_id=r.organism_id AND (r.name='%s' OR r.name2='%s')" % (
                organism, gene, gene))

    # Then wormbase. FIX: NO HARDCODED ASSEMBLY!!!
    if organism == "ce6" and lines == 0:
        lines = db.execute(
            "SELECT chrom, exonStarts, exonEnds, name, cdsStart, cdsEnd, strand, txStart, txEnd FROM sangerGene WHERE (name='%s' OR proteinID='%s')" % (
                gene, gene))

    # Then flybase. FIX: NO HARDCODED ASSEMBLY!!!
    if organism == "dm3" and lines == 0:
        lines = db.execute(
            "SELECT chrom, exonStarts, exonEnds, name, cdsStart, cdsEnd, strand, txStart, txEnd FROM flyBaseGene WHERE name='%s'" % (
                gene))

    if lines == 0:
        sys.stderr.write(
            "The gene name %s was not found in the gene sets for assembly %s. Consider trying an alternative ID (see the instruction page for supported gene identifiers) or using genomic coordinates. If you believe this type of ID should be supported for your organism contact us and we will do our best to support it. \n" % (
                gene, organism))
        sys.exit(EXIT['GENE_ERROR'])

    txInfo = []
    for i in range(lines):
        txInfo.append(db.fetchone())

    return txInfo


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


def parseTargets(target_string, genome, use_db, data, pad_size, target_region, exon_subset, ups_bp, down_bp,
                 index_dir, output_dir, use_union, make_vis, guideLen):
    targets = []
    vis_coords = []
    target_strand = "+"
    target_size = 0
    gene, isoform, gene_isoforms = (None, None, set())

    pattern = re.compile("(([.\w]+):)?([.,\d]+)-([.,\d]+)")
    is_coordinate = pattern.match(str(target_string))

    if is_coordinate:
        if ISOFORMS:
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
                sys.stderr.write("Unknown format: %s\n" % (target))
                sys.exit(EXIT['GENE_ERROR'])

    else:
        if use_db:
            if ISOFORMS:
                sys.stderr.write("--isoforms is not working with database search.\n")
                sys.exit(EXIT['ISOFORMS_ERROR'])
            txInfo = geneToCoord_db(target_string, genome, data)
            txInfo = filterRepeatingNames(txInfo)
        else:
            gene, txInfo = geneToCoord_file(target_string, data)
            txInfo = filterRepeatingNames(txInfo)
            isoform = "union" if use_union else "intersection"
            gene_isoforms = set([str(x[3]) for x in txInfo])
            if target_string in gene_isoforms:
                isoform = target_string
                gene_isoforms = get_isoforms(gene, data)

        target_chr = set([x[0] for x in txInfo])
        target_strand = set([x[6] for x in txInfo])
        isoforms = [str(x[3]) for x in txInfo]
        if len(target_strand) > 1 or len(target_chr) > 1:
            sys.stderr.write(
                "Specify which isoform you want to target as your query " + str(target_string) +
                " returns many isoforms: " + ', '.join(isoforms) +
                " which are from either inconsistent strands or chromosomes.\n")
            sys.exit(EXIT['GENE_ERROR'])
        else:
            target_strand = list(target_strand)[0]
            target_chr = list(target_chr)[0]

        for tx in txInfo:
            tx = list(tx)
            tx[4] = int(tx[4])
            tx[5] = int(tx[5])
            starts = tx[1].split(",")
            ends = tx[2].split(",")
            del starts[-1]
            del ends[-1]
            starts = map(int, starts)
            ends = map(int, ends)
            starts_v = starts[:]
            ends_v = ends[:]
            tx_vis = {"exons": [], "ATG": [], "name": tx[3]}

            if make_vis:
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
                prog = Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2> %s/twoBitToFa.err" % (
                    CONFIG["PATH"]["TWOBITTOFA"], tx[0], int(tx[4]) + 1, int(tx[5]) + 1, index_dir,
                    genome, output_dir), stdout=PIPE, shell=True)
                iso_seq = prog.communicate()
                if prog.returncode != 0:
                    sys.stderr.write("Running twoBitToFa when searching isoform sequence failed\n")
                    sys.exit(EXIT['TWOBITTOFA_ERROR'])

                iso_seq = iso_seq[0]
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

            # restrict isoforms
            coords = map(lambda x: [tx[0], x[0], x[1]], zip(starts, ends))
            if tx[6] == "-":
                coords.reverse()
            coords = subsetExons(exon_subset, coords)
            if tx[6] == "-":
                coords.reverse()

            # Truncate to region
            if target_region == "CODING":
                coords = truncateToCoding(tx[4], tx[5], coords)
            elif target_region == "UTR5":
                if tx[6] == "+":
                    coords = truncateToUTR5(tx[4], coords)
                else:
                    coords = truncateToUTR3(tx[5], coords)
            elif target_region == "PROMOTER":
                coords = truncateToPROMOTER(tx[6], coords, ups_bp, down_bp)
            elif target_region == "UTR3":
                if tx[6] == "+":
                    coords = truncateToUTR3(tx[5], coords)
                else:
                    coords = truncateToUTR5(tx[4], coords)
            elif target_region == "SPLICE":
                coords = truncateToSplice(coords)
            elif target_region != "WHOLE":
                sys.stderr.write("Unknown region: %s\n" % target_region)
                sys.exit(EXIT['PYTHON_ERROR'])

            # filter exons that are too truncated
            coords = [x for x in coords if x[1] < x[2]]
            if not coords:
                if gene_isoforms:
                    gene_isoforms.remove(tx[3])
                if vis_coords:
                    del vis_coords[-1]

            # compute intersection/union on all exons
            if txInfo[0][3] == tx[3]:  # if this is first of the isoforms
                for x in coords:
                    targets.extend(range(x[1], x[2] + 1))
                targets = set(targets)
            else:
                if not use_union:
                    targets_ = []
                    for x in coords:
                        targets_.extend(range(x[1], x[2] + 1))

                    if len(targets_) >= guideLen:  # cover cases where some transcripts provide short or none bp
                        targets &= set(targets_)

                    if len(targets) < guideLen:
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

        target_size = len(targets)
        if target_size < guideLen:
            sys.stderr.write("Search region is too small. You probably want to specify -t option as WHOLE")
            sys.exit(EXIT['GENE_ERROR'])

        starts, ends = bins(targets)
        if ISOFORMS:
            targets = map(lambda x: "%s:%s-%s" % (target_chr, x[0], x[1]), zip(starts, ends))
        else:
            targets = map(lambda x: "%s:%s-%s" % (target_chr, x[0] - pad_size, x[1] + pad_size), zip(starts, ends))

    if target_size > TARGET_MAX:
        sys.stderr.write("Search region is too large (%s nt). Maximum search region is %s nt.\n" % (
            target_size, TARGET_MAX))
        sys.exit(EXIT['GENE_ERROR'])

    return targets, vis_coords, target_strand, gene, isoform, gene_isoforms


def parseFastaTarget(fasta_file, candidate_fasta_file, target_size, eval_and_print):
    """ Parse a FASTA file as input for targeting """

    fasta_file = list(SeqIO.parse(fasta_file, 'fasta'))
    seq_name, sequence = fasta_file[0].id, str(fasta_file[0].seq)

    name = "%s:0-%s" % (seq_name, len(sequence))
    id_name = "C:" + name
    sequence = sequence.upper()
    sequence = "".join(sequence.split())

    dna_pattern = re.compile(r'([^ACGTNacgtn])')
    if dna_pattern.search(sequence):
        sys.stderr.write("Input sequence contains illegal characters.\n")
        sys.exit(EXIT['GENE_ERROR'])

    sequences = {}
    candidate_fasta_file = open(candidate_fasta_file, 'w')

    # Loop over sequence, write every k-mer into file in which k-mer ends in as PAM in fasta format
    for num in range(0, len(sequence) - (target_size - 1)):

        if (num - DOWNSTREAM_NUC) > 0:
            start5prim = num - DOWNSTREAM_NUC
        else:
            start5prim = 0

        if (num + target_size + DOWNSTREAM_NUC) > len(sequence):
            end3prim = len(sequence)
        else:
            end3prim = num + target_size + DOWNSTREAM_NUC

        downstream_5prim = sequence[start5prim:num]
        downstream_3prim = sequence[(num + target_size):end3prim]

        if eval_and_print(id_name, target_size, sequence[num:(num + target_size)], num,
                          candidate_fasta_file, downstream_5prim, downstream_3prim):
            sequences[id_name] = sequence

    return sequences, [name], [{"exons": [[seq_name, 1, len(sequence), 0, 20, "+"]],
                                "ATG": [], "name": seq_name}], sequence, "+"


#####################
##
## Functions
##

def coordToFasta(regions, fasta_file, outputDir, targetSize, evalAndPrintFunc, nonOver, indexDir, genome, strand, ext):
    """ Extracts the sequence corresponding to genomic coordinates from a FASTA file """

    ext = 0 if ISOFORMS else ext  # for genomic context for some models
    sequences = {}
    fasta_file = open(fasta_file, 'w')
    fasta_seq = ""

    if ISOFORMS and strand == "-":
        regions = regions[::-1]

    for region in regions:
        # Extracts chromosome number and region start and end
        chrom = region[0:region.rfind(':')]
        start = int(region[region.rfind(':') + 1:region.rfind('-')])
        finish = int(region[region.rfind('-') + 1:])
        start = max(start, 0)

        if ext == 0 and finish == start:
            continue

        # Run twoBitToFa program to get actual dna sequence corresponding to input genomic coordinates
        # Popen runs twoBitToFa program. PIPE pipes stdout.
        prog = Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2> %s/twoBitToFa.err" % (
            CONFIG["PATH"]["TWOBITTOFA"], chrom, start - ext, finish + ext, indexDir, genome, outputDir), stdout=PIPE,
                     shell=True)

        # Communicate converts stdout to a string
        output = prog.communicate()
        if prog.returncode != 0:
            sys.stderr.write("Running twoBitToFa failed\n")
            sys.exit(EXIT['TWOBITTOFA_ERROR'])

        output = output[0]
        exons = output.split("\n")
        dna = ''.join(exons[1:]).upper()
        ext_dna = dna
        dna = dna[ext:(len(dna) - ext)]
        if len(dna) != (finish - start):  # something is wrong with what was fetched by twoBitToFa
            continue

        if ISOFORMS and strand == "-":
            dna = str(Seq(dna).reverse_complement())

        # Write exon sequences to text file user can open in ApE. exon-intron junctions in lowercase.
        fasta_seq += dna[0].lower() + dna[1:-1] + dna[-1].lower()

        # Add 1 due to BED 0-indexing
        name = "C:%s:%d-%d" % (chrom, start, finish)

        # Loop over exon sequence, write every g-mer into file in which g-mer ends in PAM in fasta format
        positions = range(0, len(dna) - (targetSize - 1))
        while len(positions) != 0:
            num = positions.pop(0)
            downstream_5prim = ext_dna[num:(num + ext)]
            g_end = num + ext + targetSize
            downstream_3prim = ext_dna[g_end:(g_end + ext)]
            if evalAndPrintFunc(name, targetSize, dna[num:(num + targetSize)],
                                len(dna) - num - targetSize if ISOFORMS and strand == "-" else num, fasta_file,
                                downstream_5prim, downstream_3prim):
                if nonOver:  # positions overlapping those of this guide
                    for p in range(num, num + targetSize):
                        if p in positions:
                            positions.remove(p)

                if name not in sequences:
                    sequences[name] = dna

    fasta_file.close()

    if ISOFORMS and strand == "-":
        fasta_seq = str(Seq(fasta_seq).reverse_complement())

    return sequences, fasta_seq


def truncateToUTR5(cds_start, exons):
    """ Truncates the gene to only target 5' UTR """

    end_exon = 0
    for exon in range(len(exons)):
        if (cds_start > exons[exon][1]) and (cds_start < exons[exon][2]):
            exons[exon][2] = cds_start
            end_exon = exon
            break

    return exons[:end_exon + 1]


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


def truncateToUTR3(cds_end, exons):
    """ Truncates the gene to only target 3' UTR """

    start_exon = 0
    for exon in range(len(exons)):
        if (cds_end > exons[exon][1]) and (cds_end < exons[exon][2]):
            exons[exon][1] = cds_end
            start_exon = exon

    return exons[start_exon:]


def truncateToSplice(exons):
    """ Truncates the gene to only target splice sites """

    splice_sites = []
    for ind in range(0, len(exons)):
        splice_sites.append([exons[ind][0], exons[ind][1] - 1, exons[ind][1] + 1])
        splice_sites.append([exons[ind][0], exons[ind][2] - 1, exons[ind][2] + 1])
    # Remove first and last (i.e. transcription start and termination site)
    return splice_sites[1:len(splice_sites) - 1]


def truncateToCoding(cds_start, cds_end, exons):
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


def main():
    # Parse arguments
    args = parseArguments()

    # set isoforms to global as it is influencing many steps
    global ISOFORMS
    ISOFORMS = args.isoforms

    # Add TALEN length
    args.taleMin += 18
    args.taleMax += 18

    # Set mode specific parameters if not set by user
    args.scoreGC = mode_select(args.scoreGC, "SCORE_GC", args.MODE)
    args.scoreSelfComp = mode_select(args.noScoreSelfComp, "SCORE_FOLDING", args.MODE)
    args.PAM = mode_select(args.PAM, "PAM", args.MODE)
    args.guideSize = mode_select(args.guideSize, "GUIDE_SIZE", args.MODE) + len(args.PAM)
    args.maxMismatches = mode_select(args.maxMismatches, "MAX_MISMATCHES", args.MODE)
    args.maxOffTargets = mode_select(args.maxOffTargets, "MAX_OFFTARGETS", args.MODE)

    # Add TALEN length
    args.nickaseMin += args.guideSize
    args.nickaseMax += args.guideSize

    if args.scoreSelfComp:
        if args.backbone:
            tmp = args.backbone.strip().split(",")
            args.backbone = [str(Seq(el).reverse_complement()) for el in tmp]
        else:
            args.backbone = []

    # Pad each exon equal to guidesize unless
    if args.padSize != -1:
        padSize = args.padSize
    else:
        if args.MODE == ProgramMode.TALENS:
            padSize = args.taleMax
        elif args.MODE == ProgramMode.NICKASE:
            padSize = args.nickaseMax
        elif args.MODE == ProgramMode.CRISPR or args.MODE == ProgramMode.CPF1:
            padSize = args.guideSize

    # Set default functions for different modes
    evalSequence, guideClass, sortOutput = set_default_modes(args)

    # Connect to database if requested
    if args.database:
        cdb = connect_db(args.database)
        db = cdb.cursor()
        use_db = True
    else:
        db = "%s/%s.gene_table" % (
            CONFIG["PATH"]["GENE_TABLE_INDEX_DIR"] if not ISOFORMS else CONFIG["PATH"]["ISOFORMS_INDEX_DIR"],
            args.genome)
        use_db = False

    # Create output directory if it doesn't exist
    if not os.path.isdir(args.outputDir):
        os.mkdir(args.outputDir)

    candidate_fasta_file = '%s/sequence.fa' % args.outputDir
    gene, isoform, gene_isoforms = (None, None, set())
    if args.fasta:
        sequences, targets, visCoords, fastaSequence, strand = parseFastaTarget(
            args.targets, candidate_fasta_file, args.guideSize, evalSequence)
    else:
        targets, visCoords, strand, gene, isoform, gene_isoforms = parseTargets(
            args.targets, args.genome, use_db, db, padSize, args.targetRegion, args.exons,
            args.targetUpstreamPromoter, args.targetDownstreamPromoter,
            CONFIG["PATH"]["TWOBIT_INDEX_DIR"] if not ISOFORMS else CONFIG["PATH"]["ISOFORMS_INDEX_DIR"],
            args.outputDir, args.consensusUnion, args.jsonVisualize, args.guideSize)
        sequences, fastaSequence = coordToFasta(
            targets, candidate_fasta_file, args.outputDir, args.guideSize, evalSequence, args.nonOverlapping,
            CONFIG["PATH"]["TWOBIT_INDEX_DIR"] if not ISOFORMS else CONFIG["PATH"]["ISOFORMS_INDEX_DIR"],
            args.genome, strand, DOWNSTREAM_NUC)

    # Converts genomic coordinates to fasta file of all possible k-mers
    if len(sequences) == 0:
        sys.stderr.write("No target sites\n")
        sys.exit()


if __name__ == '__main__':
    main()
