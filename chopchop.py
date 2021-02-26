#!/usr/bin/env python3

import argparse
import resource
import pickle
#import featurization as feat
from Vars import *
from collections import defaultdict
from operator import itemgetter
from subprocess import Popen, PIPE
from functions.Main_Functions import *
from functions.Helper_Functions import *
from functions.TALEN_Specific_Functions import *

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

def getCoordinatesForJsonVisualization(args, visCoords, sequences, strand, resultCoords):
    # Coordinates for gene
    visCoordsFile = open('%s/viscoords.json' % args.outputDir, 'w')
    # visCoords = sorted(visCoords,  key=itemgetter(1))
    json.dump(visCoords, visCoordsFile)

    # Coordinates for sequence
    seqvis = FastaToViscoords(sequences, strand)
    seqvisFile = open('%s/seqviscoords.json' % args.outputDir, 'w')
    json.dump(seqvis, seqvisFile)

    # Coordinates for cutters
    cutCoord_file = open('%s/cutcoords.json' % args.outputDir, 'w')

    cutcoords = []
    for i in range(len(resultCoords)):
        el = []

        if args.MODE == ProgramMode.CRISPR or args.MODE == ProgramMode.CPF1:
            el.append(i + 1)
            el.extend(resultCoords[i])
        elif args.MODE == ProgramMode.TALENS or args.MODE == ProgramMode.NICKASE:
            el.extend(resultCoords[i])

        cutcoords.append(el)

    # Put bars at different heights to avoid overlap
    tiers = [0] * 23
    sortedCoords = sorted(cutcoords, key=itemgetter(1))
    for coord in sortedCoords:

        t = 0
        for j in range(len(tiers)):
            if coord[1] > tiers[j]:
                t = j
                tiers[j] = coord[1] + coord[3]
                break

        coord.append(t)

    json.dump(cutcoords, cutCoord_file)

    return cutcoords



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
