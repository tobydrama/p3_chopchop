#!/usr/bin/env python3
import logging
import subprocess
from operator import itemgetter
from typing import List, Callable, Union

import argparse
import resource

import scoring
from Vars import *

from classes.Guide import Guide
from classes.PAIR import Pair
from classes.ProgramMode import ProgramMode
import functions.Main_Functions as mainFunctions
import functions.Helper_Functions as helperFunctions
import functions.TALEN_Specific_Functions as talenFunctions


soft, HARD_LIMIT = resource.getrlimit(resource.RLIMIT_NOFILE)
resource.setrlimit(resource.RLIMIT_NOFILE, (HARD_LIMIT, HARD_LIMIT))


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-Target", "--targets", type=str, required=True,
                        help="Target genes or regions")

    parser.add_argument("-r", "--gRVD", default="NH ", const="NN ", action="store_const", dest="g_RVD",
                        help="Use RVD 'NN' instead of 'NH' for guanine nucleotides. 'NH' appears to be more specific "
                             "than 'NN' but the choice depends on assembly kit.")

    parser.add_argument("-D", "--database", metavar="DATABASE", dest="database",
                        help="Connect to a chopchop database to retrieve gene: user_name:passwd@host/database")

    parser.add_argument("-e", "--exon", metavar="EXON_NUMBER", dest="exons",
                        help="Comma separated list of exon indices. Only find sites in this subset. ")

    parser.add_argument("-TDP", "--targetDownstreamPromoter", type=int, default=200,
                        help="how many bp to target downstream of TSS")

    parser.add_argument("-TUP", "--targetUpstreamPromoter", type=int, default=200,
                        help="how many bp to target upstream of TSS")

    parser.add_argument("-G", "--genome", default="danRer7", metavar="GENOME",
                        help="The genome to search.")

    parser.add_argument("-g", "--guideSize", type=int, default=None, metavar="GUIDE_SIZE",
                        help="The size of the guide RNA.")

    parser.add_argument("-c", "--scoreGC", default=None, action="store_false",
                        help="Score GC content. True for CRISPR, False for TALENs.")

    parser.add_argument("-SC", "--noScoreSelfComp", default=None, action="store_false",
                        help="Do not penalize self-complementarity of CRISPR.")

    parser.add_argument("-BB", "--backbone", type=str, default=None,
                        help="Penalize self-complementarity versus backbone regions (comma-separated list, same strand "
                             "as guide). Requires -C.")

    # TODO FIX: AT THE MOMENT THIS IS ONLY APPLIES TO FOLDING/SELF-COMPL
    parser.add_argument("-R5", "--replace5P", default=None, metavar="REPLACE_5P",
                        help="Replace bases from 5' end (with e.g. 'GG') ")

    parser.add_argument("-t", "--target", default="CODING", dest="targetRegion",
                        help="Target the whole gene CODING/WHOLE/UTR5/UTR3/SPLICE. Default is CODING.")

    parser.add_argument("-T", "--MODE", type=int, default=1, choices=[1, 2, 3, 4],
                        help="Set mode (int): default is Cas9 = 1, Talen = 2, Cpf1 = 3, Nickase = 4")

    # 14 + 18(length of TALE) = 32
    parser.add_argument("-taleMin", "--taleMin", type=int, default=14,
                        help="Minimum distance between TALENs. Default is 14.")

    # 20 + 18(length of TALE) = 38
    parser.add_argument("-taleMax", "--taleMax", type=int, default=20,
                        help="Maximum distance between TALENs. Default is 20.")

    parser.add_argument("-nickaseMin", "--nickaseMin", type=int, default=10,
                        help="Minimum distance between TALENs. Default is 10.")

    parser.add_argument("-nickaseMax", "--nickaseMax",  type=int, default=31,
                        help="Maximum distance between TALENs. Default is 31.")

    parser.add_argument("-offtargetMaxDist", "--offtargetMaxDist", type=int, default=100,
                        help="Maximum distance between offtargets for Nickase. Default is 100.")

    parser.add_argument("-f", "--fivePrimeEnd", type=str, default="NN",
                        help="Specifies the requirement of the two nucleotides 5' end of the CRISPR guide: A/C/G/T/N. "
                             "Default: NN.")

    parser.add_argument("-n", "--enzymeCo", default="N", metavar="ENZYME_CO",
                        help="The restriction enzyme company for TALEN spacer.")

    parser.add_argument("-R", "--minResSiteLen", type=int, default=4,
                        help="The minimum length of the restriction enzyme.")

    parser.add_argument("-v", "--maxMismatches", type=int, choices=[0, 1, 2, 3], default=3, metavar="MAX_MISMATCHES",
                        help="The number of mismatches to check across the sequence.")

    parser.add_argument("-m", "--maxOffTargets", metavar="MAX_HITS",
                        help="The maximum number of off targets allowed.")

    parser.add_argument("-M", "--PAM", type=str,
                        help="The PAM motif.")

    parser.add_argument("-o", "--outputDir", default="./", metavar="OUTPUT_DIR",
                        help="The output directory. Default is the current directory.")

    parser.add_argument("-F", "--fasta", default=False, action="store_true",
                        help="Use FASTA file as input rather than gene or genomic region.")

    parser.add_argument("-p", "--padSize", type=int, default=-1,
                        help="Extra bases searched outside the exon. Defaults to the size of the guide RNA for CRISPR "
                             "and TALEN + maximum spacer for TALENS.")

    parser.add_argument("-P", "--makePrimers", default=False, action="store_true",
                        help="Designds primers using Primer3 to detect mutation.")

    parser.add_argument("-3", "--primer3options", default=None,
                        help="Options for Primer3. E.g. 'KEY1=VALUE1,KEY2=VALUE2'")

    parser.add_argument("-A", "--primerFlanks", type=int, default=300,
                        help="Size of flanking regions to search for primers.")

    parser.add_argument("-DF", "--displaySeqFlanks", type=int, default=300,
                        help="Size of flanking regions to output sequence into locusSeq_.")

    parser.add_argument("-a", "--guidePadding", type=int, default=20,
                        help="Minimum distance of primer to target site.")

    parser.add_argument("-O", "--limitPrintResults", type=int, default=(3000 if HARD_LIMIT > 3000 else HARD_LIMIT),
                        dest="limitPrintResults",
                        help="The number of results to print extended information for. Web server can handle 4k of "
                             "these.")

    parser.add_argument("-w", "--uniqueMethod_Cong", default=False, action="store_true", dest="uniqueMethod_Cong",
                        help="A method to determine how unique the site is in the genome: allows 0 mismatches in last "
                             "15 bp.")

    parser.add_argument("-J", "--jsonVisualize", default=False, action="store_true",
                        help="Create files for visualization with json.")

    parser.add_argument("-nonO", "--nonOverlapping", default=False, action="store_true",
                        help="Will not produce overlapping guides, saves time, and recommended for permissive PAMs ("
                             "e.g. Cas13d).")

    parser.add_argument("-scoringMethod", "--scoringMethod", type=str, default="G_20",
                        choices=["XU_2015", "DOENCH_2014", "DOENCH_2016", "MORENO_MATEOS_2015", "CHARI_2015", "G_20",
                                 "KIM_2018", "ALKAN_2018", "ZHANG_2019", "ALL"],
                        help="Scoring used for Cas9 and Nickase. Default is G_20. If a method fails to give scores, "
                             "CHOPCHOP will output 0 instead of terminating.")

    parser.add_argument("-repairPredictions", "--repairPredictions", type=str, default=None,
                        choices=['mESC', 'U2OS', 'HEK293', 'HCT116', 'K562'],
                        help="Use inDelphi from Shen et al 2018 to predict repair profiles for every guideRNA, "
                             "this will make .repProfile and .repStats files")

    parser.add_argument("-rm1perfOff", "--rm1perfOff", default=False, action="store_true",
                        help="For fasta input, don't score one off-target without mismatches.")

    parser.add_argument("-isoforms", "--isoforms", default=False, action="store_true",
                        help="Search for offtargets on the transcriptome.")

    parser.add_argument("-filterGCmin", "--filterGCmin", type=int, default=0,
                        help="Minimum required GC percentage. Default is 0.")

    parser.add_argument("-filterGCmax", "--filterGCmax", type=int, default=100,
                        help="Maximum allowed GC percentage. Default is 100.")

    parser.add_argument("-filterSelfCompMax", "--filterSelfCompMax", type=int, default=-1,
                        help="Maximum acceptable Self-complementarity score. Default is -1, no filter.")

    parser.add_argument("-consensusUnion", "--consensusUnion", default=False, action="store_true",
                        help="When calculating consensus sequence from multiple isoforms default uses intersection. "
                             "This option specifies union of isoforms.")

    parser.add_argument("-BED", "--BED", default=False, action="store_true",
                        help="Create results as BED file, can be used for integration with UCSC.")

    parser.add_argument("-GenBank", "--GenBank", default=False, action="store_true",
                        help="Create results as GenBank file, sequence of targeted region with introns is included.")

    parser.add_argument("-offtargetsTable", "--offtargetsTable", default=False, action="store_true",
                        help="Create .tsv table with off-targets. Not all off-targets will be reported when early "
                             "stopping will work on a guide! Limited also to CRISPR mode only and limited by "
                             "--limitPrintResults option.")

    parser.add_argument("--logLevel", type=str, default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], dest="log_level",
                        help="Set logging level.")

    args = parser.parse_args()

    # Logging
    logging.basicConfig(level=logging.getLevelName(args.log_level.upper()),
                        format="%(levelname)-8s:: %(message)s")
    logging.debug("Log level set to %s." % args.log_level)

    # Change args.MODE type from int to ProgramMode
    args.MODE = ProgramMode(args.MODE)

    # Add TALEN length
    args.taleMin += 18
    args.taleMax += 18

    # Set mode specific parameters if not set by user
    args.scoreGC = mainFunctions.mode_select(args.scoreGC, "SCORE_GC", args.MODE)
    args.scoreSelfComp = mainFunctions.mode_select(args.noScoreSelfComp, "SCORE_FOLDING", args.MODE)
    args.PAM = mainFunctions.mode_select(args.PAM, "PAM", args.MODE)
    args.guideSize = mainFunctions.mode_select(args.guideSize, "GUIDE_SIZE", args.MODE) + len(args.PAM)
    args.maxMismatches = mainFunctions.mode_select(args.maxMismatches, "MAX_MISMATCHES", args.MODE)
    args.maxOffTargets = mainFunctions.mode_select(args.maxOffTargets, "MAX_OFFTARGETS", args.MODE)

    # Add TALEN length
    args.nickaseMin += args.guideSize
    args.nickaseMax += args.guideSize

    if args.scoreSelfComp:
        if args.backbone:
            tmp = args.backbone.strip().split(",")
            args.backbone = [str(helperFunctions.Seq(el).reverse_complement()) for el in tmp]
        else:
            args.backbone = []

    logging.debug("Finished parsing arguments.")

    return args


def getCoordinatesForJsonVisualization(args, visCoords, sequences, strand, resultCoords):
    # Coordinates for gene
    visCoordsFile = open('%s/viscoords.json' % args.outputDir, 'w')
    # visCoords = sorted(visCoords,  key=itemgetter(1))
    helperFunctions.json.dump(visCoords, visCoordsFile)

    # Coordinates for sequence
    seqvis = mainFunctions.FastaToViscoords(sequences, strand)
    seqvisFile = open('%s/seqviscoords.json' % args.outputDir, 'w')
    helperFunctions.json.dump(seqvis, seqvisFile)

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

    helperFunctions.json.dump(cutcoords, cutCoord_file)

    return cutcoords


def getClusterPairsTALENS(results, sequences, args):
    pairs = talenFunctions.pairTalens(results, sequences, args.guideSize, int(args.taleMin), int(args.taleMax), args.enzymeCo,
                                      args.maxOffTargets, args.g_RVD, args.minResSiteLen)

    if (not len(pairs)):
        sys.stderr.write("No TALEN pairs could be generated for this region.\n")
        sys.exit(helperFunctions.EXIT['GENE_ERROR'])

    if args.rm1perfOff and args.fasta:
        for pair in pairs:
            if pair.diffStrandOffTarget > 0:
                pair.score = pair.score - SCORE["OFFTARGET_PAIR_DIFF_STRAND"]
            if pair.sameStrandOffTarget > 0:
                pair.score = pair.score - SCORE["OFFTARGET_PAIR_SAME_STRAND"]

    cluster, results = talenFunctions.clusterPairs(pairs)
    return cluster, results


def getClusterPairsNICKASE(results, sequences, args):
    pairs = talenFunctions.pairCas9(results, sequences, args.guideSize, int(args.nickaseMin), int(args.nickaseMax), args.enzymeCo,
                                    args.maxOffTargets, args.minResSiteLen, args.offtargetMaxDist)

    if (not len(pairs)):
        sys.stderr.write("No Cas9 nickase pairs could be generated for this region.\n")
        sys.exit(helperFunctions.EXIT['GENE_ERROR'])

    if args.rm1perfOff and args.fasta:
        for pair in pairs:
            if pair.diffStrandOffTarget > 0:
                pair.score = pair.score - SCORE["OFFTARGET_PAIR_DIFF_STRAND"]

    cluster, results = talenFunctions.clusterPairs(pairs)
    return cluster, results


def print_scores(sorted_output: Union[List[Guide], List[Pair]],
                 mode: ProgramMode = ProgramMode.CRISPR,
                 scoring_method: str = "G20",
                 isoforms: bool = False) -> None:
    if isoforms:
        print("Rank\tTarget sequence\tGenomic location\tGene\tIsoform\tGC content (%)\tSelf-complementarity\t"
              "Local structure\tMM0\tMM1\tMM2\tMM3\tConstitutive\tIsoformsMM0\tIsoformsMM1\tIsoformsMM2\tIsoformsMM3")
        for i, guide in enumerate(sorted_output):
            print("%s\t%s" % (i + 1, guide))

    else:
        if mode == ProgramMode.CRISPR:
            common_header = "Rank\tTarget sequence\tGenomic location\tStrand\tGC content (%)\tSelf-complementarity\t" \
                            "MM0\tMM1\tMM2\tMM3 "

            if scoring_method == "ALL":
                print(common_header + "\tXU_2015\tDOENCH_2014\tDOENCH_2016\tMORENO_MATEOS_2015\tCHARI_2015\tG_20\t"
                                      "ALKAN_2018\tZHANG_2019")
            else:
                print(common_header + "\tEfficiency")

            for i, guide in enumerate(sorted_output):
                print("%s\t%s" % (i + 1, guide))

        elif mode == ProgramMode.CPF1:
            print("Rank\tTarget sequence\tGenomic location\tStrand\tGC content (%)\tSelf-complementarity\tEfficiency\t"
                  "MM0\tMM1\tMM2\tMM3")

            for i, guide in enumerate(sorted_output):
                print("%s\t%s" % (1 + i, guide))

        elif mode == ProgramMode.TALENS or mode == ProgramMode.NICKASE:
            if mode == ProgramMode.TALENS:
                print("Rank\tTarget sequence\tGenomic location\tTALE 1\tTALE 2\tCluster\tOff-target pairs\t"
                      "Off-targets MM0\tOff-targets MM1\tOff-targets MM2\tOff-targets MM3\tRestriction sites\tBest ID")
            else:
                print("Rank\tTarget sequence\tGenomic location\tCluster\tOff-target pairs\t"
                      "Off-targets MM0\tOff-targets MM1\tOff-targets MM2\tOff-targets MM3\tRestriction sites\tBest ID")

            for i, guide in enumerate(sorted_output):
                print("%s\t%s\t%s" % (i + 1, guide, guide.ID))


def generate_result_coordinates(sorted_output: Union[List[Guide], List[Pair]],
                                clusters: List[List[Pair]],
                                sort_function: Callable[[Union[List[Guide], List[Pair]]], List[Guide]],
                                mode: ProgramMode = ProgramMode.CRISPR,
                                isoforms: bool = False) -> List[List[any]]:
    result_coordinates = []
    if isoforms or mode in [ProgramMode.CRISPR, ProgramMode.CPF1]:
        for i, guide in enumerate(sorted_output):
            result_coordinates.append([guide.start, guide.score, guide.guideSize, guide.strand])

    elif mode in [ProgramMode.TALENS, ProgramMode.NICKASE]:
        final_output = []
        for cluster in clusters:
            if len(cluster) > 0:
                final_output.append(cluster[0])

        sorted_output = sort_function(final_output)
        for guide in sorted_output:
            result_coordinates.append([guide.spacerStart, guide.score, guide.spacerSize, guide.strand, guide.ID,
                                       guide.tale1.start, guide.tale2.end])

    return result_coordinates


def main():
    # Parse arguments
    args = parse_arguments()

    # set isoforms to global as it is influencing many steps
    global ISOFORMS
    ISOFORMS = args.isoforms

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
    # new function
    countMM, evalSequence, guideClass, sortOutput = mainFunctions.set_default_modes(args)

    ### General ARGPARSE done, upcoming Target parsing

    # Connect to database if requested
    if args.database:
        cdb = mainFunctions.connect_db(args.database)
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
        sequences, targets, visCoords, fastaSequence, strand = mainFunctions.parseFastaTarget(
            args.targets, candidate_fasta_file, args.guideSize, evalSequence)
    else:
        targets, visCoords, strand, gene, isoform, gene_isoforms = mainFunctions.parseTargets(
            args.targets, args.genome, use_db, db, padSize, args.targetRegion, args.exons,
            args.targetUpstreamPromoter, args.targetDownstreamPromoter,
            helperFunctions.CONFIG["PATH"]["TWOBIT_INDEX_DIR"] if not ISOFORMS else CONFIG["PATH"]["ISOFORMS_INDEX_DIR"],
            args.outputDir, args.consensusUnion, args.jsonVisualize, args.guideSize)
        sequences, fastaSequence = mainFunctions.coordToFasta(
            targets, candidate_fasta_file, args.outputDir, args.guideSize, evalSequence, args.nonOverlapping,
            helperFunctions.CONFIG["PATH"]["TWOBIT_INDEX_DIR"] if not ISOFORMS else CONFIG["PATH"]["ISOFORMS_INDEX_DIR"],
            args.genome, strand, mainFunctions.DOWNSTREAM_NUC)

    # Converts genomic coordinates to fasta file of all possible k-mers
    if len(sequences) == 0:
        sys.stderr.write("No target sites\n")
        sys.exit()

    # Run bowtie and get results
    bowtieResultsFile = mainFunctions.runBowtie(len(args.PAM), args.uniqueMethod_Cong, candidate_fasta_file, args.outputDir,
                                                int(args.maxOffTargets),
                                                helperFunctions.CONFIG["PATH"]["ISOFORMS_INDEX_DIR"] if ISOFORMS else helperFunctions.CONFIG["PATH"][
                                      "BOWTIE_INDEX_DIR"],
                                                args.genome, int(args.maxMismatches))
    results = helperFunctions.parseBowtie(guideClass, bowtieResultsFile, True, args.scoreGC, args.scoreSelfComp,
                                          args.backbone, args.replace5P, args.maxOffTargets, countMM, args.PAM,
                                          args.MODE != mainFunctions.ProgramMode.TALENS,
                                          args.scoringMethod, args.genome, gene, isoform, gene_isoforms)  # TALENS: MAKE_PAIRS + CLUSTER

    # TODO this is a temporary fix, args.scoringMethod should be converted to type ScoringMethod like args.MODE
    scoring_method = scoring.ScoringMethod.G_20
    for sm in scoring.ScoringMethod:
        if args.scoringMethod == sm.name:
            scoring_method = sm
            break

    cluster_info = scoring.ClusterInfo(sequences, args.guideSize,
                                       args.taleMin if args.MODE == ProgramMode.TALENS else args.nickaseMin,
                                       args.taleMax if args.MODE == ProgramMode.TALENS else args.nickaseMax,
                                       args.enzymeCo, args.maxOffTargets, args.g_RVD, args.minResSiteLen,
                                       args.offtargetMaxDist)

    info = scoring.ScoringInfo(args.genome, args.PAM, strand, sortOutput, cluster_info, args.outputDir,
                               args.repairPredictions is not None, args.repairPredictions, args.isoforms,
                               visCoords, args.fasta, args.rm1perfOff, args.MODE, scoring_method)

    sorted_output, cluster = scoring.score_guides(results, info)

    # Write individual results to file
    listOfClusters = mainFunctions.writeIndividualResults(args.outputDir, args.maxOffTargets, sorted_output,
                                                          args.guideSize, args.MODE, cluster,
                                                          args.limitPrintResults, args.offtargetsTable)

    if args.makePrimers:
        if args.fasta:
            mainFunctions.make_primers_fasta(sorted_output, args.outputDir, args.primerFlanks,
                                             args.displaySeqFlanks, args.genome, args.limitPrintResults,
                                             CONFIG["PATH"]["BOWTIE_INDEX_DIR"], fastaSequence,
                                             args.primer3options, args.guidePadding, args.enzymeCo,
                                             args.minResSiteLen, "sequence", args.maxOffTargets)
        else:
            mainFunctions.make_primers_genome(sorted_output, args.outputDir, args.primerFlanks,
                                              args.displaySeqFlanks, args.genome, args.limitPrintResults,
                                              CONFIG["PATH"]["BOWTIE_INDEX_DIR"],
                                              CONFIG["PATH"]["TWOBIT_INDEX_DIR"] if not ISOFORMS
                                                         else CONFIG["PATH"]["ISOFORMS_INDEX_DIR"], args.primer3options,
                                              args.guidePadding, args.enzymeCo, args.minResSiteLen, strand,
                                              args.targets, args.maxOffTargets)

    #########- Print part -##########
    ## Print results
    print_scores(sorted_output, args.MODE, args.scoringMethod, args.isoforms)

    resultCoords = generate_result_coordinates(sorted_output,
                                               listOfClusters,
                                               sortOutput,
                                               args.MODE,
                                               args.isoforms)

    # Print gene annotation files
    # FASTA file
    geneFile = open('%s/gene_file.fa' % args.outputDir, 'w')
    geneFile.write(">%s\n" % args.targets)
    geneFile.write(fastaSequence)
    geneFile.close()


    # Visualize with json
    if args.jsonVisualize:

        # new function
        cutcoords = getCoordinatesForJsonVisualization(args, visCoords, sequences, strand, resultCoords)

        info = open("%s/run.info" % args.outputDir, 'w')
        info.write("%s\t%s\t%s\t%s\t%s\n" % ("".join(args.targets), args.genome, args.MODE, args.uniqueMethod_Cong,
                                             args.guideSize))
        info.close()

        if args.BED:
            mainFunctions.print_bed(args.MODE, visCoords, cutcoords, '%s/results.bed' % args.outputDir,
                                    visCoords[0]["name"] if args.fasta else args.targets)

        if args.GenBank:
            if args.fasta:
                seq = fastaSequence
                chrom = visCoords[0]["name"]
                start = 0
                finish = len(fastaSequence)
            else:
                # targets min-max (with introns)
                regions = targets
                chrom = regions[0][0:regions[0].rfind(':')]
                start = []
                finish = []
                targets = []
                for region in regions:
                    start_r = int(region[region.rfind(':') + 1:region.rfind('-')])
                    start_r = max(start_r, 0)
                    start.append(start_r)
                    finish_r = int(region[region.rfind('-') + 1:])
                    finish.append(finish_r)
                    targets.append([chrom, start_r, finish_r])
                start = min(start)
                finish = max(finish)

                prog = subprocess.Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2> %s/twoBitToFa.err" % (
                    CONFIG["PATH"]["TWOBITTOFA"], chrom, start, finish, CONFIG["PATH"]["TWOBIT_INDEX_DIR"] if not ISOFORMS else CONFIG["PATH"]["ISOFORMS_INDEX_DIR"],
                    args.genome, args.outputDir), stdout=subprocess.PIPE, shell=True)
                output = prog.communicate()
                if prog.returncode != 0:
                    sys.stderr.write("Running twoBitToFa failed when creating GenBank file\n")
                    sys.exit(EXIT['TWOBITTOFA_ERROR'])

                output = output[0]
                output = output.split("\n")
                seq = ''.join(output[1:]).upper()

            mainFunctions.print_genbank(args.MODE, chrom if args.fasta else args.targets, seq,
                                        [] if args.fasta else targets, cutcoords, chrom, start, finish,
                                        strand, '%s/results.gb' % args.outputDir, "CHOPCHOP results")


    # remove .sam files as they take up wayyy to much space
    for fl in os.listdir(args.outputDir):
        if fl.endswith(".sam"):
            os.remove(os.path.join(args.outputDir, fl))


if __name__ == '__main__':
    main()
