#!/usr/bin/env python3

import pandas
import numpy
import argparse
import resource
from Vars import *
from collections import defaultdict
from operator import itemgetter
from subprocess import Popen, PIPE

from classes.ProgramMode import ProgramModeAction
import functions.Main_Functions
from functions.Helper_Functions import *
from functions.TALEN_Specific_Functions import *

from dockers.CRISPRoff_wrapper import run_coefficient_score
from dockers.doench_2016_wrapper import run_doench_2016

soft, HARD_LIMIT = resource.getrlimit(resource.RLIMIT_NOFILE)
resource.setrlimit(resource.RLIMIT_NOFILE, (HARD_LIMIT, HARD_LIMIT))


def parse_arguments():
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

    parser.add_argument("-T", "--MODE", type=int, default=1, choices=[1, 2, 3, 4], action=ProgramModeAction,
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

    return parser.parse_args()


def getCoordinatesForJsonVisualization(args, visCoords, sequences, strand, resultCoords):
    # Coordinates for gene
    visCoordsFile = open('%s/viscoords.json' % args.outputDir, 'w')
    # visCoords = sorted(visCoords,  key=itemgetter(1))
    json.dump(visCoords, visCoordsFile)

    # Coordinates for sequence
    seqvis = functions.Main_Functions.FastaToViscoords(sequences, strand)
    seqvisFile = open('%s/seqviscoords.json' % args.outputDir, 'w')
    json.dump(seqvis, seqvisFile)

    # Coordinates for cutters
    cutCoord_file = open('%s/cutcoords.json' % args.outputDir, 'w')

    cutcoords = []
    for i in range(len(resultCoords)):
        el = []

        if args.MODE == functions.Main_Functions.ProgramMode.CRISPR or args.MODE == functions.Main_Functions.ProgramMode.CPF1:
            el.append(i + 1)
            el.extend(resultCoords[i])
        elif args.MODE == functions.Main_Functions.ProgramMode.TALENS or args.MODE == functions.Main_Functions.ProgramMode.NICKASE:
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


def scoringMethodCHARI_2005(args, results):
    try:
        # make file to score
        svmInputFile = '%s/chari_score.SVMInput.txt' % args.outputDir
        svmOutputFile = '%s/chari_score.SVMOutput.txt' % args.outputDir
        encoding = defaultdict(str)
        encoding['A'] = '0001'
        encoding['C'] = '0010'
        encoding['T'] = '0100'
        encoding['G'] = '1000'

        svmFile = open(svmInputFile, 'w')

        for guide in results:
            seq = guide.downstream5prim + guide.strandedGuideSeq[:-len(guide.PAM)]
            PAM = guide.strandedGuideSeq[-len(guide.PAM):]
            sequence = (seq[-20:] + PAM).upper()
            x = 0
            tw = '-1'
            # end index
            if len(sequence) == 27:
                endIndex = 22
            else:
                endIndex = 21

            while x < endIndex:
                y = 0
                while y < 4:
                    tw = tw + ' ' + str(x + 1) + str(y + 1) + ':' + encoding[sequence[x]][y]
                    y += 1
                x += 1
            svmFile.write(tw + '\n')
        svmFile.close()
        newScores = functions.Main_Functions.scoreChari_2015(svmInputFile, svmOutputFile, args.PAM, args.genome)

        for i, guide in enumerate(results):
            guide.CoefficientsScore["CHARI_2015"] = newScores[i]
            if args.scoringMethod == "CHARI_2015":
                guide.score -= (guide.CoefficientsScore["CHARI_2015"] / 100) * SCORE['COEFFICIENTS']
    except:
        pass


def scoringMethodZHANG_2009(args, results):
    try:
        zhangInputFile = '%s/zhang_score.txt' % args.outputDir
        zhangFile = open(zhangInputFile, 'w')
        for guide in results:
            zhangFile.write(guide.downstream5prim[-4:] + guide.strandedGuideSeq + guide.downstream3prim[:3] + '\n')
        zhangFile.close()

        prog = Popen("%s/uCRISPR/uCRISPR -on %s" % (f_p, zhangInputFile), stdout=PIPE, stderr=PIPE, shell=True)
        output = prog.communicate()
        output = output[0].splitlines()
        output = output[1:]
        # distribution calculated on 100k random guides
        output = [functions.Main_Functions.functions.Main_Functions.ss.norm.cdf(float(x.split()[1]), loc=11.92658, scale=0.2803797) for x in output]
        for i, guide in enumerate(results):
            guide.CoefficientsScore["ZHANG_2019"] = output[i] * 100
            if args.scoringMethod == "ZHANG_2019":
                guide.score -= (guide.CoefficientsScore["ZHANG_2019"] / 100) * SCORE['COEFFICIENTS']
    except:
        pass


def scoringMethodKIM_2018(results):
    # noinspection PyBroadException
    try:
        with functions.Main_Functions.warnings.catch_warnings(record=True):
            functions.Main_Functions.warnings.simplefilter("ignore")

            os.environ['KERAS_BACKEND'] = 'theano'
            stderr = sys.stderr  # keras prints welcome message to stderr! lolz!
            sys.stderr = open(os.devnull, 'w')

            from keras.models import Model
            from keras.layers import Input
            from keras.layers.merge import Multiply
            from keras.layers.core import Dense, Dropout, Activation, Flatten
            from keras.layers.convolutional import Convolution1D, AveragePooling1D
            sys.stderr = stderr

            seq_deep_cpf1_input_seq = Input(shape=(34, 4))
            seq_deep_cpf1_c1 = Convolution1D(80, 5, activation='relu')(seq_deep_cpf1_input_seq)
            seq_deep_cpf1_p1 = AveragePooling1D(2)(seq_deep_cpf1_c1)
            seq_deep_cpf1_f = Flatten()(seq_deep_cpf1_p1)
            seq_deep_cpf1_do1 = Dropout(0.3)(seq_deep_cpf1_f)
            seq_deep_cpf1_d1 = Dense(80, activation='relu')(seq_deep_cpf1_do1)
            seq_deep_cpf1_do2 = Dropout(0.3)(seq_deep_cpf1_d1)
            seq_deep_cpf1_d2 = Dense(40, activation='relu')(seq_deep_cpf1_do2)
            seq_deep_cpf1_do3 = Dropout(0.3)(seq_deep_cpf1_d2)
            seq_deep_cpf1_d3 = Dense(40, activation='relu')(seq_deep_cpf1_do3)
            seq_deep_cpf1_do4 = Dropout(0.3)(seq_deep_cpf1_d3)
            seq_deep_cpf1_output = Dense(1, activation='linear')(seq_deep_cpf1_do4)
            seq_deep_cpf1 = Model(inputs=[seq_deep_cpf1_input_seq], outputs=[seq_deep_cpf1_output])
            seq_deep_cpf1.load_weights(f_p + '/models/Seq_deepCpf1_weights.h5')

            # process data
            data_n = len(results)
            one_hot = numpy.zeros((data_n, 34, 4), dtype=int)

            for l in range(0, data_n):
                prim5 = results[l].downstream5prim[-4:]
                if len(prim5) < 4:  # cover weird genomic locations
                    prim5 = "N" * (4 - len(prim5)) + prim5
                guide_seq = results[l].strandedGuideSeq
                prim3 = results[l].downstream3prim[:6]
                if len(prim3) < 6:
                    prim5 = "N" * (6 - len(prim5)) + prim5
                seq = prim5 + guide_seq + prim3

                for i in range(34):
                    if seq[i] in "Aa":
                        one_hot[l, i, 0] = 1
                    elif seq[i] in "Cc":
                        one_hot[l, i, 1] = 1
                    elif seq[i] in "Gg":
                        one_hot[l, i, 2] = 1
                    elif seq[i] in "Tt":
                        one_hot[l, i, 3] = 1
                    elif seq[i] in "Nn":  # N will activate all nodes
                        one_hot[l, i, 0] = 1
                        one_hot[l, i, 1] = 1
                        one_hot[l, i, 2] = 1
                        one_hot[l, i, 3] = 1

            seq_deep_cpf1_score = seq_deep_cpf1.predict([one_hot], batch_size=50, verbose=0)

        for i, guide in enumerate(results):
            guide.CoefficientsScore = seq_deep_cpf1_score[i][0]
            guide.score -= (guide.CoefficientsScore / 100) * SCORE['COEFFICIENTS']
    except:
        pass

"""
def scoringMethodDOENCH_2016(args, results):
    # noinspection PyBroadException
    try:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("ignore")
            with open(f_p + '/models/Doench_2016_18.01_model_nopos.pickle', 'rb') as f:
                model = pickle.load(f)

        model, learn_options = model
        learn_options["V"] = 2

        results_ok = []
        sequences_d2016 = []
        for i, guide in enumerate(results):
            seq_d2016 = guide.downstream5prim + guide.strandedGuideSeq[:-len(guide.PAM)]
            pam_d2016 = guide.strandedGuideSeq[-len(guide.PAM):]
            tail_d2016 = guide.downstream3prim
            if len(seq_d2016) < 24 or len(pam_d2016) < 3 or len(tail_d2016) < 3:
                results_ok.append(False)
            else:
                results_ok.append(True)
                dada = seq_d2016[-24:] + pam_d2016 + tail_d2016[:3]
                sequences_d2016.append(dada)

        sequences_d2016 = numpy.array(sequences_d2016)
        xdf = pandas.DataFrame(columns=[u'30mer', u'Strand'],
                               data=zip(sequences_d2016, numpy.repeat('NA', sequences_d2016.shape[0])))
        gene_position = pandas.DataFrame(columns=[u'Percent Peptide', u'Amino Acid Cut position'],
                                         data=zip(numpy.ones(sequences_d2016.shape[0]) * -1,
                                                  numpy.ones(sequences_d2016.shape[0]) * -1))
        feature_sets = feat.featurize_data(xdf, learn_options, pandas.DataFrame(), gene_position, pam_audit=True,
                                           length_audit=False)
        inputs = concatenate_feature_sets(feature_sets)[0]
        outputs = model.predict(inputs)

        j = 0
        for i, guide in enumerate(results):
            if results_ok[i]:
                if outputs[j] > 1:
                    outputs[j] = 1
                elif outputs[j] < 0:
                    outputs[j] = 0
                guide.CoefficientsScore["DOENCH_2016"] = outputs[j] * 100
                j += 1
                if args.scoringMethod == "DOENCH_2016":
                    guide.score -= (guide.CoefficientsScore["DOENCH_2016"] / 100) * SCORE['COEFFICIENTS']
    except:
        pass


def getClusterPairsTALENS(results, sequences, args):
    pairs = pairTalens(results, sequences, args.guideSize, int(args.taleMin), int(args.taleMax), args.enzymeCo,
                       args.maxOffTargets, args.g_RVD, args.minResSiteLen)

    if (not len(pairs)):
        sys.stderr.write("No TALEN pairs could be generated for this region.\n")
        sys.exit(EXIT['GENE_ERROR'])

    if args.rm1perfOff and args.fasta:
        for pair in pairs:
            if pair.diffStrandOffTarget > 0:
                pair.score = pair.score - SCORE["OFFTARGET_PAIR_DIFF_STRAND"]
            if pair.sameStrandOffTarget > 0:
                pair.score = pair.score - SCORE["OFFTARGET_PAIR_SAME_STRAND"]

    cluster, results = clusterPairs(pairs)
    return cluster, results


def getClusterPairsNICKASE(results, sequences, args):
    pairs = pairCas9(results, sequences, args.guideSize, int(args.nickaseMin), int(args.nickaseMax), args.enzymeCo,
                     args.maxOffTargets, args.minResSiteLen, args.offtargetMaxDist)

    if (not len(pairs)):
        sys.stderr.write("No Cas9 nickase pairs could be generated for this region.\n")
        sys.exit(EXIT['GENE_ERROR'])

    if args.rm1perfOff and args.fasta:
        for pair in pairs:
            if pair.diffStrandOffTarget > 0:
                pair.score = pair.score - SCORE["OFFTARGET_PAIR_DIFF_STRAND"]

    cluster, results = clusterPairs(pairs)
    return cluster, results
"""

def main():
    # Parse arguments
    args = parse_arguments()

    # set isoforms to global as it is influencing many steps
    global ISOFORMS
    ISOFORMS = args.isoforms

    # Add TALEN length
    args.taleMin += 18
    args.taleMax += 18

    # Set mode specific parameters if not set by user
    args.scoreGC = functions.Main_Functions.mode_select(args.scoreGC, "SCORE_GC", args.MODE)
    args.scoreSelfComp = functions.Main_Functions.mode_select(args.noScoreSelfComp, "SCORE_FOLDING", args.MODE)
    args.PAM = functions.Main_Functions.mode_select(args.PAM, "PAM", args.MODE)
    args.guideSize = functions.Main_Functions.mode_select(args.guideSize, "GUIDE_SIZE", args.MODE) + len(args.PAM)
    args.maxMismatches = functions.Main_Functions.mode_select(args.maxMismatches, "MAX_MISMATCHES", args.MODE)
    args.maxOffTargets = functions.Main_Functions.mode_select(args.maxOffTargets, "MAX_OFFTARGETS", args.MODE)

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
        if args.MODE == functions.Main_Functions.ProgramMode.TALENS:
            padSize = args.taleMax
        elif args.MODE == functions.Main_Functions.ProgramMode.NICKASE:
            padSize = args.nickaseMax
        elif args.MODE == functions.Main_Functions.ProgramMode.CRISPR or args.MODE == functions.Main_Functions.ProgramMode.CPF1:
            padSize = args.guideSize

    # Set default functions for different modes
    # new function
    countMM, evalSequence, guideClass, sortOutput = functions.Main_Functions.set_default_modes(args)

    # Connect to database if requested
    if args.database:
        cdb = functions.Main_Functions.connect_db(args.database)
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
        sequences, targets, visCoords, fastaSequence, strand = functions.Main_Functions.parseFastaTarget(
            args.targets, candidate_fasta_file, args.guideSize, evalSequence)
    else:
        targets, visCoords, strand, gene, isoform, gene_isoforms = functions.Main_Functions.parseTargets(
            args.targets, args.genome, use_db, db, padSize, args.targetRegion, args.exons,
            args.targetUpstreamPromoter, args.targetDownstreamPromoter,
            CONFIG["PATH"]["TWOBIT_INDEX_DIR"] if not ISOFORMS else CONFIG["PATH"]["ISOFORMS_INDEX_DIR"],
            args.outputDir, args.consensusUnion, args.jsonVisualize, args.guideSize)
        sequences, fastaSequence = functions.Main_Functions.coordToFasta(
            targets, candidate_fasta_file, args.outputDir, args.guideSize, evalSequence, args.nonOverlapping,
            CONFIG["PATH"]["TWOBIT_INDEX_DIR"] if not ISOFORMS else CONFIG["PATH"]["ISOFORMS_INDEX_DIR"],
            args.genome, strand, functions.Main_Functions.DOWNSTREAM_NUC)

    # Converts genomic coordinates to fasta file of all possible k-mers
    if len(sequences) == 0:
        sys.stderr.write("No target sites\n")
        sys.exit()

    # Run bowtie and get results
    bowtieResultsFile = functions.Main_Functions.runBowtie(len(args.PAM), args.uniqueMethod_Cong, candidate_fasta_file, args.outputDir,
                                                           int(args.maxOffTargets),
                                                           CONFIG["PATH"]["ISOFORMS_INDEX_DIR"] if ISOFORMS else CONFIG["PATH"][
                                      "BOWTIE_INDEX_DIR"],
                                                           args.genome, int(args.maxMismatches))
    results = parseBowtie(guideClass, bowtieResultsFile, True, args.scoreGC, args.scoreSelfComp,
                          args.backbone, args.replace5P, args.maxOffTargets, countMM, args.PAM,
                          args.MODE != functions.Main_Functions.ProgramMode.TALENS,
                          args.scoringMethod, args.genome, gene, isoform, gene_isoforms)  # TALENS: MAKE_PAIRS + CLUSTER



    ########## -Scoring part- ###########
    if args.rm1perfOff and args.fasta:
        for guide in results:
            if guide.offTargetsMM[0] > 0:
                guide.score -= SINGLE_OFFTARGET_SCORE[0]

    if ISOFORMS:
        for guide in results:
            if guide.isoform in ["union", "intersection"]: # calculate base pair probabilities of folding
                # iterate all isoforms
                bpp = []
                for tx_id in guide.gene_isoforms:
                    tx_start, tx_end = functions.Main_Functions.tx_relative_coordinates(visCoords, tx_id, guide.start, guide.end)
                    if tx_start != -1:
                        bpp.append(functions.Main_Functions.rna_folding_metric(args.genome, tx_id, tx_start, tx_end))
                guide.meanBPP = 100 if len(bpp) == 0 else max(bpp) # penalize guide that has no real target!
            else:
                if not args.fasta:
                    tx_start, tx_end = functions.Main_Functions.tx_relative_coordinates(visCoords, guide.isoform, guide.start, guide.end)
                    guide.meanBPP = functions.Main_Functions.rna_folding_metric(args.genome, guide.isoform, tx_start, tx_end)

            guide.score += guide.meanBPP / 100 * SCORE['COEFFICIENTS']

            if guide.isoform in guide.gene_isoforms:
                guide.gene_isoforms.remove(guide.isoform)

            if guide.isoform in guide.offTargetsIso[0]:
                guide.offTargetsIso[0].remove(guide.isoform)

            guide.constitutive = int(guide.gene_isoforms == guide.offTargetsIso[0])

    # Scoring methods
    if (args.scoringMethod == "CHARI_2015" or args.scoringMethod == "ALL") and (args.PAM == "NGG" or args.PAM == "NNAGAAW") and (args.genome == "hg19" or args.genome == "mm10") and not ISOFORMS:
        # new function
        scoringMethodCHARI_2005(args, results)

    if (args.scoringMethod == "ZHANG_2019" or args.scoringMethod == "ALL") and (args.PAM == "NGG") and not ISOFORMS:
        # new function
        scoringMethodZHANG_2009(args, results)

    if (args.scoringMethod == "KIM_2018" or args.scoringMethod == "ALL") and args.PAM in "TTTN" \
            and not ISOFORMS and args.MODE == functions.Main_Functions.ProgramMode.CPF1:
        # new function
        scoringMethodKIM_2018(results)

    if (args.scoringMethod == "DOENCH_2016" or args.scoringMethod == "ALL") and not ISOFORMS and args.MODE == functions.Main_Functions.ProgramMode.CRISPR:
        guides = run_doench_2016(args.scoringMethod, results)
        # new function
        #scoringMethodDOENCH_2016(args, results)

    if args.repairPredictions is not None and not ISOFORMS and args.MODE == functions.Main_Functions.ProgramMode.CRISPR:
        sys.path.append(f_p + '/models/inDelphi-model/')
        with functions.Main_Functions.warnings.catch_warnings(record=True):
            functions.Main_Functions.warnings.simplefilter("ignore")
            import inDelphi
            inDelphi.init_model(celltype=args.repairPredictions)
            for i, guide in enumerate(results):
                # noinspection PyBroadException
                try:
                    left_seq = guide.downstream5prim + guide.strandedGuideSeq[:-(len(guide.PAM) + 3)]
                    left_seq = left_seq[-60:]
                    right_seq = guide.strandedGuideSeq[-(len(guide.PAM) + 3):] + guide.downstream3prim
                    right_seq = right_seq[:60]
                    seq = left_seq + right_seq
                    cutsite = len(left_seq)
                    pred_df, stats = inDelphi.predict(seq, cutsite)
                    pred_df = pred_df.sort_values(pred_df.columns[4], ascending=False)
                    guide.repProfile = pred_df
                    guide.repStats = stats
                except:
                    pass


    if args.MODE == functions.Main_Functions.ProgramMode.CRISPR or args.MODE == functions.Main_Functions.ProgramMode.CPF1 or ISOFORMS:
        cluster = 0
    elif args.MODE == functions.Main_Functions.ProgramMode.TALENS:
        pairs = pairTalens(results, sequences, args.guideSize, int(args.taleMin), int(args.taleMax), args.enzymeCo,
                           args.maxOffTargets, args.g_RVD, args.minResSiteLen)

        if (not len(pairs)):
            sys.stderr.write("No TALEN pairs could be generated for this region.\n")
            sys.exit(EXIT['GENE_ERROR'])

        if args.rm1perfOff and args.fasta:
            for pair in pairs:
                if pair.diffStrandOffTarget > 0:
                    pair.score = pair.score - SCORE["OFFTARGET_PAIR_DIFF_STRAND"]
                if pair.sameStrandOffTarget > 0:
                    pair.score = pair.score - SCORE["OFFTARGET_PAIR_SAME_STRAND"]

        cluster, results = functions.Main_Functions.clusterPairs(pairs)
        return cluster, results

    elif args.MODE == functions.Main_Functions.ProgramMode.NICKASE:
        pairs = pairCas9(results, sequences, args.guideSize, int(args.nickaseMin), int(args.nickaseMax), args.enzymeCo, args.maxOffTargets, args.minResSiteLen, args.offtargetMaxDist)

        if (not len(pairs)):
            sys.stderr.write("No Cas9 nickase pairs could be generated for this region.\n")
            sys.exit(EXIT['GENE_ERROR'])

        if args.rm1perfOff and args.fasta:
            for pair in pairs:
                if pair.diffStrandOffTarget > 0:
                    pair.score = pair.score - SCORE["OFFTARGET_PAIR_DIFF_STRAND"]

        cluster, results = functions.Main_Functions.clusterPairs(pairs)

    # Sorts pairs according to score/penalty and cluster
    if strand == "-" and not ISOFORMS:
        results.reverse()

    sortedOutput = sortOutput(results)

    # Write individual results to file
    listOfClusters = functions.Main_Functions.writeIndividualResults(args.outputDir, args.maxOffTargets, sortedOutput, args.guideSize, args.MODE, cluster, args.limitPrintResults, args.offtargetsTable)

    if args.makePrimers:
        if args.fasta:
            functions.Main_Functions.make_primers_fasta(sortedOutput, args.outputDir, args.primerFlanks, args.displaySeqFlanks, args.genome, args.limitPrintResults, CONFIG["PATH"]["BOWTIE_INDEX_DIR"], fastaSequence, args.primer3options, args.guidePadding, args.enzymeCo, args.minResSiteLen, "sequence", args.maxOffTargets)
        else:
            functions.Main_Functions.make_primers_genome(sortedOutput, args.outputDir, args.primerFlanks, args.displaySeqFlanks, args.genome, args.limitPrintResults, CONFIG["PATH"]["BOWTIE_INDEX_DIR"], CONFIG["PATH"]["TWOBIT_INDEX_DIR"] if not ISOFORMS else CONFIG["PATH"]["ISOFORMS_INDEX_DIR"], args.primer3options, args.guidePadding, args.enzymeCo, args.minResSiteLen, strand, args.targets, args.maxOffTargets)




    #########- Print part -##########
    ## Print results
    resultCoords = []

    if ISOFORMS:
        print ("Rank\tTarget sequence\tGenomic location\tGene\tIsoform\tGC content (%)\tSelf-complementarity\tLocal structure\tMM0\tMM1\tMM2\tMM3\tConstitutive\tIsoformsMM0\tIsoformsMM1\tIsoformsMM2\tIsoformsMM3")
        for i in range(len(sortedOutput)):
            print ("%s\t%s" % (i+1, sortedOutput[i]))
            resultCoords.append([sortedOutput[i].start, sortedOutput[i].score, sortedOutput[i].guideSize, sortedOutput[i].strand])
    else:
        if args.MODE == functions.Main_Functions.ProgramMode.CRISPR:
            common_header = "Rank\tTarget sequence\tGenomic location\tStrand\tGC content (%)\tSelf-complementarity\tMM0\tMM1\tMM2\tMM3"
            if args.scoringMethod == "ALL":
                print(common_header + "\tXU_2015\tDOENCH_2014\tDOENCH_2016\tMORENO_MATEOS_2015\tCHARI_2015\tG_20\tALKAN_2018\tZHANG_2019")
            else:
                print(common_header + "\tEfficiency")
            for i in range(len(sortedOutput)):
                print ("%s\t%s" % (i+1, sortedOutput[i]))
                resultCoords.append([sortedOutput[i].start, sortedOutput[i].score, sortedOutput[i].guideSize, sortedOutput[i].strand])

        elif args.MODE == functions.Main_Functions.ProgramMode.CPF1:
            print ("Rank\tTarget sequence\tGenomic location\tStrand\tGC content (%)\tSelf-complementarity\tEfficiency\tMM0\tMM1\tMM2\tMM3")
            for i in range(len(sortedOutput)):
                print ("%s\t%s" % (i+1, sortedOutput[i]))
                resultCoords.append([sortedOutput[i].start, sortedOutput[i].score, sortedOutput[i].guideSize, sortedOutput[i].strand])

        elif args.MODE == functions.Main_Functions.ProgramMode.TALENS or args.MODE == functions.Main_Functions.ProgramMode.NICKASE:

            if args.MODE == functions.Main_Functions.ProgramMode.TALENS:
                print ("Rank\tTarget sequence\tGenomic location\tTALE 1\tTALE 2\tCluster\tOff-target pairs\tOff-targets MM0\tOff-targets MM1\tOff-targets MM2\tOff-targets MM3\tRestriction sites\tBest ID")
            else:
                print ("Rank\tTarget sequence\tGenomic location\tCluster\tOff-target pairs\tOff-targets MM0\tOff-targets MM1\tOff-targets MM2\tOff-targets MM3\tRestriction sites\tBest ID")
            finalOutput = []
            for cluster in listOfClusters:  ## FIX: WHY ARE THERE EMPTY CLUSTERS???
                if len(cluster) == 0:
                    continue

                finalOutput.append(cluster[0])

            sortedFinalOutput = sortOutput(finalOutput)
            resultCoords = [[j+1, sortedFinalOutput[j].spacerStart, sortedFinalOutput[j].score, sortedFinalOutput[j].spacerSize, sortedFinalOutput[j].strand, sortedFinalOutput[j].ID, sortedFinalOutput[j].tale1.start, sortedFinalOutput[j].tale2.end] for j in range(len(sortedFinalOutput))]

            for i in range(len(sortedFinalOutput)):
                print ("%s\t%s\t%s" % (i+1,sortedFinalOutput[i], sortedFinalOutput[i].ID))

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
            functions.Main_Functions.print_bed(args.MODE, visCoords, cutcoords, '%s/results.bed' % args.outputDir,
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

                prog = Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2> %s/twoBitToFa.err" % (
                    CONFIG["PATH"]["TWOBITTOFA"], chrom, start, finish, CONFIG["PATH"]["TWOBIT_INDEX_DIR"] if not ISOFORMS else CONFIG["PATH"]["ISOFORMS_INDEX_DIR"],
                    args.genome, args.outputDir), stdout=PIPE, shell=True)
                output = prog.communicate()
                if prog.returncode != 0:
                    sys.stderr.write("Running twoBitToFa failed when creating GenBank file\n")
                    sys.exit(EXIT['TWOBITTOFA_ERROR'])

                output = output[0]
                output = output.split("\n")
                seq = ''.join(output[1:]).upper()

            functions.Main_Functions.print_genbank(args.MODE, chrom if args.fasta else args.targets, seq,
                                                   [] if args.fasta else targets, cutcoords, chrom, start, finish,
                                                   strand, '%s/results.gb' % args.outputDir, "CHOPCHOP results")


    # remove .sam files as they take up wayyy to much space
    for fl in os.listdir(args.outputDir):
        if fl.endswith(".sam"):
            os.remove(os.path.join(args.outputDir, fl))


if __name__ == '__main__':
    main()
