import argparse
import logging

from Bio.Seq import Seq

import config
from classes.ProgramMode import ProgramMode
from functions.main_functions import mode_select


def generate_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-Target", "--targets",
                        type=str, required=True,
                        help="Target genes or regions")

    parser.add_argument("-r", "--gRVD",
                        action="store_const", const="NN ", default="NH ", dest="g_rvd",
                        help="Use RVD 'NN' instead of 'NH' for guanine nucleotides. 'NH' appears to be more specific "
                             "than 'NN' but the choice depends on assembly kit.")

    parser.add_argument("-D", "--database",
                        metavar="DATABASE", dest="database",
                        help="Connect to a chopchop database to retrieve gene: user_name:passwd@host/database")

    parser.add_argument("-e", "--exon",
                        metavar="EXON_NUMBER", dest="exons",
                        help="Comma separated list of exon indices. Only find sites in this subset. ")

    parser.add_argument("-TDP", "--targetDownstreamPromoter",
                        type=int, default=200, dest="target_downstream_promoter",
                        help="how many bp to target downstream of TSS")

    parser.add_argument("-TUP", "--targetUpstreamPromoter",
                        type=int, default=200, dest="target_upstream_promoter",
                        help="how many bp to target upstream of TSS")

    parser.add_argument("-G", "--genome",
                        default="danRer7", metavar="GENOME", dest="genome",
                        help="The genome to search.")

    parser.add_argument("-g", "--guideSize",
                        type=int, default=None, metavar="GUIDE_SIZE", dest="guide_size",
                        help="The size of the guide RNA.")

    parser.add_argument("-c", "--scoreGC",
                        default=None, action="store_false", dest="score_gc",
                        help="Score GC content. True for CRISPR, False for TALENs.")

    parser.add_argument("-SC", "--noScoreSelfComp",
                        default=None, action="store_false", dest="no_score_self_comp",
                        help="Do not penalize self-complementarity of CRISPR.")

    parser.add_argument("-BB", "--backbone",
                        type=str, default=None, dest="backbone",
                        help="Penalize self-complementarity versus backbone regions (comma-separated list, same strand "
                             "as guide). Requires -C.")

    # TODO FIX: AT THE MOMENT THIS IS ONLY APPLIES TO FOLDING/SELF-COMPL
    parser.add_argument("-R5", "--replace5P",
                        default=None, metavar="REPLACE_5P", dest="replace_5_prime",
                        help="Replace bases from 5' end (with e.g. 'GG') ")

    parser.add_argument("-t", "--target",
                        default="CODING", dest="target_region",
                        help="Target the whole gene CODING/WHOLE/UTR5/UTR3/SPLICE. Default is CODING.")

    parser.add_argument("-T", "--MODE",
                        type=int, default=1, choices=[1, 2, 3, 4], dest="program_mode",
                        help="Set mode (int): default is Cas9 = 1, Talen = 2, Cpf1 = 3, Nickase = 4")

    # 14 + 18(length of TALE) = 32
    parser.add_argument("-taleMin", "--taleMin",
                        type=int, default=14, dest="tale_min",
                        help="Minimum distance between TALENs. Default is 14.")

    # 20 + 18(length of TALE) = 38
    parser.add_argument("-taleMax", "--taleMax",
                        type=int, default=20, dest="tale_max",
                        help="Maximum distance between TALENs. Default is 20.")

    parser.add_argument("-nickaseMin", "--nickaseMin",
                        type=int, default=10, dest="nickase_min",
                        help="Minimum distance between TALENs. Default is 10.")

    parser.add_argument("-nickaseMax", "--nickaseMax",
                        type=int, default=31, dest="nickase_max",
                        help="Maximum distance between TALENs. Default is 31.")

    parser.add_argument("-offtargetMaxDist", "--offtargetMaxDist",
                        type=int, default=100, dest="off_target_max_dist",
                        help="Maximum distance between offtargets for Nickase. Default is 100.")

    parser.add_argument("-f", "--fivePrimeEnd",
                        type=str, default="NN", dest="five_prime_end",
                        help="Specifies the requirement of the two nucleotides 5' end of the CRISPR guide: A/C/G/T/N. "
                             "Default: NN.")

    parser.add_argument("-n", "--enzymeCo",
                        default="N", metavar="ENZYME_CO", dest="enzyme_co",
                        help="The restriction enzyme company for TALEN spacer.")

    parser.add_argument("-R", "--minResSiteLen",
                        type=int, default=4, dest="min_res_site_len",
                        help="The minimum length of the restriction enzyme.")

    parser.add_argument("-v", "--maxMismatches",
                        type=int, choices=[0, 1, 2, 3], default=3, metavar="MAX_MISMATCHES", dest="max_mismatches",
                        help="The number of mismatches to check across the sequence.")

    parser.add_argument("-m", "--maxOffTargets",
                        metavar="MAX_HITS", dest="max_off_targets",
                        help="The maximum number of off targets allowed.")

    parser.add_argument("-M", "--PAM",
                        type=str, dest="pam",
                        help="The PAM motif.")

    parser.add_argument("-o", "--outputDir",
                        default="./", metavar="OUTPUT_DIR", dest="output_dir",
                        help="The output directory. Default is the current directory.")

    parser.add_argument("-F", "--fasta",
                        default=False, action="store_true", dest="fasta",
                        help="Use FASTA file as input rather than gene or genomic region.")

    parser.add_argument("-p", "--padSize",
                        type=int, default=-1, dest="pad_size",
                        help="Extra bases searched outside the exon. Defaults to the size of the guide RNA for CRISPR "
                             "and TALEN + maximum spacer for TALENS.")

    parser.add_argument("-P", "--makePrimers",
                        default=False, action="store_true", dest="make_primers",
                        help="Designds primers using Primer3 to detect mutation.")

    parser.add_argument("-3", "--primer3options",
                        default=None, dest="primer_3_options",
                        help="Options for Primer3. E.g. 'KEY1=VALUE1,KEY2=VALUE2'")

    parser.add_argument("-A", "--primerFlanks",
                        type=int, default=300, dest="primer_flanks",
                        help="Size of flanking regions to search for primers.")

    parser.add_argument("-DF", "--displaySeqFlanks",
                        type=int, default=300, dest="display_seq_flanks",
                        help="Size of flanking regions to output sequence into locusSeq_.")

    parser.add_argument("-a", "--guidePadding",
                        type=int, default=20, dest="guide_padding",
                        help="Minimum distance of primer to target site.")

    parser.add_argument("-O", "--limitPrintResults",
                        type=int, default=(3000 if config.limits.hard() > 3000 else config.limits.hard()),
                        dest="limit_print_results",
                        help="The number of results to print extended information for. Web server can handle 4k of "
                             "these.")

    parser.add_argument("-w", "--uniqueMethod_Cong",
                        default=False, action="store_true", dest="unique_method_cong",
                        help="A method to determine how unique the site is in the genome: allows 0 mismatches in last "
                             "15 bp.")

    parser.add_argument("-J", "--jsonVisualize",
                        default=False, action="store_true", dest="json_visualize",
                        help="Create files for visualization with json.")

    parser.add_argument("-nonO", "--nonOverlapping",
                        default=False, action="store_true", dest="non_overlapping",
                        help="Will not produce overlapping guides, saves time, and recommended for permissive PAMs ("
                             "e.g. Cas13d).")

    parser.add_argument("-scoringMethod", "--scoringMethod",
                        type=str, default="G_20", dest="scoring_method",
                        choices=["XU_2015", "DOENCH_2014", "DOENCH_2016", "MORENO_MATEOS_2015", "CHARI_2015", "G_20",
                                 "KIM_2018", "ALKAN_2018", "ZHANG_2019", "ALL"],
                        help="Scoring used for Cas9 and Nickase. Default is G_20. If a method fails to give scores, "
                             "CHOPCHOP will output 0 instead of terminating.")

    parser.add_argument("-repairPredictions", "--repairPredictions",
                        type=str, default=None, dest="repair_predictions",
                        choices=['mESC', 'U2OS', 'HEK293', 'HCT116', 'K562'],
                        help="Use inDelphi from Shen et al 2018 to predict repair profiles for every guideRNA, "
                             "this will make .repProfile and .repStats files")

    parser.add_argument("-rm1perfOff", "--rm1perfOff",
                        default=False, action="store_true", dest="rm_1_perf_off",
                        help="For fasta input, don't score one off-target without mismatches.")

    parser.add_argument("-isoforms", "--isoforms",
                        default=False, action="store_true", dest="isoforms",
                        help="Search for offtargets on the transcriptome.")

    parser.add_argument("-filterGCmin", "--filterGCmin",
                        type=int, default=0, dest="filter_gc_min",
                        help="Minimum required GC percentage. Default is 0.")

    parser.add_argument("-filterGCmax", "--filterGCmax",
                        type=int, default=100, dest="filter_gc_max",
                        help="Maximum allowed GC percentage. Default is 100.")

    parser.add_argument("-filterSelfCompMax", "--filterSelfCompMax",
                        type=int, default=-1, dest="filter_self_comp_max",
                        help="Maximum acceptable Self-complementarity score. Default is -1, no filter.")

    parser.add_argument("-consensusUnion", "--consensusUnion",
                        default=False, action="store_true", dest="consensus_union",
                        help="When calculating consensus sequence from multiple isoforms default uses intersection. "
                             "This option specifies union of isoforms.")

    parser.add_argument("-BED", "--BED",
                        default=False, action="store_true", dest="bed",
                        help="Create results as BED file, can be used for integration with UCSC.")

    parser.add_argument("-GenBank", "--GenBank",
                        default=False, action="store_true", dest="gen_bank",
                        help="Create results as GenBank file, sequence of targeted region with introns is included.")

    parser.add_argument("-offtargetsTable", "--offtargetsTable",
                        default=False, action="store_true", dest="off_targets_table",
                        help="Create .tsv table with off-targets. Not all off-targets will be reported when early "
                             "stopping will work on a guide! Limited also to CRISPR mode only and limited by "
                             "--limitPrintResults option.")

    parser.add_argument("--logLevel",
                        type=str, default="ERROR", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        dest="log_level",
                        help="Set logging level.")

    return parser


def parse_arguments() -> argparse.Namespace:
    args = generate_parser().parse_args()

    # Logging
    logging.basicConfig(level=logging.getLevelName(args.log_level.upper()),
                        format="%(levelname)-8s:: %(message)s")
    logging.debug("Log level set to %s." % args.log_level)

    # Change args.MODE type from int to ProgramMode
    args.program_mode = ProgramMode(args.program_mode)

    # Add TALEN length
    args.tale_min += 18
    args.tale_max += 18

    # Set mode specific parameters if not set by user
    args.score_gc = mode_select(args.score_gc, "SCORE_GC", args.program_mode)
    args.score_self_comp = mode_select(args.no_score_self_comp, "SCORE_FOLDING", args.program_mode)
    args.pam = mode_select(args.pam, "PAM", args.program_mode)
    args.guide_size = mode_select(args.guide_size, "GUIDE_SIZE", args.program_mode) + len(args.pam)
    args.max_mismatches = mode_select(args.max_mismatches, "MAX_MISMATCHES", args.program_mode)
    args.max_off_targets = mode_select(args.max_off_targets, "MAX_OFFTARGETS", args.program_mode)

    # Add TALEN length
    args.nickase_min += args.guide_size
    args.nickase_max += args.guide_size

    if args.score_self_comp:
        if args.backbone:
            tmp = args.backbone.strip().split(",")
            args.backbone = [str(Seq(el).reverse_complement()) for el in tmp]
        else:
            args.backbone = []

    logging.debug("Finished parsing arguments.")

    return args
