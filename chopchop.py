#!/usr/bin/env python3
import json
import os
import subprocess
import sys
from operator import itemgetter
from typing import List, Callable, Union

import config
import scoring
from classes.Guide import Guide
from classes.PAIR import Pair
from classes.ProgramMode import ProgramMode
from constants import *
from functions.arguments import parse_arguments
from functions.main_functions import coord_to_fasta, run_bowtie
from functions.main_functions import print_bed, print_genbank, fasta_to_viscoords
from functions.main_functions import write_individual_results, parse_fasta_target, connect_db
from functions.make_primers import make_primers_fasta, make_primers_genome, parse_bowtie
from functions.parse_target import parse_targets
from functions.set_default_modes import set_default_modes


def get_coordinates_for_json_visualization(args, vis_coords, sequences, strand, result_coords):
    # Coordinates for gene
    vis_coords_file = open('%s/viscoords.json' % args.outputDir, 'w')
    # visCoords = sorted(visCoords,  key=itemgetter(1))
    json.dump(vis_coords, vis_coords_file)

    # Coordinates for sequence
    seqvis = fasta_to_viscoords(sequences, strand)
    seqvis_file = open('%s/seqviscoords.json' % args.outputDir, 'w')
    # TODO dont use hack fix with list her
    json.dump(list(seqvis), seqvis_file)

    # Coordinates for cutters
    cut_coord_file = open('%s/cutcoords.json' % args.outputDir, 'w')

    cutcoords = []
    for i in range(len(result_coords)):
        el = []

        if args.MODE == ProgramMode.CRISPR or args.MODE == ProgramMode.CPF1:
            el.append(i + 1)
            el.extend(result_coords[i])
        elif args.MODE == ProgramMode.TALENS or args.MODE == ProgramMode.NICKASE:
            el.extend(result_coords[i])

        cutcoords.append(el)

    # Put bars at different heights to avoid overlap
    tiers = [0] * 23
    sorted_coords = sorted(cutcoords, key=itemgetter(1))
    for coord in sorted_coords:
        t = 0
        for j in range(len(tiers)):
            if coord[1] > tiers[j]:
                t = j
                tiers[j] = coord[1] + coord[3]
                break

        coord.append(t)

    json.dump(cutcoords, cut_coord_file)

    return cutcoords


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
            if len(cluster) == 0:
                continue
            final_output.append(cluster[0])

        sorted_output = sort_function(final_output)
        for i, guide in enumerate(sorted_output):
            result_coordinates.append([i + 1, guide.spacerStart, guide.score, guide.spacerSize, guide.strand, guide.ID,
                                       guide.tale1.start, guide.tale2.end])

    return result_coordinates


def visualize_with_json(args, vis_coords, sequences, strand, result_coords, fasta_sequence, targets):
    # new function
    cutcoords = get_coordinates_for_json_visualization(args, vis_coords, sequences, strand, result_coords)

    info = open("%s/run.info" % args.outputDir, 'w')
    info.write("%s\t%s\t%s\t%s\t%s\n" % ("".join(args.targets), args.genome, args.MODE, args.uniqueMethod_Cong,
                                         args.guideSize))
    info.close()

    if args.BED:
        print_bed(args.MODE, vis_coords, cutcoords, '%s/results.bed' % args.outputDir,
                  vis_coords[0]["name"] if args.fasta else args.targets)

    if args.GenBank:
        if args.fasta:
            seq = fasta_sequence
            chrom = vis_coords[0]["name"]
            start = 0
            finish = len(fasta_sequence)
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
                config.path("TWOBITTOFA"), chrom, start, finish,
                config.path("TWOBIT_INDEX_DIR") if not config.isoforms else config.path("ISOFORMS_INDEX_DIR"),
                args.genome, args.outputDir), stdout=subprocess.PIPE, shell=True)
            output = prog.communicate()
            if prog.returncode != 0:
                sys.stderr.write("Running twoBitToFa failed when creating GenBank file\n")
                sys.exit(EXIT['TWOBITTOFA_ERROR'])

            output = output[0].decode()
            output = output.split("\n")
            seq = ''.join(output[1:]).upper()

        print_genbank(args.MODE, chrom if args.fasta else args.targets, seq,
                      [] if args.fasta else targets, cutcoords, chrom, start, finish,
                      strand, '%s/results.gb' % args.outputDir, "CHOPCHOP results")


def main():
    # Parse arguments
    args = parse_arguments()

    # set isoforms to global as it is influencing many steps
    config.isoforms = args.isoforms

    # Pad each exon equal to guidesize unless
    if args.padSize != -1:
        pad_size = args.padSize
    else:
        if args.MODE == ProgramMode.TALENS:
            pad_size = args.taleMax
        elif args.MODE == ProgramMode.NICKASE:
            pad_size = args.nickaseMax
        elif args.MODE == ProgramMode.CRISPR or args.MODE == ProgramMode.CPF1:
            pad_size = args.guideSize

    # Set default functions for different modes
    # new function
    count_MM, eval_sequence, guide_class, sort_output = set_default_modes(args)

    # General ARGPARSE done, upcoming Target parsing

    # Connect to database if requested
    if args.database:
        cdb = connect_db(args.database)
        db = cdb.cursor()
        use_db = True
    else:
        db = "%s/%s.gene_table" % (
            config.path("GENE_TABLE_INDEX_DIR") if not config.isoforms else config.path("ISOFORMS_INDEX_DIR"),
            args.genome)
        use_db = False

    # Create output directory if it doesn't exist
    if not os.path.isdir(args.outputDir):
        os.mkdir(args.outputDir)

    candidate_fasta_file = '%s/sequence.fa' % args.outputDir
    gene, isoform, gene_isoforms = (None, None, set())
    if args.fasta:
        sequences, targets, vis_coords, fasta_sequence, strand = parse_fasta_target(
            args.targets, candidate_fasta_file, args.guideSize, eval_sequence)

    else:
        targets, vis_coords, strand, gene, isoform, gene_isoforms = parse_targets(
            args.targets, args.genome, use_db, db, pad_size, args.targetRegion, args.exons,
            args.targetUpstreamPromoter, args.targetDownstreamPromoter,
            config.path("TWOBIT_INDEX_DIR") if not config.isoforms else config.path("ISOFORMS_INDEX_DIR"),
            args.outputDir, args.consensusUnion, args.jsonVisualize, args.guideSize)

        sequences, fasta_sequence = coord_to_fasta(
            targets, candidate_fasta_file, args.outputDir, args.guideSize, eval_sequence, args.nonOverlapping,
            config.path("TWOBIT_INDEX_DIR") if not config.isoforms else config.path("ISOFORMS_INDEX_DIR"),
            args.genome, strand, DOWNSTREAM_NUC)

    # Converts genomic coordinates to fasta file of all possible k-mers
    if len(sequences) == 0:
        sys.stderr.write("No target sites\n")
        sys.exit()

    # Run bowtie and get results
    bowtie_results_file = run_bowtie(len(args.PAM), args.uniqueMethod_Cong, candidate_fasta_file, args.outputDir,
                                   int(args.maxOffTargets),
                                   config.path("ISOFORMS_INDEX_DIR") if config.isoforms else config.path(
                                       "BOWTIE_INDEX_DIR"),
                                   args.genome, int(args.maxMismatches))

    results = parse_bowtie(guide_class, bowtie_results_file, True, args.scoreGC, args.scoreSelfComp,
                           args.backbone, args.replace5P, args.maxOffTargets, count_MM, args.PAM,
                           args.MODE != ProgramMode.TALENS,
                           args.scoringMethod, args.genome, gene, isoform,
                           gene_isoforms)  # TALENS: MAKE_PAIRS + CLUSTER

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

    info = scoring.ScoringInfo(args.genome, args.PAM, strand, sort_output, cluster_info, args.outputDir,
                               args.repairPredictions is not None, args.repairPredictions, args.isoforms,
                               vis_coords, args.fasta, args.rm1perfOff, args.MODE, scoring_method)
    sorted_output, cluster = scoring.score_guides(results, info)

    # Write individual results to file
    list_of_clusters = write_individual_results(args.outputDir, args.maxOffTargets, sorted_output,
                                              args.guideSize, args.MODE, cluster,
                                              args.limitPrintResults, args.offtargetsTable)

    if args.makePrimers:
        if args.fasta:
            make_primers_fasta(sorted_output, args.outputDir, args.primerFlanks,
                               args.displaySeqFlanks, args.genome, args.limitPrintResults,
                               config.path("BOWTIE_INDEX_DIR"), fasta_sequence,
                               args.primer3options, args.guidePadding, args.enzymeCo,
                               args.minResSiteLen, "sequence", args.maxOffTargets)
        else:
            make_primers_genome(sorted_output, args.outputDir, args.primerFlanks,
                                args.displaySeqFlanks, args.genome, args.limitPrintResults,
                                config.path("BOWTIE_INDEX_DIR"),
                                config.path("TWOBIT_INDEX_DIR") if not config.isoforms else config.path(
                                    "ISOFORMS_INDEX_DIR"), args.primer3options,
                                args.guidePadding, args.enzymeCo, args.minResSiteLen, strand,
                                args.targets, args.maxOffTargets)

    # Print results
    print_scores(sorted_output, args.MODE, args.scoringMethod, args.isoforms)

    result_coords = generate_result_coordinates(sorted_output,
                                               list_of_clusters,
                                               sort_output,
                                               args.MODE,
                                               args.isoforms)

    # Print gene annotation files
    # FASTA file
    gene_file = open('%s/gene_file.fa' % args.outputDir, 'w')
    gene_file.write(">%s\n" % args.targets)
    gene_file.write(fasta_sequence)
    gene_file.close()

    # Visualize with json
    if args.jsonVisualize:
        # new function
        visualize_with_json(args, vis_coords, sequences, strand, result_coords, fasta_sequence, targets)

    # remove .sam files as they take up wayyy to much space
    for fl in os.listdir(args.outputDir):
        if fl == "primer_results.sam" or fl == "output.sam":
            os.remove(os.path.join(args.outputDir, fl))


if __name__ == '__main__':
    main()
