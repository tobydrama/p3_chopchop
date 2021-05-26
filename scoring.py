import logging
import math
import os
import subprocess
import sys
import warnings
from enum import Enum
from typing import Callable, List, Tuple, Dict, Union

import numpy
import pandas
import scipy.stats as ss

import config
import constants
from classes.Cas9 import Cas9
from classes.Guide import Guide
from classes.ProgramMode import ProgramMode
from dockers.chari_2015_wrapper import run_chari_2015
from dockers.kim_2018_wrapper import run_kim_2018
from dockers.doench_2016_wrapper import run_doench_2016
from dockers.inDelphi_wrapper import run_repair_prediction
from functions.TALEN_specific_functions import pair_talens, pair_cas9, cluster_pairs


class ScoringMethod(Enum):
    ALL = 0
    XU_2015 = 1
    DOENCH_2014 = 2
    DOENCH_2016 = 3
    MORENO_MATEOS_2015 = 4
    CHARI_2015 = 5
    G_20 = 6
    KIM_2018 = 7
    ALKAN_2018 = 8
    ZHANG_2019 = 9


class ClusterInfo:
    def __init__(self,
                 sequences: Dict[str, str],
                 guide_size: int,
                 min_distance: int,
                 max_distance: int,
                 enzyme_company: str,
                 max_off_targets: int,
                 g_rvd: str,
                 min_res_site_len: int,
                 off_target_max_distance: int):
        self.sequences = sequences
        self.guide_size = guide_size
        self.min_distance = min_distance
        self.max_distance = max_distance
        self.enzyme_company = enzyme_company
        self.max_off_targets = max_off_targets
        self.g_rvd = g_rvd
        self.min_res_site_len = min_res_site_len
        self.off_target_max_distance = off_target_max_distance


class ScoringInfo:
    def __init__(self,
                 genome: str,
                 pam: str,
                 strand: str,
                 sort_fn: Callable[[List[Guide]], List[Guide]],
                 cluster_info: ClusterInfo,
                 output_dir: str,
                 use_repair_predictions: bool = False,
                 repair_predictions: str = None,
                 use_isoforms: bool = False,
                 vis_coords: List[Dict[str, Union[List[List[Union[int, str]]], List]]] = None,
                 use_fasta: bool = False,
                 rm1_perf_off: bool = False,
                 program_mode: ProgramMode = ProgramMode.CRISPR,
                 scoring_method: ScoringMethod = ScoringMethod.G_20):
        self.genome = genome
        self.pam = pam
        self.strand = strand
        self.cluster_info = cluster_info
        self.output_dir = output_dir
        self.sort_fn = sort_fn
        self.use_repair_predictions = use_repair_predictions
        self.repair_predictions = repair_predictions
        self.use_isoforms = use_isoforms
        self.vis_coords = vis_coords
        self.use_fasta = use_fasta
        self.rm1_perf_off = rm1_perf_off
        self.program_mode = program_mode
        self.scoring_method = scoring_method


def score_guides(guides: List[Guide], info: ScoringInfo) -> Tuple[List[Guide], int]:
    logging.info("Scoring %d guides with scoring method %s." % (len(guides), info.scoring_method.name))

    if info.rm1_perf_off and info.use_fasta:
        for guide in guides:
            if guide.off_targets_mm[0] > 0:
                # TODO do something about this constant
                guide.score -= constants.SINGLE_OFFTARGET_SCORE[0]

    if len(guides) > 0 and type(guides[0]) == Cas9:
        guides = score_cas9(guides)

    if info.scoring_method == ScoringMethod.ALKAN_2018 or info.scoring_method == ScoringMethod.ALL:
        guides = score_alkan_2018(guides)

    if info.use_isoforms:
        guides = score_isoforms(guides, info)
    else:
        if (info.scoring_method == ScoringMethod.CHARI_2015 or info.scoring_method == ScoringMethod.ALL) \
                and info.pam in ["NGG", "NNAGAAW"]:
            if info.genome in ["hg19", "mm10"]:
                # TODO Type hinting is confused here.
                guides = score_chari_2015(guides, info)
            else:
                logging.warning(f"Scoring method CHARI 2015 is incompatible with '{info.genome}', skipping. CHARI 2015 "
                                "compatible genomes are 'hg19' & 'mm10'.")

        if (info.scoring_method == ScoringMethod.ZHANG_2019 or info.scoring_method == ScoringMethod.ALL) \
                and info.pam == "NGG":
            guides = score_zhang_2019(guides, info)

        if (info.scoring_method == ScoringMethod.KIM_2018 or info.scoring_method == ScoringMethod.ALL) \
                and info.pam in "TTTN" \
                and info.program_mode == ProgramMode.CPF1:
            guides = score_kim_2018(guides)

        if info.program_mode == ProgramMode.CRISPR:
            if info.scoring_method == ScoringMethod.DOENCH_2016 or info.scoring_method == ScoringMethod.ALL:
                guides = score_doench_2016(guides, info.scoring_method)

            if info.use_repair_predictions:
                guides = run_repair_predictions(guides, info.repair_predictions)

    cluster = 0
    if info.program_mode in [ProgramMode.TALENS, ProgramMode.NICKASE]:
        cluster, guides = get_cluster_pairs(guides, info, info.program_mode)

    if info.strand == '-' and not info.use_isoforms:
        guides.reverse()

    sorted_guides = info.sort_fn(guides)

    return sorted_guides, cluster


def tx_relative_coordinates(vis_coords, tx_id, start, end):
    tx_start, tx_end = -1, -1
    exons = [e["exons"] for e in vis_coords if e["name"] == tx_id][0]
    e_id = -1
    for i, e in enumerate(exons):
        if e[1] <= (start - 1) and e[2] >= (end - 1):
            e_id = i
            break

    if e_id != -1:
        for i in range(0, e_id) if exons[0][5] == "+" else range(e_id + 1, len(exons)):
            tx_start += exons[i][2] - exons[i][1]

        tx_start += (exons[e_id][1] - start - 1) if exons[0][5] == "+" else (exons[e_id][2] - end - 1)
        tx_end = tx_start + end - start

    return tx_start, tx_end


def rna_folding_metric(specie, tx_id, tx_start, tx_end):
    mean_bpp = 0
    file_path = config.path("ISOFORMS_MT_DIR") + "/" + specie + "/" + tx_id + ".mt"
    if os.path.isfile(file_path):
        mt = pandas.read_csv(file_path, sep="\t", header=None, skiprows=tx_start, nrows=tx_end - tx_start)
        mean_bpp = numpy.mean(mt[1].tolist())

    return mean_bpp


def score_isoforms(guides: List[Guide], info: ScoringInfo) -> List[Guide]:
    logging.info("Scoring isoforms.")

    for guide in guides:
        # Calculate base pair probabilities of folding
        if guide.isoform in ["union", "intersection"]:
            bpp = []

            for tx_id in guide.gene_isoforms:
                tx_start, tx_end = tx_relative_coordinates(info.vis_coords, tx_id, guide.start, guide.end)

                if tx_start != -1:
                    bpp.append(rna_folding_metric(info.genome, tx_id, tx_start, tx_end))

            guide.mean_bpp = 100 if len(bpp) == 0 else max(bpp)
        else:
            if not info.use_fasta:
                tx_start, tx_end = tx_relative_coordinates(info.vis_coords, guide.isoform,
                                                           guide.start, guide.end)

                guide.mean_bpp = rna_folding_metric(info.genome, guide.isoform, tx_start, tx_end)

        guide.score += guide.mean_bpp / 100 * config.score("COEFFICIENTS")

        if guide.isoform in guide.gene_isoforms:
            guide.gene_isoforms.remove(guide.isoform)

        if guide.isoform in guide.off_targets_iso[0]:
            guide.off_targets_iso[0].remove(guide.isoform)

        guide.constitutive = int(guide.gene_isoforms == guide.off_targets_iso[0])

    return guides


def score_cas9(guides: List[Cas9]):
    logging.info("Scoring coefficients for Cas9 guides.")
    for guide in guides:
        if guide.scoring_method not in ["CHARI_2015", "DOENCH_2016", "ALKAN_2018", "ZHANG_2019", "ALL"]:
            guide.coefficients_score[guide.scoring_method] = score_grna(
                guide.downstream_5_prim + guide.stranded_guide_seq[:-len(guide.pam)],
                guide.stranded_guide_seq[-len(guide.pam):],
                guide.downstream_3_prim,
                constants.scores[guide.scoring_method]
            )
            guide.score -= guide.coefficients_score[guide.scoring_method] * config.score('COEFFICIENTS')
        elif guide.scoring_method == "ALL":
            for met in ["XU_2015", "DOENCH_2014", "MORENO_MATEOS_2015", "G_20"]:
                guide.coefficients_score[met] = score_grna(
                    guide.downstream_5_prim + guide.stranded_guide_seq[:-len(guide.pam)],
                    guide.stranded_guide_seq[-len(guide.pam):],
                    guide.downstream_3_prim,
                    constants.scores[met]
                )

    logging.debug("Finished scoring coefficients.")

    return guides


def score_alkan_2018(guides: List[Cas9]) -> List[Guide]:
    logging.info("Running Alkan 2018 scoring method.")

    from dockers.CRISPRoff_wrapper import run_coefficient_score

    for guide in guides:
        guide.coefficients_score['ALKAN_2018'] = run_coefficient_score(guide.stranded_guide_seq)
        guide.score -= guide.coefficients_score['ALKAN_2018'] * config.score('COEFFICIENTS')

    logging.debug("Finished running Alkan 2018.")

    return guides


def score_grna(seq, pam, tail, lookup):
    """ Calculate score from model coefficients. score is 0-1, higher is better """
    score = 0
    if "Intercept" in lookup:
        score = lookup["Intercept"]

    seq = seq[::-1]  # we calculate from PAM in a way: 321PAM123

    if "gc_low" in lookup:
        gc = seq[:20].count('G') + seq[:20].count('C')
        if gc < 10:
            score = score + (abs(gc - 10) * lookup["gc_low"])
        elif gc > 10:
            score = score + ((gc - 10) * lookup["gc_high"])

    for i in range(len(seq)):
        key = seq[i] + str(i + 1)
        if key in lookup:
            score += lookup[key]

        if i + 1 < len(seq):
            double_key = seq[i] + seq[i + 1] + str(i + 1)
            if double_key in lookup:
                score += lookup[double_key]

        if i == 0:
            double_key = pam[0] + seq[0] + str(0)
            if double_key in lookup:
                score += lookup[double_key]

    for i in range(len(pam)):
        key = 'PAM' + pam[i] + str(i + 1)
        if key in lookup:
            score += lookup[key]

        if i + 1 < len(pam):
            double_key = 'PAM' + pam[i] + pam[i + 1] + str(i + 1)
            if double_key in lookup:
                score += lookup[double_key]

    for i in range(len(tail)):
        key = str(i + 1) + tail[i]
        if key in lookup:
            score += lookup[key]

        if i + 1 < len(tail):
            double_key = str(i + 1) + tail[i] + tail[i + 1]
            if double_key in lookup:
                score += lookup[double_key]

    score = 1 / (1 + math.e ** -score)
    return score


def score_chari_2015(guides: List[Cas9], info: ScoringInfo) -> List[Guide]:
    logging.info("Running Chari 2015 scoring method.")

    return run_chari_2015(guides, info)


def score_zhang_2019(guides: List[Cas9], info: ScoringInfo) -> List[Guide]:
    logging.info("Running Zhang 2019 scoring method.")

    try:
        zhang_input_file = '%s/zhang_score.txt' % info.output_dir
        zhang_file = open(zhang_input_file, 'w')
        for guide in guides:
            zhang_file.write(guide.downstream_5_prim[-4:] + guide.stranded_guide_seq +
                             guide.downstream_3_prim[:3] + '\n')
        zhang_file.close()

        # TODO kind of hacky file path stuff here
        prog = subprocess.run(["%s/uCRISPR/uCRISPR" % config.file_path(), "-on", config.file_path() + '/' +
                               zhang_input_file], capture_output=True, check=True)

        # Convert from bytes, split on newline, remove header
        output = prog.stdout.decode().splitlines()[1:]

        # distribution calculated on 100k random guides
        output = [ss.norm.cdf(float(x.split()[1]), loc=11.92658, scale=0.2803797) for x in output]
        for i, guide in enumerate(guides):
            guide.coefficients_score["ZHANG_2019"] = output[i] * 100

            if info.scoring_method == ScoringMethod.ZHANG_2019:
                guide.score -= (guide.coefficients_score["ZHANG_2019"] / 100) * config.score('COEFFICIENTS')

        logging.debug("Finished running Zhang 2019.")

    except:  # TODO this is most likely a subprocess exception? Add proper catch
        logging.warning("Zhang 2019 failed: %s" % sys.exc_info()[0])  # ew
        pass

    return guides


def score_kim_2018(guides: List[Guide]) -> List[Guide]:
    logging.info("Running Kim 2018 scoring method.")

    try:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("ignore")

            guides = run_kim_2018(guides)

            logging.debug("Finished running Kim 2018.")

    except:  # TODO what exceptions are caught here?
        logging.warning("Kim 2018 failed!")
        pass
    return guides


def score_doench_2016(guides: List[Guide], scoring_method: ScoringMethod) -> List[Guide]:
    logging.info("Running Doench 2016 scoring method.")

    return run_doench_2016(scoring_method.name, guides)


def run_repair_predictions(guides: List[Guide], repair_predictions: str) -> List[Guide]:
    logging.info("Running inDelphi repair predictions on %d guides" % len(guides))

    return run_repair_prediction(repair_predictions, guides)


def get_cluster_pairs(guides: List[Guide], info: ScoringInfo, program_mode: ProgramMode) -> Tuple[int, List[Guide]]:
    logging.info("Acquiring cluster pairs from guides.")

    pairs = []

    cluster_info = info.cluster_info
    if program_mode == ProgramMode.TALENS:
        pairs = pair_talens(guides,
                            cluster_info.sequences, cluster_info.guide_size,
                            int(cluster_info.min_distance), int(cluster_info.max_distance),
                            cluster_info.enzyme_company, cluster_info.max_off_targets,
                            cluster_info.g_rvd, cluster_info.min_res_site_len)
    elif program_mode == ProgramMode.NICKASE:
        pairs = pair_cas9(guides,
                          cluster_info.sequences, cluster_info.guide_size,
                          int(cluster_info.min_distance), int(cluster_info.max_distance),
                          cluster_info.enzyme_company, cluster_info.max_off_targets,
                          cluster_info.min_res_site_len, cluster_info.off_target_max_distance)

    if not len(pairs):
        # TODO better error handling (error logging, exceptions, etc)
        sys.stderr.write("No " + ("TALEN" if program_mode == ProgramMode.TALENS else "Cas9 NICKASE")
                         + " pairs could be generated for this region.\n")
        sys.exit(constants.EXIT['GENE_ERROR'])

    if info.rm1_perf_off and info.use_fasta:
        for pair in pairs:
            if pair.diff_strand_off_target > 0:
                pair.score -= config.score("OFFTARGET_PAIR_DIFF_STRAND")
            if program_mode == ProgramMode.TALENS and pair.same_strand_off_target > 0:
                pair.score -= config.score("OFFTARGET_PAIR_SAME_STRAND")

    cluster, guides = cluster_pairs(pairs)

    logging.debug("Got %d clusters." % cluster)

    return cluster, guides
