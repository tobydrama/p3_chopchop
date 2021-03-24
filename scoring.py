import logging
import os
import subprocess
import sys
import warnings
from collections import defaultdict
from enum import Enum
from typing import Callable, List, Tuple, Dict, Union

import numpy
import scipy.stats as ss

import Vars
from classes.Cas9 import Cas9
from classes.Guide import Guide
from classes.ProgramMode import ProgramMode
import functions.Main_Functions as mainFunctions
import functions.TALEN_Specific_Functions as talenFunctions
import dockers.doench_2016_wrapper as doench_2016


class ScoringMethod(Enum):
    ALL = 0,
    XU_2015 = 1,
    DOENCH_2014 = 2,
    DOENCH_2016 = 3,
    MORENO_MATEOS_2015 = 4,
    CHARI_2015 = 5,
    G_20 = 6,
    KIM_2018 = 7,
    ALKAN_2018 = 8,
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
            if guide.offTargetsMM[0] > 0:
                # TODO do something about this constant
                guide.score -= Vars.SINGLE_OFFTARGET_SCORE[0]

    if info.use_isoforms:
        guides = score_isoforms(guides, info)
    else:
        if (info.scoring_method == ScoringMethod.CHARI_2015 or info.scoring_method == ScoringMethod.ALL) \
                and info.pam in ["NGG", "NNAGAAW"] and info.genome in ["hg19", "mm10"]:
            # TODO Type hinting is confused here.
            guides = score_chari_2015(guides, info)

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


def score_isoforms(guides: List[Guide], info: ScoringInfo) -> List[Guide]:
    logging.info("Scoring isoforms.")

    for guide in guides:
        # Calculate base pair probabilities of folding
        if guide.isoform in ["union", "intersection"]:
            bpp = []

            for tx_id in guide.gene_isoforms:
                tx_start, tx_end = mainFunctions.tx_relative_coordinates(info.vis_coords, tx_id, guide.start, guide.end)

                if tx_start != -1:
                    bpp.append(mainFunctions.rna_folding_metric(info.genome, guide.isoform, guide.start, guide.end))

            guide.meanBPP = 100 if len(bpp) == 0 else max(bpp)
        else:
            if not info.use_fasta:
                tx_start, tx_end = mainFunctions.tx_relative_coordinates(info.vis_coords, guide.isoform,
                                                                         guide.start, guide.end)

                guide.meanBPP = mainFunctions.rna_folding_metric(info.genome, guide.isoform, tx_start, tx_end)

        guide.score += guide.meanBPP / 100 * Vars.SCORE["COEFFICIENTS"]

        if guide.isoform in guide.gene_isoforms:
            guide.gene_isoforms.remove(guide.isoform)

        if guide.isoform in guide.offTargetsIso[0]:
            guide.offTargetsIso[0].remove(guide.isoform)

        guide.constitutive = int(guide.gene_isoforms == guide.offTargetsIso[0])

    return guides


def score_chari_2015(guides: List[Cas9], info: ScoringInfo) -> List[Guide]:
    logging.info("Running Chari 2015 scoring method.")

    try:
        svm_input_file = "%s/chari_score.SVMInput.txt" % info.output_dir
        svm_output_file = "%s/chari_score.SVMOutput.txt" % info.output_dir

        encoding = defaultdict(str)
        encoding['A'] = '0001'
        encoding['C'] = '0010'
        encoding['T'] = '0100'
        encoding['G'] = '1000'

        svm_file = open(svm_input_file, 'w')

        for guide in guides:
            seq = guide.downstream5prim + guide.strandedGuideSeq[:-len(guide.PAM)]
            pam = guide.strandedGuideSeq[-len(guide.PAM):]
            sequence = (seq[-20:] + pam).upper()

            # end index
            if len(sequence) == 27:
                end_index = 22
            else:
                end_index = 21

            x = 0
            tw = '-1'
            while x < end_index:
                y = 0
                while y < 4:
                    tw = tw + ' ' + str(x + 1) + str(y + 1) + ':' + encoding[sequence[x]][y]
                    y += 1
                x += 1
            svm_file.write(tw + '\n')
        svm_file.close()
        new_scores = mainFunctions.scoreChari_2015(svm_input_file, svm_output_file, info.pam, info.genome)

        for i, guide in enumerate(guides):
            guide.CoefficientsScore["CHARI_2015"] = new_scores[i]
            if info.scoring_method == ScoringMethod.CHARI_2015:
                guide.score -= (guide.CoefficientsScore["CHARI_2015"] / 100) * Vars.SCORE['COEFFICIENTS']

        logging.debug("Finished running Chari 2015.")

    except:  # TODO what exceptions are caught here?
        logging.warning("Chari 2015 failed!")
        pass

    return guides


def score_zhang_2019(guides: List[Cas9], info: ScoringInfo) -> List[Guide]:
    logging.info("Running Zhang 2019 scoring method.")

    try:
        zhang_input_file = '%s/zhang_score.txt' % info.output_dir
        zhang_file = open(zhang_input_file, 'w')
        for guide in guides:
            zhang_file.write(guide.downstream5prim[-4:] + guide.strandedGuideSeq + guide.downstream3prim[:3] + '\n')
        zhang_file.close()

        # TODO kind of hacky file path stuff here
        prog = subprocess.run(["%s/uCRISPR/uCRISPR" % Vars.f_p, "-on", Vars.f_p + '/' + zhang_input_file],
                              capture_output=True, check=True)

        # Convert from bytes, split on newline, remove header
        output = prog.stdout.decode().splitlines()[1:]

        # distribution calculated on 100k random guides
        output = [ss.norm.cdf(float(x.split()[1]), loc=11.92658, scale=0.2803797) for x in output]
        for i, guide in enumerate(guides):
            guide.CoefficientsScore["ZHANG_2019"] = output[i] * 100

            if info.scoring_method == ScoringMethod.ZHANG_2019:
                guide.score -= (guide.CoefficientsScore["ZHANG_2019"] / 100) * Vars.SCORE['COEFFICIENTS']

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

            os.environ['KERAS_BACKEND'] = 'theano'

            # Keras prints welcome message to stderr
            stderr = sys.stderr
            sys.stderr = open(os.devnull)

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
            seq_deep_cpf1.load_weights(Vars.f_p + '/models/Seq_deepCpf1_weights.h5')

            # process data
            data_n = len(guides)
            one_hot = numpy.zeros((data_n, 34, 4), dtype=int)

            for l in range(0, data_n):
                prim5 = guides[l].downstream5prim[-4:]
                if len(prim5) < 4:  # cover weird genomic locations
                    prim5 = "N" * (4 - len(prim5)) + prim5
                guide_seq = guides[l].strandedGuideSeq
                prim3 = guides[l].downstream3prim[:6]
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

        for i, guide in enumerate(guides):
            guide.CoefficientsScore = seq_deep_cpf1_score[i][0]
            guide.score -= (guide.CoefficientsScore / 100) * Vars.SCORE["COEFFICIENTS"]

        logging.debug("Finished running Kim 2018.")

    except:  # TODO what exceptions are caught here?
        logging.warning("Kim 2018 failed!")
        pass

    return guides


def score_doench_2016(guides: List[Guide], scoring_method: ScoringMethod) -> List[Guide]:
    logging.info("Running Doench 2016 scoring method.")

    return doench_2016.run_doench_2016(scoring_method.name, guides)


def run_repair_predictions(guides: List[Guide], repair_predictions: str) -> List[Guide]:
    logging.info("Running inDelphi repair predictions on %d guides" % len(guides))

    # TODO 'Vars.f_p' isn't very pretty, change it
    sys.path.append(Vars.f_p + "/models/inDelphi-model/")

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("ignore")

        import inDelphi
        inDelphi.init_model(celltype=repair_predictions)

        for i, guide in enumerate(guides):
            try:
                left_seq = guide.downstream5prim + guide.strandedGuideSeq[:-(len(guide.PAM) + 3)]
                left_seq = left_seq[-60:]

                right_seq = guide.strandedGuideSeq[-(len(guide.PAM) + 3):] + guide.downstream3prim
                right_seq = right_seq[:60]

                seq = left_seq + right_seq
                cut_site = len(left_seq)

                pred_df, stats = inDelphi.predict(seq, cut_site)
                pred_df = pred_df.sort_values(pred_df.columns[4], ascending=False)

                guide.repProfile = pred_df
                guide.repStats = stats
            except:  # TODO what exceptions are caught here?
                pass

    logging.debug("Finished running inDelphi repair predictions.")

    return guides


def get_cluster_pairs(guides: List[Guide], info: ScoringInfo, program_mode: ProgramMode) -> Tuple[int, List[Guide]]:
    logging.info("Acquiring cluster pairs from guides.")

    pairs = []

    cluster_info = info.cluster_info
    if program_mode == ProgramMode.TALENS:
        pairs = talenFunctions.pairTalens(guides,
                                          cluster_info.sequences, cluster_info.guide_size,
                                          int(cluster_info.min_distance), int(cluster_info.max_distance),
                                          cluster_info.enzyme_company, cluster_info.max_off_targets,
                                          cluster_info.g_rvd, cluster_info.min_res_site_len)
    elif program_mode == ProgramMode.NICKASE:
        pairs = talenFunctions.pairCas9(guides,
                                        cluster_info.sequences, cluster_info.guide_size,
                                        int(cluster_info.min_distance), int(cluster_info.max_distance),
                                        cluster_info.enzyme_company, cluster_info.max_off_targets,
                                        cluster_info.min_res_site_len, cluster_info.off_target_max_distance)

    if not len(pairs):
        # TODO better error handling (error logging, exceptions, etc)
        sys.stderr.write("No " + ("TALEN" if program_mode == ProgramMode.TALENS else "Cas9 NICKASE")
                         + " pairs could be generated for this region.\n")
        sys.exit(Vars.EXIT['GENE_ERROR'])

    if info.rm1_perf_off and info.use_fasta:
        for pair in pairs:
            if pair.diffStrandOffTarget > 0:
                pair.score -= Vars.SCORE["OFFTARGET_PAIR_DIFF_STRAND"]
            if program_mode == ProgramMode.TALENS and pair.sameStrandOffTarget > 0:
                pair.score -= Vars.SCORE["OFFTARGET_PAIR_SAME_STRAND"]

    cluster, guides = talenFunctions.clusterPairs(pairs)

    logging.debug("Got %d clusters." % cluster)

    return cluster, guides
