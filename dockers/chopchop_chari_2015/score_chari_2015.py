#!/usr/bin/env python2.7

import argparse
# Base64 encoding for guides pickle
import codecs
import pickle
# Receive guides through STDIN
import sys
import Cas9Emulation
from collections import defaultdict
from subprocess import Popen
import scipy.stats as ss


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--coefficientScore", required=True, type=int,
                        help="The coefficient score multiplier")
    parser.add_argument("-p", "--pam", required=True,
                        help="The pam")
    parser.add_argument("-g", "--genome", required=True,
                        help="The genome")
    parser.add_argument("-s", "--scoring_method", required=True,
                        help="The scoring_method")
    return parser.parse_args()


def recv_tuples():
    """
    Receives a Base64 encoded pickled list of tuples containing arguments for cas9Emulation objects.

    Source: https://stackoverflow.com/q/30469575

    :return a list of tuples
    """
    pickled = sys.stdin.read()

    tuples = pickle.loads(codecs.decode(pickled.encode(), 'base64'))

    return tuples


def tuple_to_cas9(t):
    """
    Takes a tuple containing all fields of the cas9Emulation object and returns a cas9Emmulation object.

    :param t: The tuple representation of a cas9Emulation object.
    :return: The corresponding cas9Emulation object.
    """
    return Cas9Emulation.Cas9Emulation(*t)


def cas9_to_reduced_tuple(guide):
    """
    Returns a reduced tuple containing only the key & modiied values of the cas9Emulation objects after running
    chari 2015.

    :param guide: The guide to reduce into a tuple.
    :return: A tuple containing the inout guide's key, score & "chari 2015" coefficient score.
    """
    return guide.key, guide.score, guide.CoefficientsScore


def run_svm_classify(svm_input_file, svm_output_file, PAM, genome):  # Only one use in main
    f_p = sys.path[0]
    """ Calculate score from SVM model as in Chari 2015 20-NGG or 20-NNAGAAW, only for hg19 and mm10"""

    model = f_p + '/models/293T_HiSeq_SP_Nuclease_100_SVM_Model.txt'
    dist = f_p + '/models/Hg19_RefFlat_Genes_75bp_NoUTRs_SPSites_SVMOutput.txt'

    if PAM == 'NGG' and genome == 'mm10':
        model = f_p + '/models/293T_HiSeq_SP_Nuclease_100_SVM_Model.txt'
        dist = f_p + '/models/Mm10_RefFlat_Genes_75bp_NoUTRs_SPSites_SVMOutput.txt'
    elif PAM == 'NNAGAAW' and genome == 'hg19':
        model = f_p + '/models/293T_HiSeq_ST1_Nuclease_100_V2_SVM_Model.txt'
        dist = f_p + '/models/Hg19_RefFlat_Genes_75bp_NoUTRs_ST1Sites_SVMOutput.txt'
    elif PAM == 'NNAGAAW' and genome == 'mm10':
        model = f_p + '/models/293T_HiSeq_ST1_Nuclease_100_V2_SVM_Model.txt'
        dist = f_p + '/models/Mm10_RefFlat_Genes_75bp_NoUTRs_ST1Sites_SVMOutput.txt'

    prog = Popen("%s/svm_light/svm_classify -v 0 %s %s %s" % (f_p, svm_input_file, model, svm_output_file), shell=True)
    prog.communicate()

    svm_all = open(dist, 'r')
    svm_this = open(svm_output_file, 'r')

    # first through go all scores and get the max and min
    all_data = []
    for line in svm_all:
        line = line.rstrip('\r\n')
        all_data.append(float(line))
    svm_all.close()

    score_array = []
    for line in svm_this:
        line = line.rstrip('\r\n')
        score_array.append(float(line))

    return [ss.percentileofscore(all_data, i) for i in score_array]



def score_chari_2015(guides, args) :
    svm_input_file = "chari_score.SVMInput.txt"
    svm_output_file = "chari_score.SVMOutput.txt"

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
    new_scores = run_svm_classify(svm_input_file, svm_output_file, args.pam, args.genome)

    for i, guide in enumerate(guides):
        guide.CoefficientsScore["CHARI_2015"] = new_scores[i]
        if args.scoring_method == "CHARI_2015": #ScoringMethod.CHARI_2015:
            guide.score -= (guide.CoefficientsScore["CHARI_2015"] / 100) * args.coefficientScore



    return guides

def main():
    args = parse_args()
    guides = []
    for t in recv_tuples():
        guides.append(tuple_to_cas9(t))

    scored_guides = score_chari_2015(guides, args)

    # If chari 2015 did not run, return exit code 1
    if not scored_guides:
        exit(1)

    tuples = []
    for guide in scored_guides:
        tuples.append(cas9_to_reduced_tuple(guide))

    # Encode & print the pickled tuples to STDOUT for the main script to catch.
    print codecs.encode(pickle.dumps(tuples), 'base64').decode()

    for t in tuples:
        sys.stderr.write("KEY: %d\tSCORE: %0.20f\tCOEFFICIENT %0.20f\n" % (t[0], t[-2], t[-1]))


if __name__ == '__main__':
    main()
