#!/usr/bin/env python2.7
import argparse
# Base64 encoding for guides pickle
import codecs
import pickle
# Receive guides through STDIN
import sys
import warnings

import numpy
import pandas

import Cas9Emulation
import featurization as feat


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--scoringMethod", choices=["DOENCH_2016", "ALL"], required=True,
                        help="The scoring method being used in the main CHOPCHOP script.")
    parser.add_argument("-c", "--coefficientScore", required=True, type=int,
                        help="The coefficient score multiplier")
    return parser.parse_args()


def recv_tuples():
    """
    Receives a Base64 encoded pickled list of tuples containing arguments for Cas9Emulation objects.

    Source: https://stackoverflow.com/q/30469575

    :return a list of tuples
    """
    pickled = sys.stdin.read()

    tuples = pickle.loads(codecs.decode(pickled.encode(), 'base64'))

    return tuples


def tuple_to_cas9(t):
    """
    Takes a tuple containing all fields of the Cas9Emulation object and returns a Cas9Emulation object.

    :param t: The tuple representation of a Cas9Emulation object.
    :return: The corresponding Cas9Emulation object.
    """
    return Cas9Emulation.Cas9Emulation(*t)


def cas9_to_reduced_tuple(guide):
    """
    Returns a reduced tuple containing only the key & modiied values of the Cas9Emulation objects after running
    Doench 2016.

    :param guide: The guide to reduce into a tuple.
    :return: A tuple containing the inout guide's key, score & "DOENCH_"2016" coefficient score.
    """
    return guide.key, guide.score, guide.CoefficientsScore["DOENCH_2016"]


def concatenate_feature_sets(feature_sets):
    '''
    Given a dictionary of sets of features, each in a Pandas.DataFrame,
    concatenate them together to form one big np.array, and get the dimension
    of each set
    Returns: inputs, dim
    Source: Doench 2016
    '''
    assert feature_sets != {}, "no feature sets present"
    F = feature_sets[feature_sets.keys()[0]].shape[0]
    for fset in feature_sets.keys():
        F2 = feature_sets[fset].shape[0]
        assert F == F2, "not same # individuals for features %s and %s" % (feature_sets.keys()[0], fset)

    N = feature_sets[feature_sets.keys()[0]].shape[0]
    inputs = numpy.zeros((N, 0))
    feature_names = []
    dim = {}
    dimsum = 0
    for fset in feature_sets.keys():
        inputs_set = feature_sets[fset].values
        dim[fset] = inputs_set.shape[1]
        dimsum += dim[fset]
        inputs = numpy.hstack((inputs, inputs_set))
        feature_names.extend(feature_sets[fset].columns.tolist())

    return inputs, dim, dimsum, feature_names


def doench_2016_score(args, results):
    """
    Minimally modified Doench 2016 algorithm from CHOPCHOP 2.7 lines 3120-3168

    :param args: argparse.ArgumentParser instance containing field `scoringMethod`, used to set scores.
    :param results: The guides to evaluate.
    :return: The evaluated guides.
    """
    try:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("ignore")
            with open('models/Doench_2016_18.01_model_nopos.pickle', 'rb') as f:  # TODO check that file path is correct
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
                    guide.score -= (guide.CoefficientsScore["DOENCH_2016"] / 100) * args.coefficientScore

        return results

    except:
        pass


def main():
    args = parse_args()

    guides = []
    for t in recv_tuples():
        guides.append(tuple_to_cas9(t))

    scored_guides = doench_2016_score(args, guides)

    # If Doench 2016 did not run, return exit code 1
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
