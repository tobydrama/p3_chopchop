#!/usr/bin/env python2.7
import logging

import os

import argparse
# Base64 encoding for guides pickle
import codecs
import pickle
# Receive guides through STDIN
import sys
import numpy
import Cpf1Emulation


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--coefficientScore", required=True, type=int,
                        help="The coefficient score multiplier")

    return parser.parse_args()


def recv_tuples():
    """
    Receives a Base64 encoded pickled list of tuples containing arguments for Cpf1Emulation objects.

    Source: https://stackoverflow.com/q/30469575

    :return a list of tuples
    """
    pickled = sys.stdin.read()

    tuples = pickle.loads(codecs.decode(pickled.encode(), 'base64'))

    return tuples


def tuple_to_cpf1(t):
    """
    Takes a tuple containing all fields of the Cpf1Emulation object and returns a Cpf1Emmulation object.

    :param t: The tuple representation of a Cpf1Emulation object.
    :return: The corresponding Cpf1Emulation object.
    """
    return Cpf1Emulation.Cpf1Emulation(*t)


def cpf1_to_reduced_tuple(guide):
    """
    Returns a reduced tuple containing only the key & modiied values of the Cpf1Emulation objects after running
    KIM 2018.

    :param guide: The guide to reduce into a tuple.
    :return: A tuple containing the inout guide's key, score & "KIM_"2018" coefficient score.
    """
    return guide.key, guide.score, guide.CoefficientsScore


def score_kim_2018(args, guides):
    os.environ['KERAS_BACKEND'] = 'theano'

    from keras.models import Model
    from keras.layers import Input
    from keras.layers.merge import Multiply
    from keras.layers.core import Dense, Dropout, Activation, Flatten
    from keras.layers.convolutional import Convolution1D, AveragePooling1D

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
    seq_deep_cpf1.load_weights('./models/Seq_deepCpf1_weights.h5')

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
        guide.score -= (guide.CoefficientsScore / 100) * args.coefficientScore

    return guides


def main():
    args = parse_args()
    guides = []
    for t in recv_tuples():
        guides.append(tuple_to_cpf1(t))

    scored_guides = score_kim_2018(args, guides)

    # If KIM 2018 did not run, return exit code 1
    if not scored_guides:
        exit(1)

    tuples = []
    for guide in scored_guides:
        tuples.append(cpf1_to_reduced_tuple(guide))

    # Encode & print the pickled tuples to STDOUT for the main script to catch.
    print codecs.encode(pickle.dumps(tuples), 'base64').decode()

    for t in tuples:
        sys.stderr.write("KEY: %d\tSCORE: %0.20f\tCOEFFICIENT %0.20f\n" % (t[0], t[-2], t[-1]))


if __name__ == '__main__':
    main()
