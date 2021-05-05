#!/usr/bin/env python2.7

import argparse
import codecs
import pickle

from CRISPRoff.CRISPRoff_specificity import CRISPRoff_score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('guide_seq', help="The stranded_guide_seq to use for the CRISRoff score function.")
    args = parser.parse_args()

    if args.guide_seq is None:
        exit(1)

    result = CRISPRoff_score(args.guide_seq)

    # Print result to stdout so that it can be caught by CHOPCHOP
    # TODO make sure floating point precision is the same as in CHOPCHOP 2.7
    print codecs.encode(pickle.dumps(result).decode(), 'base64').decode()


if __name__ == '__main__':
    main()
