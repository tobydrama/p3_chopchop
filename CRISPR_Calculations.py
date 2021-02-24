#!/usr/bin/env python2.7
import argparse

from CRISPRoff.CRISPRoff_specificity import CRISPRoff_score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('guide_seq', help="The stranded_guide_seq to use for the CRISRoff score function.")
    args = parser.parse_args()

    if args.guide_seq is None:
        exit(1)

    result = CRISPRoff_score(args.guide_seq)

    print(result)


if __name__ == '__main__':
    main()
