#!/usr/bin/env python3

import os
import sys
import json
import string
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("-o", "--outputDir", default="./", metavar="OUTPUT_DIR", help="The output directory. Default is the current directory.")

args = parser.parse_args()


