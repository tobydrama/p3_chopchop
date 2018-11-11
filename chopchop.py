#!/usr/bin/env python3

import os
import sys
import json
import string
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("-o", "--outputDir", default="./", metavar="OUTPUT_DIR", help="The output directory. Default is the current directory.")

args = parser.parse_args()

data = {'a list': [1, 42, 3.141, 1337, 'help', u'€'],
        'a string': 'bla',
        'another dict': {'foo': 'bar',
                         'key': 'value',
                         'the answer': 42}}

with open(os.path.join(args.outputDir + 'results.json'), 'w') as f:
  json.dump(data, f, ensure_ascii=False)


