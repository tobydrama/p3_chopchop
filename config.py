import json
import logging
import os
import sys

# File path
_file_path = sys.path[0]


# Main config file
_config = {}

if os.path.isfile(_file_path + "/config_local.json"):
    with open(_file_path + "/config_local.json") as f:
        _config = json.load(f)
elif os.path.isfile(_file_path + "/config.json"):
    with open(_file_path + "/config.json") as f:
        _config = json.load(f)
else:
    logging.critical("Could not find 'config_local.json' or 'config.json'; exiting.")
    sys.exit(1)


# Lib / directory paths
_config_paths = _config['PATH']

for _location in ["PRIMER3", "BOWTIE", "TWOBITTOFA", "TWOBIT_INDEX_DIR", "BOWTIE_INDEX_DIR", "ISOFORMS_INDEX_DIR",
                  "ISOFORMS_MT_DIR", "GENE_TABLE_INDEX_DIR"]:
    if _location not in _config_paths:
        logging.warning("%s is missing from config/paths!")


# Bowtie multithreading config
_config_threads = _config['THREADS']
if type(_config_threads) != int:
    logging.error("Config/THREADS misconfigured, setting to '1'.")
    _config_threads = 1

# Scoring config TODO add this to config.json
_config_scoring = _config['SCORING']


def path_to(loc: str) -> str:
    return _config_paths[loc]


@property
def file_path() -> str:
    return _file_path


@property
def thread_count() -> int:
    return _config_threads

_isoforms = False

@property
def use_isoforms() -> bool:
    return _isoforms