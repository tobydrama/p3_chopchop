import json
import logging
import os
import resource
import sys


class _Limits:
    def __init__(self):
        self._soft, self._hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        resource.setrlimit(resource.RLIMIT_NOFILE, (self._hard, self._hard))

    def soft(self) -> int:
        return self._soft

    def hard(self) -> int:
        return self._hard


# Resource limits
limits = _Limits()

# Absolute file path
_file_path = sys.path[0]


def file_path() -> str:
    """Returns the current absolute file path."""
    return _file_path


def get_config_from_file(*paths) -> dict:
    for p in paths:
        if os.path.isfile(p):
            with open(p) as f:
                return json.load(f)

    logging.critical(f"Could not find config file with supplied path(s): {paths}! Exiting.")
    sys.exit(1)


_config = get_config_from_file(file_path() + '/config_local.json', file_path() + '/config.json')
_paths = _config['PATH']

# Path safety check
for _location in ["PRIMER3", "BOWTIE", "TWOBITTOFA", "TWOBIT_INDEX_DIR", "BOWTIE_INDEX_DIR", "ISOFORMS_INDEX_DIR",
                  "ISOFORMS_MT_DIR", "GENE_TABLE_INDEX_DIR"]:
    if _location not in _paths:
        logging.warning("%s is missing from config/paths!")


def path(loc: str) -> str:
    """Returns path to directories defined in config.json"""
    return _paths[loc]


# Bowtie threads
_threads = _config['THREADS']

if type(_threads) != int:
    logging.error("Config/THREADS misconfigured, setting to '1'.")
    _threads = 1


def threads() -> int:
    """Returns config Bowtie thread count setting """
    return _threads


# Use isoforms
_isoforms = False


@property
def isoforms() -> bool:
    """Returns use-isoforms"""
    return _isoforms


_scoring = _config["SCORING"]


def score(key: str) -> any:
    return _scoring[key]
