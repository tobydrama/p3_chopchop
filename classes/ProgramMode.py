from enum import Enum
import argparse


class ProgramModeAction(argparse.Action):
    """
    Action class to convert from integer value to enum in the command-line parsing step.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, ProgramMode(values))


class ProgramMode(Enum):
    CRISPR = 1
    TALENS = 2
    CPF1 = 3
    NICKASE = 4
