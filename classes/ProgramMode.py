import argparse
from enum import Enum


class ProgramModeAction(argparse.Action):
    """
    Action class to convert from integer value to enum in the command-line parsing step.
    """

    def __init__(self, option_strings, dest=None, nargs=0, default=None, required=False, type=None, choices=None,
                 metavar=None, help=None):
        super(ProgramModeAction, self).__init__(option_strings, dest=dest, nargs=nargs, default=ProgramMode(default),
                                                required=required, type=ProgramMode, choices=choices, metavar=metavar,
                                                help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        if values in self.choices:
            setattr(namespace, self.dest, ProgramMode(values[0]))
        print("YEP")


class ProgramMode(Enum):
    CRISPR = 1
    TALENS = 2
    CPF1 = 3
    NICKASE = 4
