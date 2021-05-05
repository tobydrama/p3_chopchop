import codecs
import pickle
import subprocess


def run_coefficient_score(stranded_guide_seq: str) -> float:
    """
    **Docker wrapper for CRISPRoff.**

    Runs the `CRISPRoff_specificity.CRISPRoff_score()` function with supplied input and returns a float.

    WARNING: this function may raise an exception if the main python script does not have either 'dockers' group
    permissions, or root (sudo) permissions; dockers requires either group or root permissions to run.

    :param stranded_guide_seq: The input sequence for CRISPRoff_score()
    :return: A float value. Coefficient score maybe?
    :raises subprocess.CalledProcessError: if dockers returns a non-zero exit value.
    """
    command = ['docker', 'run', 'chopchop_crisproff', stranded_guide_seq]

    """Error handling is done by subprocess.run flag 'check=True',
    if the dockers container returns a non-zero return value,
    the subprocess will raise a 'CalledProcessError' with captured STDERR"""
    crisproff = subprocess.run(command, capture_output=True, text=True, check=True)

    return pickle.loads(codecs.decode(crisproff.stdout.encode(), 'base64'), encoding='latin1')
