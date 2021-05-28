import codecs
import pickle
import subprocess

import config
from classes.Cas9 import Cas9


def convert_cas9_to_tuple(key: int, guide: Cas9) -> (int, str, str, str, str, float, float):
    return (key,
            guide.downstream_5_prim,
            guide.downstream_3_prim,
            guide.stranded_guide_seq,
            guide.pam,
            guide.score,
            guide.coefficients_score)


def run_chari_2015(guides: [Cas9], info) -> [Cas9]:
    """
    Runs chopchop_chari_2015 docker image using the supplied guides & scoring method.

    :param scoring_method: The scoring method to use. Accepted values are "chari_2015" & "ALL".
    :param guides: A list of Cas9 objects to score.
    :return: Returns a list of Cas9 scored objects.
    """
    keyed_tuples = []
    for key, guide in enumerate(guides):
        keyed_tuples.append(convert_cas9_to_tuple(key, guide))

    encoded = codecs.encode(pickle.dumps(keyed_tuples, protocol=2), 'base64').decode()

    command = ['docker', 'run', '-i', 'chopchop_chari_2015', '-c', str(config.score('COEFFICIENTS')), '-p',
               str(info.pam), '-g', str(info.genome), '-s', str(info.scoring_method.name)]
    chari_2015 = subprocess.run(command, capture_output=True, text=True, input=encoded)

    # encoding='latin1' is for backwards compatibility.
    results = pickle.loads(codecs.decode(chari_2015.stdout.encode(), 'base64'), encoding='latin1')

    # TODO currently we loop through key, guide pairs to set results, as we did for the keyed_tuples. We might want to
    #  save keyed_tuple & guide relationships to be sure the different guides get the correct scores.
    for key, guide in enumerate(guides):
        for t in results:
            if t[0] == key:
                _, guide.score, guide.coefficients_score = t

    return guides
