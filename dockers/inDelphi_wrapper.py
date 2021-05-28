import codecs
import pickle
import subprocess

import config
from classes.Cas9 import Cas9


def convert_cas9_to_tuple(key: int, guide: Cas9) -> (int, str, str, str, str, any, any):
    return (key,
            guide.downstream_5_prim,
            guide.downstream_3_prim,
            guide.stranded_guide_seq,
            guide.pam,
            guide.repair_profile,
            guide.repair_stats)


def run_repair_prediction(repair_prediction: str, guides: [Cas9]) -> [Cas9]:
    """
    Runs chopchop_inDelphi docker image using the supplied guides & repair prediction mode.

    :param repair_prediction: The repair mode to use. Accepted values are 'mESC', 'U2OS', 'HEK293', 'HCT116' & 'K562'.
    :param guides: A list of Cas9 objects to score.
    :return: Returns a list of Cas9 scored objects.
    """
    keyed_tuples = []
    for key, guide in enumerate(guides):
        keyed_tuples.append(convert_cas9_to_tuple(key, guide))

    encoded = codecs.encode(pickle.dumps(keyed_tuples, protocol=2), 'base64').decode()

    decoded = pickle.loads(codecs.decode(encoded.encode(), 'base64'))

    command = ['docker', 'run', '-i', 'chopchop_indelphi', '-r', repair_prediction]

    repair_prediction = subprocess.run(command, capture_output=True, text=True, input=encoded)

    # encoding='latin1' is for backwards compatibility.
    results = pickle.loads(codecs.decode(repair_prediction.stdout.encode(), 'base64'), encoding='latin1')

    # TODO currently we loop through key, guide pairs to set results, as we did for the keyed_tuples. We might want to
    #  save keyed_tuple & guide relationships to be sure the different guides get the correct scores.
    for key, guide in enumerate(guides):
        for t in results:
            if t[0] == key:
                _, guide.repair_profile, guide.repair_stats = t

    return guides