import codecs
import pickle
import subprocess

import config
from classes.CPF1 import Cpf1


def convert_CPF1_to_tuple(key: int, guide: Cpf1) -> (int, str, str, str, float, float):
    return (key,
            guide.downstream_5_prim,
            guide.downstream_3_prim,
            guide.stranded_guide_seq,
            guide.score,
            guide.coefficients_score)


def run_kim_2018(guides: [Cpf1]) -> [Cpf1]:
    """
    Runs chopchop_KIM_2018 docker image using the supplied guides & scoring method.

    :param scoring_method: The scoring method to use. Accepted values are "KIM_2018" & "ALL".
    :param guides: A list of Cpf1 objects to score.
    :return: Returns a list of Cpf1 scored objects.
    """
    keyed_tuples = []
    for key, guide in enumerate(guides):
        keyed_tuples.append(convert_CPF1_to_tuple(key, guide))

    encoded = codecs.encode(pickle.dumps(keyed_tuples, protocol=2), 'base64').decode()

    command = ['docker', 'run', '-i', 'chopchop_kim_2018', '-c', str(config.score('COEFFICIENTS'))]

    kim_2018 = subprocess.run(command, capture_output=True, text=True, input=encoded)

    # encoding='latin1' is for backwards compatibility.
    results = pickle.loads(codecs.decode(kim_2018.stdout.encode(), 'base64'), encoding='latin1')
    # TODO currently we loop through key, guide pairs to set results, as we did for the keyed_tuples. We might want to
    #  save keyed_tuple & guide relationships to be sure the different guides get the correct scores.
    for key, guide in enumerate(guides):
        for t in results:
            if t[0] == key:
                _, guide.score, guide.coefficients_score = t

    return guides
