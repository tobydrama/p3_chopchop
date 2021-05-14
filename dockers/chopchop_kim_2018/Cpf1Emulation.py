class Cpf1Emulation(object):
    """
    `Cpf1` emulation.

    This class emulates the `Cpf1` class from the main CHOPCHOP script, guides are transmitted through the argument
    parser using the `Cpf1EmulationAction` class.
    """

    def __init__(self, key, downstream_5_prim, downstream_3_prim, stranded_guide_seq, score, coefficient_score):
        self.CoefficientsScore = {}

        # Unique identifier for each Cpf1 object
        self.key = key

        # Cpf1 fields used by the KIM 2018 algorithm
        self.downstream5prim = downstream_5_prim
        self.downstream3prim = downstream_3_prim
        self.strandedGuideSeq = stranded_guide_seq
        self.score = score
        self.CoefficientsScore = coefficient_score

