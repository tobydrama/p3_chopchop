class Cas9Emulation(object):
    """
    `Cas9` emulation.

    This class emulates the `Cas9` class from the main CHOPCHOP script, guides are transmitted through the argument
    parser using the `Cas9EmulationAction` class.
    """

    def __init__(self, key, downstream_5_prim, downstream_3_prim, stranded_guide_seq, pam, score, coefficient_score):
        self.CoefficientsScore = {}

        # Unique identifier for each Cas9 object
        self.key = key
        # Cas9 fields used by the Doench 2016 algorithm
        self.downstream5prim = downstream_5_prim
        self.downstream3prim = downstream_3_prim
        self.strandedGuideSeq = stranded_guide_seq
        self.PAM = pam
        self.score = score
        self.CoefficientsScore["DOENCH_2016"] = coefficient_score

