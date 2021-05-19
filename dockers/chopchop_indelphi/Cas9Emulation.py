class Cas9Emulation(object):
    """
    `Cas9` emulation.

    This class emulates the `Cas9` class from the main CHOPCHOP script, guides are transmitted through the argument
    parser using the `Cas9EmulationAction` class.
    """

    def __init__(self, key, downstream_5_prim, downstream_3_prim, stranded_guide_seq, pam, rep_profile, rep_stats):
        self.CoefficientsScore = {}

        # Unique identifier for each Cas9 object
        self.key = key

        # Cas9 fields used by repair predictions
        self.downstream5prim = downstream_5_prim.encode('ascii', 'ignore')
        self.downstream3prim = downstream_3_prim.encode('ascii', 'ignore')
        self.strandedGuideSeq = stranded_guide_seq.encode('ascii', 'ignore')
        self.PAM = pam.encode('ascii', 'ignore')
        self.repProfile = rep_profile
        self.repStats = rep_stats

