from classes import Guide


class Cpf1(Guide):
    def __init__(self, *args, **kwargs):
        super(Cpf1, self).__init__(*args, **kwargs)
        self.CoefficientsScore = 0  # KIM_2018

    def __str__(self):
        self.sort_offTargets()
        return "%s\t%s:%s\t%s\t%.0f\t%s\t%.0f\t%s\t%s\t%s\t%s" % (self.strandedGuideSeq, self.chrom, self.start,
                                                                  self.strand, self.GCcontent, self.folding,
                                                                  self.CoefficientsScore,
                                                                  self.offTargetsMM[0], self.offTargetsMM[1],
                                                                  self.offTargetsMM[2],
                                                                  ">=" + str(self.offTargetsMM[3]) if self.isKmaxed else
                                                                  self.offTargetsMM[3])
