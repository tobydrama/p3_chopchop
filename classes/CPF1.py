from classes.Guide import Guide


class Cpf1(Guide):
    def __init__(self, *args, **kwargs):
        super(Cpf1, self).__init__(*args, **kwargs)
        self.coefficients_score = 0  # KIM_2018

    def __str__(self):
        self.sort_off_targets()
        return "%s\t%s:%s\t%s\t%.0f\t%s\t%.0f\t%s\t%s\t%s\t%s" % (
            self.stranded_guide_seq, self.chrom, self.start, self.strand, self.gc_content, self.folding,
            self.coefficients_score, self.off_targets_mm[0], self.off_targets_mm[1], self.off_targets_mm[2],
            ">=" + str(self.off_targets_mm[3]) if self.is_kmaxed else self.off_targets_mm[3])
