class Hit:
    """Creates class for each hit from bowtie."""

    def __init__(self, line):
        self.flag_sum = int(line[1])
        self.chrom = line[2]
        self.start = int(line[3])
        self.matchSeq = line[9]
        self.mismatch = line[-1]
        self.mismatch_pos = line[-2]
        self.opts = line[11:(len(line))]
        self.mismatch_corrected = False

    def calc_mismatch_pos(self):
        """ Updates the sequence parsed from the SAM output to include the mismatches """

        last_digit = len(self.mismatch_pos) - 1
        guide_size = len(self.matchSeq)
        guide_curr = ""

        # MD:Z:GUIDESIZE means that there are no mismatches
        if not(self.mismatch_pos == "MD:Z:%s" % guide_size):
            guide_index = 0
            curr_total = 0

            for c in range(5, last_digit+1):

                # If the character is a digit, check if the next character is a digit (>9) and add number to total
                if self.mismatch_pos[c].isdigit():

                    if c != last_digit and self.mismatch_pos[c + 1].isdigit():
                        curr_total += (int(self.mismatch_pos[c]) * 10)
                    else:
                        curr_total += int(self.mismatch_pos[c])
                        guide_curr += self.matchSeq[guide_index:curr_total]
                        guide_index = curr_total

                # if character is a letter, add one to total
                else:
                    guide_curr += self.mismatch_pos[c].lower()
                    curr_total += 1
                    guide_index += 1

            self.matchSeq = guide_curr

    # Specifying how to print items in list of off-targets
    def __str__(self):
        if not self.mismatch_corrected:
            self.calc_mismatch_pos()
            self.mismatch_corrected = True

        return f"{self.chrom}:{self.start}\t{self.matchSeq}"

    def as_off_target_string(self, label, max_off_targets):
        if self.mismatch == "XM:i:%s" % max_off_targets:
            return "%s,>%s across the genome,0-3,n/a " % (label, max_off_targets)
        else:
            if not self.mismatch_corrected:
                self.calc_mismatch_pos()
                self.mismatch_corrected = True

        return f"{label},{self.chrom + ':' + str(self.start)},{self.mismatch[-1]},{self.matchSeq}"
