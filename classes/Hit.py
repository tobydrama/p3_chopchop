class Hit:
    """Creates class for each hit from bowtie."""

    def __init__(self, line):
        self.flagSum = int(line[1])
        self.chrom = line[2]
        self.start = int(line[3])
        self.matchSeq = line[9]
        self.mismatch = line[-1]
        self.mismatchPos = line[-2]
        self.opts = line[11:(len(line))]
        self.mismatchCorrected = False

    def calc_mismatchPos (self):
        """ Updates the sequence parsed from the SAM output to include the mismatches """

        lastDigit = len(self.mismatchPos)-1
        guideSize = len(self.matchSeq)
        guideCurr = ""

        ## MD:Z:GUIDESIZE means that there are no mismatches
        if not(self.mismatchPos =="MD:Z:%s" % guideSize):
            guideIndex = 0
            currTotal = 0

            for c in range(5, lastDigit+1):

                # If the character is a digit, check if the next character is a digit (>9) and add number to total
                if self.mismatchPos[c].isdigit():

                    if c != lastDigit and self.mismatchPos[c+1].isdigit():
                        currTotal += (int(self.mismatchPos[c])*10)
                    else:
                        currTotal += int(self.mismatchPos[c])
                        guideCurr += self.matchSeq[guideIndex:currTotal]
                        guideIndex = currTotal

                # if character is a letter, add one to total
                else:
                    guideCurr += self.mismatchPos[c].lower()
                    currTotal += 1
                    guideIndex += 1


            self.matchSeq = guideCurr


    # Specifying how to print items in list of off-targets
    def __str__(self):
        if not self.mismatchCorrected:
            self.calc_mismatchPos()
            self.mismatchCorrected = True

        return "%s:%s\t%s" % (self.chrom, self.start, self.matchSeq)

    def asOffTargetString(self, label, maxOffTargets):
        if self.mismatch == "XM:i:%s" % maxOffTargets:
            return "%s,>%s across the genome,0-3,n/a " % (label, maxOffTargets)
        else:
            if not self.mismatchCorrected:
                self.calc_mismatchPos()
                self.mismatchCorrected = True

        return "%s,%s,%s,%s" % (label, self.chrom + ":" + str(self.start), self.mismatch[-1], self.matchSeq)

