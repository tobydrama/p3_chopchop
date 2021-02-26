
#####################
##
## TALEN SPECIFIC FUNCTIONS
##
from Vars import TALEN_OFF_TARGET_MIN, TALEN_OFF_TARGET_MAX, PRIMER_OFF_TARGET_MIN
from classes import Nickase
from operator import attrgetter
from classes.PAIR import Pair


def pairTalens(taleList, fastaSeq, guideSize, taleMinDistance, taleMaxDistance,
               enzymeCo, maxOffTargets, g_RVD, minResSiteLen):
    pairs = []

    for i in range(len(taleList)-1):
        tale1 = taleList[i]

        # FIX: Only start looking for pair when > 30 - 36 spacer+length of i-TALE (modified for 17-mers and 18-mers)
        for j in range(i+1, len(taleList)-1):
            tale2 = taleList[j]

            # This will finish the search for more pairs if we are out of range
            if tale1.start + taleMaxDistance < tale2.start:
                break

            elif tale1.start + taleMinDistance < tale2.start and tale1.guideSeq[0] == "T" and \
                            tale2.guideSeq[guideSize-1] == "A":

                # EDV: Are all these find calls faster than a regular expression?
                pos = tale1.name.find('_')
                exon1 = tale1.name[:pos]
                exonSeq = fastaSeq[exon1]

                # Make sure the two TALENs are on the same "slice", only a problem for overlapping padding regions
                pos2 = tale2.name.find('_')
                exon2 = tale2.name[:pos2]
                if exon1 != exon2:
                    continue

                # The coordinates of the tale within the exon e.g. 128-143
                tale1coords = tale1.name[pos+1:]

                # Just the second coordinate, corresponding to the end of the first tale e.g. 143
                tale1End = int(tale1coords[tale1coords.find('-')+1:])

                # The coordinates of the tale within the exon e.g. 160-175
                tale2coords = tale2.name[tale2.name.find('_')+1:]

                # Just the first coordinate, corresponding to the beginning of the second tale e.g. 160
                tale2Start = int(tale2coords[:tale2coords.find('-')])

                # sequence of spacer between end of tale1 and beginning of tale2
                spacerSeq = exonSeq[tale1End:tale2Start]

                spacerSize = len(spacerSeq)

                # if spacerSize < 3:
                #      sys.stderr.write("(%s)  (%s)\n" % (tale1.name, tale2.name))
                #      sys.stderr.write("(%s)  (%s)\n" % (e1, e2))
                #      sys.stderr.write("%s-%s\n" % (tale1End, tale2Start))
                #      sys.stderr.write("%s\t%s\t%s\n" % (tale1.guideSeq, spacerSeq, tale2.guideSeq))
                #      sys.exit()

                # Calculates off-target pairs for tale1 and tale2 (see below)
                offTargetPairs = has_Off_targets(tale1, tale2, TALEN_OFF_TARGET_MIN, TALEN_OFF_TARGET_MAX)

                # Makes tale1 and tale2 into a Pair object, and adds to list of Pair objects
                pairs.append(Pair(tale1, tale2, spacerSeq, spacerSize, offTargetPairs, enzymeCo, maxOffTargets, g_RVD,
                                  minResSiteLen))

    return pairs


def pairCas9(taleList, fastaSeq, guideSize, taleMinDistance, taleMaxDistance, enzymeCo, maxOffTargets, minResSiteLen,
             offtargetMaxDist):
    pairs = []

    for i in range(len(taleList)-1):
        tale1 = taleList[i]

        # FIX: Only start looking for pair when > 30 - 36 spacer+length of i-TALE (modified for 17-mers and 18-mers)
        for j in range(i+1, len(taleList)-1):
            tale2 = taleList[j]

            if tale1.start + taleMaxDistance < tale2.start:
                continue

            elif tale1.start + taleMinDistance < tale2.start and tale1.strand != tale2.strand:

                # EDV: Are all these find calls faster than a regular expression?
                pos = tale1.name.rfind('_')
                exon1 = tale1.name[:pos]
                exonSeq = fastaSeq[exon1]

                # Make sure the two TALENs are on the same "slice", only a problem for overlapping padding regions
                pos2 = tale2.name.rfind('_')
                exon2 = tale2.name[:pos2]
                if exon1 != exon2:
                    continue

                # The coordinates of the tale within the exon e.g. 128-143
                tale1coords = tale1.name[pos+1:]

                # Just the second coordinate, corresponding to the end of the first tale e.g. 143
                tale1End = int(tale1coords[tale1coords.rfind('-')+1:])

                # The coordinates of the tale within the exon e.g. 160-175
                tale2coords = tale2.name[tale2.name.rfind('_')+1:]

                # Just the first coordinate, corresponding to the beginning of the second tale e.g. 160
                tale2Start = int(tale2coords[:tale2coords.rfind('-')])

                # sequence of spacer between end of tale1 and beginning of tale2
                spacerSeq = exonSeq[tale1End:tale2Start]

                spacerSize = len(spacerSeq)

                # Calculates off-target pairs for tale1 and tale2 (see below)
                offTargetPairs = has_Off_targets(tale1, tale2, taleMinDistance-guideSize, offtargetMaxDist)

                # Makes tale1 and tale2 into a Pair object, and adds to list of Pair objects
                pairs.append(Nickase(tale1, tale2, spacerSeq, spacerSize, offTargetPairs, enzymeCo, maxOffTargets,
                                     minResSiteLen))

    return pairs


def has_Off_targets(tale1, tale2, offTargetMin, offTargetMax):
    """ Returns the number of off-targets for a pair of TALENs (10-24bp apart) """

    offTargetPairs = []

    # Calls sort function to sort off-targets by chromosome and chromosome position.
    # Bowtie ranks them according to quality of hit
    tale1.sort_offTargets()
    tale2.sort_offTargets()

    ### FIX: Eivind to write this code properly. Include a way to step backwards, so as not to miss any hits.
    # Need to make a queue..?
    for i in range(len(tale1.offTargets)):
        hit1 = tale1.offTargets[i]

        for j in range(len(tale2.offTargets)):
            hit2 = tale2.offTargets[j]

            # Determines whether 2 tales are on the same chromosome and 10-24 bp apart.
            if hit2.chrom == hit1.chrom and offTargetMin <= abs(hit2.start-hit1.start) <= offTargetMax:
                offTargetPairs.append([hit1, hit2])

    return offTargetPairs


def clusterPairs(pairs):
    """ Clusters paired sequences according to overlap, so user knows which TALE pairs are redundant """

    # Sets the starting pair of TALEs to be compared to
    first = pairs[0]
    cluster = 1
    first.cluster = cluster
    inCluster = 0

    # Compares each TALE pair to previous pair in list to see whether redundant. Assigns cluster number accordingly
    for i in range(1,len(pairs)):
        cur = pairs[i]
        prev = pairs[i-1]

        # Specifically, compares location of spacer (by comparing location of tales) to see whether there is overlap,
        # and therefore TALE pairs are redundant
        if ((cur.spacerStart <= prev.spacerEnd) and (cur.spacerEnd >= prev.spacerStart) and
                    inCluster < PRIMER_OFF_TARGET_MIN):

            cur.cluster = cluster
            inCluster += 1
        else:
            # If not redundant, increase cluster number
            cluster += 1
            cur.cluster = cluster
            inCluster = 0

    return (cluster, pairs)


def eval_TALENS_sequence(name, targetSize, dna, num, fastaFile, downstream5prim, downstream3prim):
    """ Evaluates an N-mer as a potential TALENs target site """
    del downstream5prim, downstream3prim
    found = False
    if dna[0] == "T":
        # dna = Seq(dna).reverse_complement()
        fastaFile.write('>%s_%d-%d\n%s\n' % (name, num, num+targetSize, dna))
        found = True
    elif dna[-1] == "A":
        fastaFile.write('>%s_%d-%d\n%s\n' % (name, num, num+targetSize, dna))
        found = True

    return found


def sort_TALEN_pairs(pairs):
    """ Sort pairs according to score and cluster """

    return sorted(pairs, key=attrgetter('score', 'cluster'))
