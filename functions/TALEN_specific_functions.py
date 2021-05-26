from classes.Nickase import Nickase
from classes.PAIR import Pair
from constants import TALEN_OFF_TARGET_MIN, TALEN_OFF_TARGET_MAX, PRIMER_OFF_TARGET_MIN
from functions.make_primers import has_off_targets


def pair_talens(tale_list, fasta_seq, guide_size, tale_min_distance, tale_max_distance,
                enzyme_co, max_off_targets, g_rvd, min_res_site_len):
    pairs = []

    for i in range(len(tale_list) - 1):
        tale1 = tale_list[i]

        # FIX: Only start looking for pair when > 30 - 36 spacer+length of i-TALE (modified for 17-mers and 18-mers)
        for j in range(i + 1, len(tale_list) - 1):
            tale2 = tale_list[j]

            # This will finish the search for more pairs if we are out of range
            if tale1.start + tale_max_distance < tale2.start:
                break

            elif tale1.start + tale_min_distance < tale2.start and tale1.guide_seq[0] == "T" and \
                    tale2.guide_seq[guide_size - 1] == "A":

                # EDV: Are all these find calls faster than a regular expression?
                pos = tale1.name.find('_')
                exon1 = tale1.name[:pos]
                exon_seq = fasta_seq[exon1]

                # Make sure the two TALENs are on the same "slice", only a problem for overlapping padding regions
                pos2 = tale2.name.find('_')
                exon2 = tale2.name[:pos2]
                if exon1 != exon2:
                    continue

                # The coordinates of the tale within the exon e.g. 128-143
                tale1_coords = tale1.name[pos + 1:]

                # Just the second coordinate, corresponding to the end of the first tale e.g. 143
                tale1_end = int(tale1_coords[tale1_coords.find('-') + 1:])

                # The coordinates of the tale within the exon e.g. 160-175
                tale2_coords = tale2.name[tale2.name.find('_') + 1:]

                # Just the first coordinate, corresponding to the beginning of the second tale e.g. 160
                tale2_start = int(tale2_coords[:tale2_coords.find('-')])

                # sequence of spacer between end of tale1 and beginning of tale2
                spacer_seq = exon_seq[tale1_end:tale2_start]

                spacer_size = len(spacer_seq)

                # if spacer_size < 3:
                #      sys.stderr.write("(%s)  (%s)\n" % (tale1.name, tale2.name))
                #      sys.stderr.write("(%s)  (%s)\n" % (e1, e2))
                #      sys.stderr.write("%s-%s\n" % (tale1_end, tale2_start))
                #      sys.stderr.write("%s\t%s\t%s\n" % (tale1.guideSeq, spacer_seq, tale2.guideSeq))
                #      sys.exit()

                # Calculates off-target pairs for tale1 and tale2 (see below)
                off_target_pairs = has_off_targets(tale1, tale2, TALEN_OFF_TARGET_MIN, TALEN_OFF_TARGET_MAX)

                # Makes tale1 and tale2 into a Pair object, and adds to list of Pair objects
                pairs.append(Pair(tale1, tale2, spacer_seq, spacer_size, off_target_pairs, enzyme_co, max_off_targets,
                                  g_rvd, min_res_site_len))

    return pairs


def pair_cas9(tale_list, fasta_seq, guide_size, tale_min_distance, tale_max_distance, enzyme_co, max_off_targets,
              min_res_site_len, offtarget_max_dist):
    pairs = []

    for i in range(len(tale_list) - 1):
        tale1 = tale_list[i]

        # FIX: Only start looking for pair when > 30 - 36 spacer+length of i-TALE (modified for 17-mers and 18-mers)
        for j in range(i + 1, len(tale_list) - 1):
            tale2 = tale_list[j]

            if tale1.start + tale_max_distance < tale2.start:
                continue

            elif tale1.start + tale_min_distance < tale2.start and tale1.strand != tale2.strand:

                # EDV: Are all these find calls faster than a regular expression?
                pos = tale1.name.rfind('_')
                exon1 = tale1.name[:pos]
                exon_seq = fasta_seq[exon1]

                # Make sure the two TALENs are on the same "slice", only a problem for overlapping padding regions
                pos2 = tale2.name.rfind('_')
                exon2 = tale2.name[:pos2]
                if exon1 != exon2:
                    continue

                # The coordinates of the tale within the exon e.g. 128-143
                tale1_coords = tale1.name[pos + 1:]

                # Just the second coordinate, corresponding to the end of the first tale e.g. 143
                tale1_end = int(tale1_coords[tale1_coords.rfind('-') + 1:])

                # The coordinates of the tale within the exon e.g. 160-175
                tale2_coords = tale2.name[tale2.name.rfind('_') + 1:]

                # Just the first coordinate, corresponding to the beginning of the second tale e.g. 160
                tale2_start = int(tale2_coords[:tale2_coords.rfind('-')])

                # sequence of spacer between end of tale1 and beginning of tale2
                spacer_seq = exon_seq[tale1_end:tale2_start]

                spacer_size = len(spacer_seq)

                # Calculates off-target pairs for tale1 and tale2 (see below)
                off_target_pairs = has_off_targets(tale1, tale2, tale_min_distance - guide_size, offtarget_max_dist)

                # Makes tale1 and tale2 into a Pair object, and adds to list of Pair objects
                pairs.append(
                    Nickase(tale1, tale2, spacer_seq, spacer_size, off_target_pairs, enzyme_co, max_off_targets,
                            min_res_site_len))

    return pairs


def cluster_pairs(pairs):
    """ Clusters paired sequences according to overlap, so user knows which TALE pairs are redundant """

    # Sets the starting pair of TALEs to be compared to
    first = pairs[0]
    cluster = 1
    first.cluster = cluster
    in_cluster = 0

    # Compares each TALE pair to previous pair in list to see whether redundant. Assigns cluster number accordingly
    for i in range(1, len(pairs)):
        cur = pairs[i]
        prev = pairs[i - 1]

        # Specifically, compares location of spacer (by comparing location of tales) to see whether there is overlap,
        # and therefore TALE pairs are redundant
        if ((cur.spacer_start <= prev.spacer_end) and (cur.spacer_end >= prev.spacer_start) and
                in_cluster < PRIMER_OFF_TARGET_MIN):

            cur.cluster = cluster
            in_cluster += 1
        else:
            # If not redundant, increase cluster number
            cluster += 1
            cur.cluster = cluster
            in_cluster = 0

    return cluster, pairs
