from Bio.Seq import Seq
from Bio.SeqUtils import GC

import config
from constants import codes, STEM_LEN


def perm_PAM(PAM):
    PAM = PAM.upper()
    new_comb = [""]  # in case no PAM
    if len(PAM) == 1:
        new_comb = codes[PAM]

    for i in range(len(PAM) - 1):
        if i == 0:
            comb = codes[PAM[0]]
            new_comb = []
        else:
            comb = new_comb
            new_comb = []

        for j in codes[PAM[i + 1]]:
            for c in comb:
                new_comb.append(c + j)

    return new_comb


def compare_PAM(basePAM, baseDNA):
    if basePAM == "N":
        return True

    if basePAM == baseDNA:
        return True

    if basePAM == "W" and (baseDNA == "A" or baseDNA == "T"):
        return True

    if basePAM == "S" and (baseDNA == "C" or baseDNA == "G"):
        return True

    if basePAM == "M" and (baseDNA == "A" or baseDNA == "C"):
        return True

    if basePAM == "K" and (baseDNA == "G" or baseDNA == "T"):
        return True

    if basePAM == "R" and (baseDNA == "A" or baseDNA == "G"):
        return True

    if basePAM == "Y" and (baseDNA == "C" or baseDNA == "T"):
        return True

    if basePAM == "B" and baseDNA != "A":
        return True

    if basePAM == "D" and baseDNA != "C":
        return True

    if basePAM == "H" and baseDNA != "G":
        return True

    if basePAM == "V" and baseDNA != "T":
        return True

    return False


# Used in Cas9 and selfComp
def gc_content(seq):
    gc = 0
    for i in seq:
        if i == 'G' or i == 'g' or i == 'C' or i == 'c':
            gc += 1
    return float(gc) / float(len(seq))


def self_comp(fwd, backbone):
    rvs = str(fwd.reverse_complement())
    fwd = str(fwd)
    L = len(fwd) - STEM_LEN - 1
    folding = 0
    for i in range(0, len(fwd) - STEM_LEN):
        if gc_content(fwd[i:i + STEM_LEN]) >= 0.5:
            if fwd[i:i + STEM_LEN] in rvs[0:(L - i)] or any(
                    [fwd[i:i + STEM_LEN] in item for item in backbone]):
                folding += 1

    return folding


def eval_CPF1_sequence(name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, PAM,
                       filterGCmin, filterGCmax, filterSelfCompMax, replace5prime=None, backbone=None):
    """ Evaluates an k-mer as a potential Cpf1 target site """

    gLen = guideSize-len(PAM)
    revCompPAM = str(Seq(PAM).reverse_complement())
    dna = Seq(dna)

    if replace5prime:
        fwd = dna[len(PAM):-len(replace5prime)] + replace5prime  # Replace the 2 first bases with e.g. "GG"
    else:
        fwd = dna[len(PAM):]  # Do not include PAM motif in folding calculations

    add = True
    for pos in range(len(PAM)):
        if compare_PAM(PAM[pos], dna[pos]):
            continue
        else:
            add = False
            break

    if add and (filterGCmin != 0 or filterGCmax != 100):
        gc = GC(dna[len(PAM):])
        if gc < filterGCmin or gc > filterGCmax:
            add = False

    if add and filterSelfCompMax != -1:
        if replace5prime:
            fwd = replace5prime + dna[len(PAM):-len(replace5prime)]
        else:
            fwd = dna[len(PAM):]
        folding = self_comp(fwd, backbone)
        if folding > filterSelfCompMax:
            add = False

    if add:
        if config.isoforms:
            pam_comb = perm_PAM(PAM)
            for p in pam_comb:
                fastaFile.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                    name, num, num + guideSize, downstream5prim, downstream3prim,
                    dna, p, p + dna[len(PAM):]))
        else:
            dna = dna.reverse_complement()
            pam_comb = perm_PAM(revCompPAM)
            for p in pam_comb:
                fastaFile.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                                name, num, num+guideSize, downstream5prim, downstream3prim,
                                dna, p, dna[:gLen] + p))
        return True

    add = True and not config.isoforms

    for pos in range(len(PAM)):
        if compare_PAM(revCompPAM[pos], dna[gLen + pos]):
            continue
        else:
            add = False
            break

    if add and (filterGCmin != 0 or filterGCmax != 100):
        gc = GC(dna.reverse_complement()[len(PAM):])
        if gc < filterGCmin or gc > filterGCmax:
            add = False

    if add and filterSelfCompMax != -1:
        if replace5prime:
            fwd = replace5prime + dna.reverse_complement()[len(PAM):-len(replace5prime)]
        else:
            fwd = dna.reverse_complement()[len(PAM):]
        folding = self_comp(fwd, backbone)
        if folding > filterSelfCompMax:
            add = False

    if add:
        pam_comb = perm_PAM(revCompPAM)
        for p in pam_comb:
            # on the reverse strand seq of 5' downstream becomes 3' downstream and vice versa
            fastaFile.write('>%s_%d-%d:%s:%s:-:%s:%s\n%s\n' % (
                            name, num, num+guideSize,
                            Seq(downstream3prim).reverse_complement(),
                            Seq(downstream5prim).reverse_complement(),
                            dna, p, dna[:gLen] + p))
        return True

    return False


def eval_CRISPR_sequence(name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, allowed, PAM,
                         filterGCmin, filterGCmax, filterSelfCompMax, replace5prime=None, backbone=None):
    """ Evaluates an k-mer as a potential CRISPR target site """

    gLen = guideSize-len(PAM)
    revCompPAM = str(Seq(PAM).reverse_complement())
    dna = Seq(dna)

    if str(dna[0:2]) in allowed:
        add = True
        for pos in range(len(PAM)):
            if compare_PAM(PAM[pos], dna[gLen + pos]):
                continue
            else:
                add = False
                break

        if add and (filterGCmin != 0 or filterGCmax != 100):
            # TODO EVERYWHERE GC content does not assumes 5' replacement
            gc = GC(dna[0:(None if PAM == "" else -len(PAM))])
            if gc < filterGCmin or gc > filterGCmax:
                add = False

        if add and filterSelfCompMax != -1:
            if replace5prime:
                fwd = replace5prime + dna[len(replace5prime):(None if PAM == "" else -len(PAM))]
            else:
                fwd = dna[0:(None if PAM == "" else -len(PAM))]
            folding = self_comp(fwd, backbone)
            if folding > filterSelfCompMax:
                add = False

        # in order to control the number of mismatches to search in the last 8 or 3 bps,
        # need to reverse complement so the seed region can be at the start
        # rather than end of the sequence
        # not in isoforms case as we don't search reverse complement
        if add:
            if config.isoforms:
                pam_comb = perm_PAM(PAM)
                for p in pam_comb:
                    fastaFile.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                        name, num, num + guideSize, downstream5prim, downstream3prim,
                        dna, p, dna[:gLen] + p))
                return True
            else:
                # all combinations of possible PAMs
                dna = dna.reverse_complement()
                pam_comb = perm_PAM(revCompPAM)
                for p in pam_comb:
                    fastaFile.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                                    name, num, num+guideSize, downstream5prim, downstream3prim,
                                    dna, p, p + dna[len(revCompPAM):]))
                return True

    if str(dna[-2:].reverse_complement()) in allowed and not config.isoforms:
        add = True

        for pos in range(len(PAM)):
            if compare_PAM(revCompPAM[pos], dna[pos]):
                continue
            else:
                add = False
                break

        if add and (filterGCmin != 0 or filterGCmax != 100):
            gc = GC(dna[len(PAM):])
            if gc < filterGCmin or gc > filterGCmax:
                add = False

        if add and filterSelfCompMax != -1:
            if replace5prime:
                fwd = replace5prime + dna.reverse_complement()[len(PAM):-len(replace5prime)]
            else:
                fwd = dna.reverse_complement()[len(PAM):]
            folding = self_comp(fwd, backbone)
            if folding > filterSelfCompMax:
                add = False

        if add:
            pam_comb = perm_PAM(revCompPAM)
            for p in pam_comb:
                # on the reverse strand seq of 5' downstream becomes 3' downstream and vice versa
                fastaFile.write('>%s_%d-%d:%s:%s:-:%s:%s\n%s\n' % (
                                name, num, num+guideSize,
                                Seq(downstream3prim).reverse_complement(),
                                Seq(downstream5prim).reverse_complement(),
                                dna, p, p + dna[len(revCompPAM):]))
            return True

    return False


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


__all__ = ["eval_TALENS_sequence", "eval_CRISPR_sequence", "eval_CPF1_sequence", "gc_content", "compare_PAM",
           "perm_PAM"]
