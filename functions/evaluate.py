from Bio.Seq import Seq
from Bio.SeqUtils import GC

import config
from constants import codes, STEM_LEN


def perm_pam(pam):
    pam = pam.upper()
    new_comb = [""]  # in case no PAM
    if len(pam) == 1:
        new_comb = codes[pam]

    for i in range(len(pam) - 1):
        if i == 0:
            comb = codes[pam[0]]
            new_comb = []
        else:
            comb = new_comb
            new_comb = []

        for j in codes[pam[i + 1]]:
            for c in comb:
                new_comb.append(c + j)

    return new_comb


def compare_pam(base_pam, base_dna):
    if base_pam == "N":
        return True

    if base_pam == base_dna:
        return True

    if base_pam == "W" and (base_dna == "A" or base_dna == "T"):
        return True

    if base_pam == "S" and (base_dna == "C" or base_dna == "G"):
        return True

    if base_pam == "M" and (base_dna == "A" or base_dna == "C"):
        return True

    if base_pam == "K" and (base_dna == "G" or base_dna == "T"):
        return True

    if base_pam == "R" and (base_dna == "A" or base_dna == "G"):
        return True

    if base_pam == "Y" and (base_dna == "C" or base_dna == "T"):
        return True

    if base_pam == "B" and base_dna != "A":
        return True

    if base_pam == "D" and base_dna != "C":
        return True

    if base_pam == "H" and base_dna != "G":
        return True

    if base_pam == "V" and base_dna != "T":
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
    l = len(fwd) - STEM_LEN - 1
    folding = 0
    for i in range(0, len(fwd) - STEM_LEN):
        if gc_content(fwd[i:i + STEM_LEN]) >= 0.5:
            if fwd[i:i + STEM_LEN] in rvs[0:(l - i)] or any(
                    [fwd[i:i + STEM_LEN] in item for item in backbone]):
                folding += 1

    return folding


def eval_cpf1_sequence(name, guide_size, dna, num, fasta_file, downstream_5prim, downstream_3prim, pam,
                       filter_gc_min, filter_gc_max, filter_self_comp_max, replace_5prime=None, backbone=None):
    """ Evaluates an k-mer as a potential Cpf1 target site """

    g_len = guide_size - len(pam)
    rev_comp_pam = str(Seq(pam).reverse_complement())
    dna = Seq(dna)

    if replace_5prime:
        fwd = dna[len(pam):-len(replace_5prime)] + replace_5prime  # Replace the 2 first bases with e.g. "GG"
    else:
        fwd = dna[len(pam):]  # Do not include PAM motif in folding calculations

    add = True
    for pos in range(len(pam)):
        if compare_pam(pam[pos], dna[pos]):
            continue
        else:
            add = False
            break

    if add and (filter_gc_min != 0 or filter_gc_max != 100):
        gc = GC(dna[len(pam):])
        if gc < filter_gc_min or gc > filter_gc_max:
            add = False

    if add and filter_self_comp_max != -1:
        if replace_5prime:
            fwd = replace_5prime + dna[len(pam):-len(replace_5prime)]
        else:
            fwd = dna[len(pam):]
        folding = self_comp(fwd, backbone)
        if folding > filter_self_comp_max:
            add = False

    if add:
        if config.isoforms:
            pam_comb = perm_pam(pam)
            for p in pam_comb:
                fasta_file.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                    name, num, num + guide_size, downstream_5prim, downstream_3prim,
                    dna, p, p + dna[len(pam):]))
        else:
            dna = dna.reverse_complement()
            pam_comb = perm_pam(rev_comp_pam)
            for p in pam_comb:
                fasta_file.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                    name, num, num + guide_size, downstream_5prim, downstream_3prim,
                    dna, p, dna[:g_len] + p))
        return True

    add = True and not config.isoforms

    for pos in range(len(pam)):
        if compare_pam(rev_comp_pam[pos], dna[g_len + pos]):
            continue
        else:
            add = False
            break

    if add and (filter_gc_min != 0 or filter_gc_max != 100):
        gc = GC(dna.reverse_complement()[len(pam):])
        if gc < filter_gc_min or gc > filter_gc_max:
            add = False

    if add and filter_self_comp_max != -1:
        if replace_5prime:
            fwd = replace_5prime + dna.reverse_complement()[len(pam):-len(replace_5prime)]
        else:
            fwd = dna.reverse_complement()[len(pam):]
        folding = self_comp(fwd, backbone)
        if folding > filter_self_comp_max:
            add = False

    if add:
        pam_comb = perm_pam(rev_comp_pam)
        for p in pam_comb:
            # on the reverse strand seq of 5' downstream becomes 3' downstream and vice versa
            fasta_file.write('>%s_%d-%d:%s:%s:-:%s:%s\n%s\n' % (
                name, num, num + guide_size,
                Seq(downstream_3prim).reverse_complement(),
                Seq(downstream_5prim).reverse_complement(),
                dna, p, dna[:g_len] + p))
        return True

    return False


def eval_crispr_sequence(name, guide_size, dna, num, fasta_file, downstream_5prim, downstream_3prim, allowed, pam,
                         filter_gc_min, filter_gc_max, filter_self_comp_max, replace_5prime=None, backbone=None):
    """ Evaluates an k-mer as a potential CRISPR target site """

    g_len = guide_size - len(pam)
    rev_comp_pam = str(Seq(pam).reverse_complement())
    dna = Seq(dna)

    if str(dna[0:2]) in allowed:
        add = True
        for pos in range(len(pam)):
            if compare_pam(pam[pos], dna[g_len + pos]):
                continue
            else:
                add = False
                break

        if add and (filter_gc_min != 0 or filter_gc_max != 100):
            # TODO EVERYWHERE GC content does not assumes 5' replacement
            gc = GC(dna[0:(None if pam == "" else -len(pam))])
            if gc < filter_gc_min or gc > filter_gc_max:
                add = False

        if add and filter_self_comp_max != -1:
            if replace_5prime:
                fwd = replace_5prime + dna[len(replace_5prime):(None if pam == "" else -len(pam))]
            else:
                fwd = dna[0:(None if pam == "" else -len(pam))]
            folding = self_comp(fwd, backbone)
            if folding > filter_self_comp_max:
                add = False

        # in order to control the number of mismatches to search in the last 8 or 3 bps,
        # need to reverse complement so the seed region can be at the start
        # rather than end of the sequence
        # not in isoforms case as we don't search reverse complement
        if add:
            if config.isoforms:
                pam_comb = perm_pam(pam)
                for p in pam_comb:
                    fasta_file.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                        name, num, num + guide_size, downstream_5prim, downstream_3prim,
                        dna, p, dna[:g_len] + p))
                return True
            else:
                # all combinations of possible PAMs
                dna = dna.reverse_complement()
                pam_comb = perm_pam(rev_comp_pam)
                for p in pam_comb:
                    fasta_file.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                        name, num, num + guide_size, downstream_5prim, downstream_3prim,
                        dna, p, p + dna[len(rev_comp_pam):]))
                return True

    if str(dna[-2:].reverse_complement()) in allowed and not config.isoforms:
        add = True

        for pos in range(len(pam)):
            if compare_pam(rev_comp_pam[pos], dna[pos]):
                continue
            else:
                add = False
                break

        if add and (filter_gc_min != 0 or filter_gc_max != 100):
            gc = GC(dna[len(pam):])
            if gc < filter_gc_min or gc > filter_gc_max:
                add = False

        if add and filter_self_comp_max != -1:
            if replace_5prime:
                fwd = replace_5prime + dna.reverse_complement()[len(pam):-len(replace_5prime)]
            else:
                fwd = dna.reverse_complement()[len(pam):]
            folding = self_comp(fwd, backbone)
            if folding > filter_self_comp_max:
                add = False

        if add:
            pam_comb = perm_pam(rev_comp_pam)
            for p in pam_comb:
                # on the reverse strand seq of 5' downstream becomes 3' downstream and vice versa
                fasta_file.write('>%s_%d-%d:%s:%s:-:%s:%s\n%s\n' % (
                    name, num, num + guide_size,
                    Seq(downstream_3prim).reverse_complement(),
                    Seq(downstream_5prim).reverse_complement(),
                    dna, p, p + dna[len(rev_comp_pam):]))
            return True

    return False


def eval_talens_sequence(name, target_size, dna, num, fasta_file, downstream_5_prim, downstream_3_prim):
    """ Evaluates an N-mer as a potential TALENs target site """
    if dna[0] == "T" or dna[-1] == "A":
        fasta_file.write('>%s_%d-%d\n%s\n' % (name, num, num + target_size, dna))
        return True

    return False


__all__ = ["eval_talens_sequence", "eval_crispr_sequence", "eval_cpf1_sequence", "gc_content", "compare_pam",
           "perm_pam"]
