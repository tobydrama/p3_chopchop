import uuid
from typing import Union

from Bio.Seq import Seq
from Bio.SeqUtils import GC
from typing.io import IO

import config
from classes.CPF1 import Cpf1
import classes.Cas9
from classes.Guide import Guide
from classes.ProgramMode import ProgramMode
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


def compare_PAM(base_PAM, base_DNA):
    if base_PAM == "N":
        return True

    if base_PAM == base_DNA:
        return True

    if base_PAM == "W" and (base_DNA == "A" or base_DNA == "T"):
        return True

    if base_PAM == "S" and (base_DNA == "C" or base_DNA == "G"):
        return True

    if base_PAM == "M" and (base_DNA == "A" or base_DNA == "C"):
        return True

    if base_PAM == "K" and (base_DNA == "G" or base_DNA == "T"):
        return True

    if base_PAM == "R" and (base_DNA == "A" or base_DNA == "G"):
        return True

    if base_PAM == "Y" and (base_DNA == "C" or base_DNA == "T"):
        return True

    if base_PAM == "B" and base_DNA != "A":
        return True

    if base_PAM == "D" and base_DNA != "C":
        return True

    if base_PAM == "H" and base_DNA != "G":
        return True

    if base_PAM == "V" and base_DNA != "T":
        return True

    return []


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


def generate_guide(guide_class, fasta_file: IO, permutation: str, id: str, name: str, start: int, guide_size: int,
                   downstream_5_prime: str, downstream_3_prime: str, strand: str, guide_seq: str,
                   flag_sum: Union[int, None], score_gc: bool, score_self_comp: bool, backbone: str, pam: str,
                   replace_5_prime=None, scoring_method=None, genome=None, gene=None, isoform=None, gene_isoforms=None,
                   is_kmaxed=None, talens: bool = False) -> Guide:
    if talens:
        g_name = f"{name}_{start}-{start + guide_size}"
    else:
        g_name = f"{name}_{start}-{start + guide_size}:{downstream_5_prime}:{downstream_3_prime}:{strand}:{guide_seq}"

    guide = guide_class(
        g_name, flag_sum, guide_size, str(guide_seq), score_gc, score_self_comp, backbone, pam, replace_5_prime,
        scoring_method,
        genome, gene, isoform, gene_isoforms, is_kmaxed, f"{name}_{start}-{start + guide_size}"
    )

    print(f"guide {guide.guideSeq}")

    if talens:
        # fasta_file.write(f">{guide.uuid}\n{str(guide_seq)}\n")
        fasta_file.write(f">{name}_{start}-{start + guide_size}\n{str(guide_seq)}\n")
    else:
        # fasta_file.write(f">{guide.uuid}:{permutation}\n{str(guide_seq)}\n")
        fasta_file.write(f">{name}_{start}-{start + guide_size}:{permutation}\n{str(guide_seq)}\n")

    return guide


def eval_CPF1_sequence(name: str, guide_size: int, dna: str, num: int, fasta_file: IO, downstream_5_prime: str,
                       downstream_3_prime: str, pam: str, filter_gc_min: int, filter_gc_max: int,
                       filter_self_comp_max: int, score_gc: bool, score_self_comp: bool, backbone: str,
                       replace_5_prime: str, scoring_method: ProgramMode, genome: str, gene=None,
                       isoform=False, gene_isoforms=None) -> [Guide]:
    """ Evaluates an k-mer as a potential Cpf1 target site """

    g_len = guide_size - len(pam)
    rev_comp_pam = str(Seq(pam).reverse_complement())
    dna = Seq(dna)

    if replace_5_prime:
        fwd = dna[len(pam):-len(replace_5_prime)] + replace_5_prime  # Replace the 2 first bases with e.g. "GG"
    else:
        fwd = dna[len(pam):]  # Do not include PAM motif in folding calculations

    add = True
    for pos in range(len(pam)):
        if compare_PAM(pam[pos], dna[pos]):
            continue
        else:
            add = False
            break

    if add and (filter_gc_min != 0 or filter_gc_max != 100):
        gc = GC(dna[len(pam):])
        if gc < filter_gc_min or gc > filter_gc_max:
            add = False

    if add and filter_self_comp_max != -1:
        if replace_5_prime:
            fwd = replace_5_prime + dna[len(pam):-len(replace_5_prime)]
        else:
            fwd = dna[len(pam):]
        folding = self_comp(fwd, backbone)
        if folding > filter_self_comp_max:
            add = False

    if add:
        guides = []
        id = uuid.uuid4()
        if config.isoforms:
            pam_comb = perm_PAM(pam)
            for p in pam_comb:
                guides.append(generate_guide(Guide, fasta_file, p, id, name, num, guide_size, downstream_5_prime,
                                             downstream_3_prime, '+', p + dna[len(pam):], None, score_gc,
                                             score_self_comp, backbone, pam, replace_5_prime, scoring_method,
                                             genome, gene, isoform, gene_isoforms, None))
        else:
            dna = dna.reverse_complement()
            pam_comb = perm_PAM(rev_comp_pam)
            for p in pam_comb:
                guides.append(generate_guide(Cpf1, fasta_file, p, id, name, num, guide_size, downstream_5_prime,
                                             downstream_3_prime, '+', dna[:g_len] + p, None, score_gc,
                                             score_self_comp, backbone, pam, replace_5_prime, scoring_method,
                                             genome, gene, isoform, gene_isoforms, None))
        return guides

    add = not config.isoforms
    for pos in range(len(pam)):
        if compare_PAM(rev_comp_pam[pos], dna[g_len + pos]):
            continue
        else:
            add = False
            break

    if add and (filter_gc_min != 0 or filter_gc_max != 100):
        gc = GC(dna.reverse_complement()[len(pam):])
        if gc < filter_gc_min or gc > filter_gc_max:
            add = False

    if add and filter_self_comp_max != -1:
        if replace_5_prime:
            fwd = replace_5_prime + dna.reverse_complement()[len(pam):-len(replace_5_prime)]
        else:
            fwd = dna.reverse_complement()[len(pam):]
        folding = self_comp(fwd, backbone)
        if folding > filter_self_comp_max:
            add = False

    if add:
        guides = []
        id = uuid.uuid4()
        pam_comb = perm_PAM(rev_comp_pam)
        for p in pam_comb:
            # on the reverse strand seq of 5' downstream becomes 3' downstream and vice versa
            guides.append(generate_guide(Cpf1, fasta_file, p, id, name, num, guide_size,
                                         Seq(downstream_3_prime).reverse_complement(),
                                         Seq(downstream_5_prime).reverse_complement(), '-', dna[:g_len] + p, None,
                                         score_gc, score_self_comp, backbone, pam, replace_5_prime, scoring_method,
                                         genome, gene, isoform, gene_isoforms, None))

        return guides

    return []


def eval_CRISPR_sequence(name: str, guide_size: int, dna: str, num: int, fasta_file: IO, downstream_5_prime: str,
                         downstream_3_prime: str, allowed: dict, pam: str, filter_gc_min: int, filter_gc_max: int,
                         filter_self_comp_max: int, score_gc: bool, score_self_comp: bool, backbone: str,
                         replace_5_prime: str, scoring_method: ProgramMode, genome: str, gene=None,
                         isoform=False, gene_isoforms=None) -> [Guide]:
    """ Evaluates an k-mer as a potential CRISPR target site """

    g_len = guide_size - len(pam)
    rev_comp_pam = str(Seq(pam).reverse_complement())
    dna = Seq(dna)

    if str(dna[0:2]) in allowed:
        add = True
        for pos in range(len(pam)):
            if compare_PAM(pam[pos], dna[g_len + pos]):
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
            if replace_5_prime:
                fwd = replace_5_prime + dna[len(replace_5_prime):(None if pam == "" else -len(pam))]
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
            guides = []
            id = uuid.uuid4()
            if config.isoforms:
                pam_comb = perm_PAM(pam)
                for p in pam_comb:
                    print(f"1dna[:g_len]: {dna[:g_len]}")
                    print(f"1p: {p}")
                    guides.append(generate_guide(Guide, fasta_file, p, None, name, num, guide_size, downstream_5_prime,
                                                 downstream_3_prime, '+', dna[:g_len] + p, None, score_gc,
                                                 score_self_comp, backbone, pam, replace_5_prime, scoring_method,
                                                 genome, gene, isoform, gene_isoforms, None))
            else:
                # all combinations of possible PAMs
                dna = dna.reverse_complement()
                pam_comb = perm_PAM(rev_comp_pam)
                for p in pam_comb:
                    print(f"2dna[len(rev_comp_pam):]: {dna[len(rev_comp_pam):]}")
                    print(f"2p: {p}")
                    guides.append(generate_guide(classes.Cas9.Cas9, fasta_file, p, None, name, num, guide_size,
                                                 downstream_5_prime, downstream_3_prime, '+',
                                                 p + dna[len(rev_comp_pam):], None, score_gc, score_self_comp, backbone,
                                                 pam, replace_5_prime, scoring_method, genome, gene, isoform,
                                                 gene_isoforms, None))

            return guides

    if str(dna[-2:].reverse_complement()) in allowed and not config.isoforms:
        add = True

        for pos in range(len(pam)):
            if compare_PAM(rev_comp_pam[pos], dna[pos]):
                continue
            else:
                add = False
                break

        if add and (filter_gc_min != 0 or filter_gc_max != 100):
            gc = GC(dna[len(pam):])
            if gc < filter_gc_min or gc > filter_gc_max:
                add = False

        if add and filter_self_comp_max != -1:
            if replace_5_prime:
                fwd = replace_5_prime + dna.reverse_complement()[len(pam):-len(replace_5_prime)]
            else:
                fwd = dna.reverse_complement()[len(pam):]
            folding = self_comp(fwd, backbone)
            if folding > filter_self_comp_max:
                add = False

        if add:
            guides = []
            id = uuid.uuid4()
            pam_comb = perm_PAM(rev_comp_pam)
            for p in pam_comb:
                # on the reverse strand seq of 5' downstream becomes 3' downstream and vice versa
                print(f"3dna[len(rev_comp_pam):]: {dna[len(rev_comp_pam):]}")
                print(f"3p: {p}")
                guides.append(generate_guide(classes.Cas9.Cas9, fasta_file, p, None, name, num, guide_size,
                                             Seq(downstream_3_prime).reverse_complement(),
                                             Seq(downstream_5_prime).reverse_complement(), '-',
                                             p + dna[len(rev_comp_pam):], None, score_gc, score_self_comp, backbone,
                                             pam, replace_5_prime, scoring_method, genome, gene, isoform, gene_isoforms,
                                             None))
            return guides

    return []


def eval_TALENS_sequence(name: str, guide_size: int, dna: str, num: int, fasta_file: IO, downstream_5_prime: str,
                         downstream_3_prime: str, pam: str, score_gc: bool, score_self_comp: bool, backbone: str,
                         replace_5_prime: str, scoring_method: ProgramMode, genome: str, gene=None,
                         isoform=False, gene_isoforms=None) -> [Guide]:
    """ Evaluates an N-mer as a potential TALENs target site """
    if dna[0] == "T" or dna[-1] == "A":
        return [generate_guide(classes.Cas9.Cas9, fasta_file, None, None, name, num, guide_size,
                               downstream_5_prime, downstream_3_prime, '+', dna, None, score_gc, score_self_comp,
                               backbone, pam, replace_5_prime, scoring_method, genome, gene, isoform, gene_isoforms,
                               None, talens=True)]
    return []


__all__ = ["eval_TALENS_sequence", "eval_CRISPR_sequence", "eval_CPF1_sequence", "gc_content", "compare_PAM",
           "perm_PAM"]
