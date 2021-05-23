import logging
import sys
from operator import attrgetter

import config
from classes.CPF1 import Cpf1
from classes.Cas9 import Cas9
from classes.Guide import Guide
from classes.ProgramMode import ProgramMode
from functions.evaluate import eval_talens_sequence, eval_cpf1_sequence, eval_crispr_sequence


# Used in set_default_modes
def get_allowed_five_prime(allowed):
    new_allowed = []
    for el in allowed.split(","):
        if el[0] == 'N' and el[1] == 'N':
            return "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"
        elif el[0] == 'N':
            new_allowed.extend(["A" + el[1], "C" + el[1], "G" + el[1], "T" + el[1]])
        elif el[1] == 'N':
            new_allowed.extend([el[0] + "A", el[0] + "C", el[0] + "G", el[0] + "T"])
        else:
            new_allowed.append(el)
    return dict(zip(new_allowed, [True] * len(new_allowed)))


# Used in set_default_modes and tests
def get_mismatch_vectors(pam, g_length, cong):
    allowed = [True] * (g_length - len(pam))
    count = [True] * (g_length - len(pam))

    if cong:
        allowed = [True] * 9 + [False] * (g_length - len(pam) - 9)

    for char in pam:
        count.append(False)
        if char == "N":
            allowed.append(True)
        else:
            allowed.append(False)

    return allowed, count


# Used in set_default_modes
def get_cpf1_mismatch_vectors(pam, g_length):
    allowed = [True] * (g_length - len(pam))
    count = [True] * (g_length - len(pam))

    for char in pam[::-1]:
        count.insert(0, False)
        if char == "N":
            allowed.insert(0, True)
        else:
            allowed.insert(0, False)

    return allowed, count


# Used in set_default_modes
def sort_crispr_guides(guides):
    """ Sort pairs according to score  """
    return sorted(guides, key=attrgetter('score'))


# Used in set_default_modes
def sort_talen_pairs(pairs):
    """ Sort pairs according to score and cluster """

    return sorted(pairs, key=attrgetter('score', 'cluster'))


# Used in main
def set_default_modes(args):
    if args.program_mode in [ProgramMode.CRISPR, ProgramMode.NICKASE]:
        # Set mismatch checking policy
        (allowed_mm, count_mm) = get_mismatch_vectors(args.pam, args.guide_size, args.unique_method_cong)
        allowed = get_allowed_five_prime(args.five_prime_end)
        eval_sequence = lambda name, guide_size, dna, num, fasta_file, downstream_5_prime, downstream_3_prime: \
            eval_crispr_sequence(name, guide_size, dna, num, fasta_file, downstream_5_prime, downstream_3_prime,
                                 allowed=allowed, pam=args.pam, filter_gc_min=args.filter_gc_min,
                                 filter_gc_max=args.filter_gc_max, filter_self_comp_max=args.filter_self_comp_max,
                                 replace_5prime=args.replace_5_prime, backbone=args.backbone)
        if args.program_mode == ProgramMode.CRISPR:
            guide_class = Cas9 if not config.isoforms else Guide
            sort_output = sort_crispr_guides
        elif args.program_mode == ProgramMode.NICKASE:
            guide_class = Cas9
            sort_output = sort_talen_pairs

    elif args.program_mode == ProgramMode.CPF1:
        (allowed_mm, count_mm) = get_cpf1_mismatch_vectors(args.pam, args.guide_size)
        eval_sequence = lambda name, guide_size, dna, num, fasta_file, downstream_5_prime, downstream_3_prime: \
            eval_cpf1_sequence(name, guide_size, dna, num, fasta_file, downstream_5_prime, downstream_3_prime,
                               pam=args.pam, filter_gc_min=args.filter_gc_min, filter_gc_max=args.filter_gc_max,
                               filter_self_comp_max=args.filter_self_comp_max, replace_5prime=args.replace_5_prime,
                               backbone=args.backbone)
        guide_class = Cpf1 if not config.isoforms else Guide
        sort_output = sort_crispr_guides

    elif args.program_mode == ProgramMode.TALENS:
        (allowed_mm, count_mm) = get_mismatch_vectors(args.pam, args.guide_size, None)
        guide_class = Guide
        eval_sequence = eval_talens_sequence
        sort_output = sort_talen_pairs
    else:
        logging.critical("set_default_modes: unknown program mode selected, exiting.")
        sys.exit(1)

    logging.debug(f"set_default_modes: Using guide class \"{guide_class.__name__}\", "
                  f"sorting function \"{sort_output.__name__}\".")

    return count_mm, eval_sequence, guide_class, sort_output


__all__ = ["set_default_modes", "get_allowed_five_prime", "get_mismatch_vectors", "get_cpf1_mismatch_vectors"]
