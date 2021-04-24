import config
from classes.ProgramMode import ProgramMode
from classes.Guide import Guide
from classes.Cas9 import Cas9
from classes.CPF1 import Cpf1
from operator import attrgetter
from functions.evaluate import eval_TALENS_sequence, eval_CPF1_sequence, eval_CRISPR_sequence


# Used in set_default_modes
def get_allowed_five_prime(allowed):
    new_allowed = []
    for el in allowed.split(","):
        if el[0] == 'N' and el[1] == 'N':
            return "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"
        elif el[0] == 'N':
            new_allowed.extend(["A"+el[1], "C"+el[1], "G"+el[1], "T"+el[1]])
        elif el[1] == 'N':
            new_allowed.extend([el[0]+"A", el[0]+"C", el[0]+"G", el[0]+"T"])
        else:
            new_allowed.append(el)
    return dict(zip(new_allowed, [True] * len(new_allowed)))


# Used in set_default_modes and tests
def get_mismatch_vectors(pam, gLength, cong):

    allowed = [True] * (gLength -len(pam))
    count = [True] * (gLength -len(pam))

    if cong:
        allowed = [True] * 9 + [False] * (gLength -len(pam) -9)

    for char in pam:
        count.append(False)
        if char == "N":
            allowed.append(True)
        else:
            allowed.append(False)

    return allowed, count


# Used in set_default_modes
def get_CPF1_mismatch_vectors(pam, gLength):

    allowed = [True] * (gLength -len(pam))
    count = [True] * (gLength -len(pam))

    for char in pam[::-1]:
        count.insert(0, False)
        if char == "N":
            allowed.insert(0,True)
        else:
            allowed.insert(0,False)

    return allowed, count


# Used in set_default_modes
def sort_CRISPR_guides(guides):
    """ Sort pairs according to score  """
    return sorted(guides, key=attrgetter('score'))


# Used in set_default_modes
def sort_TALEN_pairs(pairs):
    """ Sort pairs according to score and cluster """

    return sorted(pairs, key=attrgetter('score', 'cluster'))


# Used in main
def set_default_modes(args):
    if args.MODE == ProgramMode.CRISPR or args.MODE == ProgramMode.NICKASE:
        # Set mismatch checking policy
        (allowedMM, countMM) = get_mismatch_vectors(args.PAM, args.guideSize, args.uniqueMethod_Cong)
        allowed = get_allowed_five_prime(args.fivePrimeEnd)
        evalSequence = lambda name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim: eval_CRISPR_sequence(
            name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, allowed=allowed, PAM=args.PAM,
            filterGCmin=args.filterGCmin, filterGCmax=args.filterGCmax,
            filterSelfCompMax=args.filterSelfCompMax, replace5prime=args.replace5P, backbone=args.backbone)
        if args.MODE == ProgramMode.CRISPR:
            guideClass = Cas9 if not config.use_isoforms else Guide
            sortOutput = sort_CRISPR_guides
        elif args.MODE == ProgramMode.NICKASE:
            guideClass = Cas9
            sortOutput = sort_TALEN_pairs

    elif args.MODE == ProgramMode.CPF1:
        (allowedMM, countMM) = get_CPF1_mismatch_vectors(args.PAM, args.guideSize)
        evalSequence = lambda name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim: eval_CPF1_sequence(
            name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, PAM=args.PAM,
            filterGCmin=args.filterGCmin, filterGCmax=args.filterGCmax,
            filterSelfCompMax=args.filterSelfCompMax, replace5prime=args.replace5P, backbone=args.backbone)
        guideClass = Cpf1 if not config.use_isoforms else Guide
        sortOutput = sort_CRISPR_guides

    elif args.MODE == ProgramMode.TALENS:
        (allowedMM, countMM) = get_mismatch_vectors(args.PAM, args.guideSize, None)
        guideClass = Guide
        evalSequence = eval_TALENS_sequence
        sortOutput = sort_TALEN_pairs

    return countMM, evalSequence, guideClass, sortOutput


__all__ = ["set_default_modes"]
