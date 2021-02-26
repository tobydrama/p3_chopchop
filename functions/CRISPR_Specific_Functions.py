#####################
##
## CRISPR SPECIFIC FUNCTIONS
##

from Functions import gccontent
from Bio.SeqUtils import GC
from Bio.Seq import Seq


def eval_CRISPR_sequence(name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, allowed, PAM,
                         filterGCmin, filterGCmax, filterSelfCompMax, replace5prime=None, backbone=None):
    """ Evaluates an k-mer as a potential CRISPR target site """

    gLen = guideSize-len(PAM)
    revCompPAM = str(Seq(PAM).reverse_complement())
    dna = Seq(dna)

    if str(dna[0:2]) in allowed:
        add = True
        for pos in range(len(PAM)):
            if comaprePAM(PAM[pos], dna[gLen + pos]):
                continue
            else:
                add = False
                break

        if add and (filterGCmin != 0 or filterGCmax != 100):
            gc = GC(dna[0:(None if PAM == "" else -len(PAM))]) #FIX EVERYWHERE GC content does not assumes 5' replacement
            if gc < filterGCmin or gc > filterGCmax:
                add = False

        if add and filterSelfCompMax != -1:
            if replace5prime:
                fwd = replace5prime + dna[len(replace5prime):(None if PAM == "" else -len(PAM))]
            else:
                fwd = dna[0:(None if PAM == "" else -len(PAM))]
            folding = selfComp(fwd, backbone)
            if folding > filterSelfCompMax:
                add = False

        # in order to control the number of mismatches to search in the last 8 or 3 bps,
        # need to reverse complement so the seed region can be at the start
        # rather than end of the sequence
        # not in isoforms case as we don't search reverse complement
        if add:
            if ISOFORMS:
                pam_comb = permPAM(PAM)
                for p in pam_comb:
                    fastaFile.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                        name, num, num + guideSize, downstream5prim, downstream3prim,
                        dna, p, dna[:gLen] + p))
                return True
            else:
                # all combinations of possible PAMs
                dna = dna.reverse_complement()
                pam_comb = permPAM(revCompPAM)
                for p in pam_comb:
                    fastaFile.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                                    name, num, num+guideSize, downstream5prim, downstream3prim,
                                    dna, p, p + dna[len(revCompPAM):]))
                return True

    if str(dna[-2:].reverse_complement()) in allowed and not ISOFORMS:
        add = True

        for pos in range(len(PAM)):
            if comaprePAM(revCompPAM[pos], dna[pos]):
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
            folding = selfComp(fwd, backbone)
            if folding > filterSelfCompMax:
                add = False

        if add:
            pam_comb = permPAM(revCompPAM)
            for p in pam_comb:
                #on the reverse strand seq of 5' downstream becomes 3' downstream and vice versa
                fastaFile.write('>%s_%d-%d:%s:%s:-:%s:%s\n%s\n' % (
                                name, num, num+guideSize,
                                Seq(downstream3prim).reverse_complement(),
                                Seq(downstream5prim).reverse_complement(),
                                dna, p, p + dna[len(revCompPAM):]))
            return True

    return False


def sort_CRISPR_guides(guides):
    """ Sort pairs according to score  """
    return sorted(guides, key=attrgetter('score'))

