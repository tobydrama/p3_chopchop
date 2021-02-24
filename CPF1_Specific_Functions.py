#####################
##
## CPF1 SPECIFIC FUNCTIONS
##
from Bio.Seq import Seq
from Functions import comaprePAM

def eval_CPF1_sequence(name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, PAM,
    filterGCmin, filterGCmax, filterSelfCompMax, replace5prime = None, backbone = None):
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
        if comaprePAM(PAM[pos], dna[pos]):
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
        folding = selfComp(fwd, backbone)
        if folding > filterSelfCompMax:
            add = False

    if add:
        if ISOFORMS:
            pam_comb = permPAM(PAM)
            for p in pam_comb:
                fastaFile.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                    name, num, num + guideSize, downstream5prim, downstream3prim,
                    dna, p, p + dna[len(PAM):]))
        else:
            dna = dna.reverse_complement()
            pam_comb = permPAM(revCompPAM)
            for p in pam_comb:
                fastaFile.write('>%s_%d-%d:%s:%s:+:%s:%s\n%s\n' % (
                                name, num, num+guideSize, downstream5prim, downstream3prim,
                                dna, p, dna[:gLen] + p))
        return True


    add = True and not ISOFORMS

    for pos in range(len(PAM)):
        if comaprePAM(revCompPAM[pos], dna[gLen + pos]):
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
                            dna, p, dna[:gLen] + p))
        return True

    return False
