import re
import sys
import json
import numpy
import scipy.stats as ss
import warnings

from functions.Helper_Functions import *
from functions.CRISPR_Specific_Functions import eval_CRISPR_sequence, sort_CRISPR_guides
from functions.CPF1_Specific_Functions import eval_CPF1_sequence
from functions.TALEN_Specific_Functions import sort_TALEN_pairs, eval_TALENS_sequence
import classes.Cas9 as Cas9
import classes.Guide as Guide
from Vars import CONFIG, EXIT, ISOFORMS, ProgramMode, TARGET_MAX, NICKASE_DEFAULT
from Vars import DOWNSTREAM_NUC, CPF1_DEFAULT, TALEN_DEFAULT, CRISPR_DEFAULT

#Used in main
def set_default_modes(args):
    if args.MODE == ProgramMode.CRISPR or ProgramMode.NICKASE:
        # Set mismatch checking policy
        (allowedMM, countMM) = getMismatchVectors(args.PAM, args.guideSize, args.uniqueMethod_Cong)
        allowed = getAllowedFivePrime(args.fivePrimeEnd)
        evalSequence = lambda name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim: eval_CRISPR_sequence(
            name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, allowed=allowed, PAM=args.PAM,
            filterGCmin=args.filterGCmin, filterGCmax=args.filterGCmax,
            filterSelfCompMax=args.filterSelfCompMax, replace5prime=args.replace5P, backbone=args.backbone)
        if args.MODE == ProgramMode.CRISPR:
            guideClass = Cas9 if not ISOFORMS else Guide
            sortOutput = sort_CRISPR_guides
        elif args.MODE == ProgramMode.NICKASE:
            guideClass = Cas9
            sortOutput = sort_TALEN_pairs

    elif args.MODE == ProgramMode.CPF1:
        (allowedMM, countMM) = getCpf1MismatchVectors(args.PAM, args.guideSize)
        evalSequence = lambda name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim: eval_CPF1_sequence(
            name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, PAM=args.PAM,
            filterGCmin=args.filterGCmin, filterGCmax=args.filterGCmax,
            filterSelfCompMax=args.filterSelfCompMax, replace5prime=args.replace5P, backbone=args.backbone)
        guideClass = ProgramMode.Cpf1 if not ISOFORMS else Guide
        sortOutput = sort_CRISPR_guides

    elif args.MODE == ProgramMode.TALENS:
        (allowedMM, countMM) = getMismatchVectors(args.PAM, args.guideSize, None)
        guideClass = Guide
        evalSequence = eval_TALENS_sequence
        sortOutput = sort_TALEN_pairs

    return countMM, evalSequence, guideClass, sortOutput


#Used in main
def scoreChari_2015(svmInputFile, svmOutputFile, PAM, genome):  #Only one use in main
    f_p = sys.path[0]
    """ Calculate score from SVM model as in Chari 2015 20-NGG or 20-NNAGAAW, only for hg19 and mm10"""

    model = f_p + '/models/293T_HiSeq_SP_Nuclease_100_SVM_Model.txt'
    dist = f_p + '/models/Hg19_RefFlat_Genes_75bp_NoUTRs_SPSites_SVMOutput.txt'

    if PAM == 'NGG' and genome == 'mm10':
        model = f_p + '/models/293T_HiSeq_SP_Nuclease_100_SVM_Model.txt'
        dist = f_p + '/models/Mm10_RefFlat_Genes_75bp_NoUTRs_SPSites_SVMOutput.txt'
    elif PAM == 'NNAGAAW' and genome == 'hg19':
        model = f_p + '/models/293T_HiSeq_ST1_Nuclease_100_V2_SVM_Model.txt'
        dist = f_p + '/models/Hg19_RefFlat_Genes_75bp_NoUTRs_ST1Sites_SVMOutput.txt'
    elif PAM == 'NNAGAAW' and genome == 'mm10':
        model = f_p + '/models/293T_HiSeq_ST1_Nuclease_100_V2_SVM_Model.txt'
        dist = f_p + '/models/Mm10_RefFlat_Genes_75bp_NoUTRs_ST1Sites_SVMOutput.txt'

    prog = Popen("%s/svm_light/svm_classify -v 0 %s %s %s" % (f_p, svmInputFile, model, svmOutputFile), shell=True)
    prog.communicate()

    svmAll = open(dist,'r')
    svmThis = open(svmOutputFile, 'r')

    # first through go all scores and get the max and min
    allData = []
    for line in svmAll:
        line = line.rstrip('\r\n')
        allData.append(float(line))
    svmAll.close()

    scoreArray = []
    for line in svmThis:
        line = line.rstrip('\r\n')
        scoreArray.append(float(line))

    return [ss.percentileofscore(allData, i) for i in scoreArray]

#Used in main
def concatenate_feature_sets(feature_sets):
    '''
    Given a dictionary of sets of features, each in a Pandas.DataFrame,
    concatenate them together to form one big np.array, and get the dimension
    of each set
    Returns: inputs, dim
    Source: Doench 2016
    '''
    assert feature_sets != {}, "no feature sets present"
    F = feature_sets[feature_sets.keys()[0]].shape[0]
    for fset in feature_sets.keys():
        F2 = feature_sets[fset].shape[0]
        assert F == F2, "not same # individuals for features %s and %s" % (feature_sets.keys()[0], fset)

    N = feature_sets[feature_sets.keys()[0]].shape[0]
    inputs = numpy.zeros((N, 0))
    feature_names = []
    dim = {}
    dimsum = 0
    for fset in feature_sets.keys():
        inputs_set = feature_sets[fset].values
        dim[fset] = inputs_set.shape[1]
        dimsum += dim[fset]
        inputs = numpy.hstack((inputs, inputs_set))
        feature_names.extend(feature_sets[fset].columns.tolist())

    return inputs, dim, dimsum, feature_names

#Used in main
def coordToFasta(regions, fasta_file, outputDir, targetSize, evalAndPrintFunc, nonOver, indexDir, genome, strand, ext):
    """ Extracts the sequence corresponding to genomic coordinates from a FASTA file """

    ext = 0 if ISOFORMS else ext # for genomic context for some models
    sequences = {}
    fasta_file = open(fasta_file, 'w')
    fasta_seq = ""

    if ISOFORMS and strand == "-":
        regions = regions[::-1]

    for region in regions:
        # Extracts chromosome number and region start and end
        chrom = region[0:region.rfind(':')]
        start = int(region[region.rfind(':')+1:region.rfind('-')])
        finish = int(region[region.rfind('-')+1:])
        start = max(start, 0)

        if ext == 0 and finish == start:
            continue

        # Run twoBitToFa program to get actual dna sequence corresponding to input genomic coordinates
        # Popen runs twoBitToFa program. PIPE pipes stdout.
        prog = Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2> %s/twoBitToFa.err" % (
            CONFIG["PATH"]["TWOBITTOFA"], chrom, start - ext, finish + ext, indexDir, genome, outputDir), stdout=PIPE, shell=True)

        # Communicate converts stdout to a string
        output = prog.communicate()
        if prog.returncode != 0:
            sys.stderr.write("Running twoBitToFa failed\n")
            sys.exit(EXIT['TWOBITTOFA_ERROR'])

        output = output[0]
        exons = output.split("\n")
        dna = ''.join(exons[1:]).upper()
        ext_dna = dna
        dna = dna[ext:(len(dna)-ext)]
        if len(dna) != (finish - start):  # something is wrong with what was fetched by twoBitToFa
            continue

        if ISOFORMS and strand == "-":
            dna = str(Seq(dna).reverse_complement())

        # Write exon sequences to text file user can open in ApE. exon-intron junctions in lowercase.
        fasta_seq += dna[0].lower()+dna[1:-1]+dna[-1].lower()

        # Add 1 due to BED 0-indexing
        name = "C:%s:%d-%d" % (chrom, start, finish)

        # Loop over exon sequence, write every g-mer into file in which g-mer ends in PAM in fasta format
        positions = range(0, len(dna)-(targetSize-1))
        while len(positions) != 0:
            num = positions.pop(0)
            downstream_5prim = ext_dna[num:(num + ext)]
            g_end = num + ext + targetSize
            downstream_3prim = ext_dna[g_end:(g_end + ext)]
            if evalAndPrintFunc(name, targetSize, dna[num:(num + targetSize)],
                                len(dna) - num - targetSize if ISOFORMS and strand == "-" else num, fasta_file,
                                downstream_5prim, downstream_3prim):
                if nonOver:  # positions overlapping those of this guide
                    for p in range(num, num + targetSize):
                        if p in positions:
                            positions.remove(p)

                if name not in sequences:
                    sequences[name] = dna

    fasta_file.close()

    if ISOFORMS and strand == "-":
        fasta_seq = str(Seq(fasta_seq).reverse_complement())

    return sequences, fasta_seq

#Used in main
def runBowtie(PAMlength, unique_method_cong, fasta_file, output_dir,
              max_off_targets, index_dir, genome, max_mismatches):

    bwt_results_file = '%s/output.sam' % output_dir
    if unique_method_cong and not ISOFORMS:
        # When ISOFORMS dna string is not reverse complemented and Cong can't be used
        # the -l alignment mode specifies a seed region to search for the number of mismatches specified with the
        # -n option. Outside of that seed, up to 2 mismatches are searched.
        # E.g. -l 15 -n 0 will search the first 15 bases with no mismatches, and the rest with up to 3 mismatches
        command = "%s -p %s -l %d -n %d -m %d --sam-nohead -k %d %s/%s -f %s -S %s " % (
            CONFIG["PATH"]["BOWTIE"], CONFIG["THREADS"], (PAMlength + 11), max_mismatches, max_off_targets, max_off_targets, index_dir,
            genome, fasta_file, bwt_results_file)
    else:
        command = "%s -p %s -v %d --sam-nohead -k %d %s/%s -f %s -S %s " % (
            CONFIG["PATH"]["BOWTIE"], CONFIG["THREADS"], max_mismatches, max_off_targets, index_dir, genome, fasta_file, bwt_results_file)

    if ISOFORMS: # When ISFORMS we don't check reverse complement
        command += "--norc "

    command += "2> %s/bowtie.err" % output_dir

    prog = Popen(command, shell=True)
    prog.wait()

    if prog.returncode != 0:
        sys.stderr.write("Running bowtie failed\n")
        sys.exit(EXIT['BOWTIE_ERROR'])

    return bwt_results_file

#Used in main
def make_primers_fasta(targets, outputDir, flanks, displayFlanks, genome, limitPrintResults, bowtieIndexDir,
                       fastaSequence, primer3options, guidePadding, enzymeCo, minResSiteLen, geneID, maxOffTargets):
    primers = {}
    primerOpt = get_primer_options(primer3options)

    primerFastaFileName = '%s/primers.fa' % outputDir
    primerFastaFile = open(primerFastaFileName, 'w')
    for i in range(min(limitPrintResults-1, len(targets))):
        target = targets[i]
        seq, seqLenBeforeTarget = get_primer_query_sequence_fasta(target, outputDir, flanks, fastaSequence)
        primer3_output = make_primer_for_target(target, outputDir, seq, seqLenBeforeTarget, primerOpt, guidePadding)
        region = "%s:%s-%s" % (target.chrom, max(0, target.start-flanks), min(len(fastaSequence), target.end+flanks))
        target_primers, primerPos = parse_primer3_output(target, region, primer3_output, primerFastaFile)
        primers[target.ID] = target_primers

        # Restriction sites
        restSites = dump_restriction_sites(target, seq, flanks, enzymeCo, outputDir, minResSiteLen)
        # Sequence for visualization of locus
        seq2, seqLenBeforeTarget2 = get_primer_query_sequence_fasta(target, outputDir, displayFlanks, fastaSequence)
        dump_locus_sequence(target, outputDir, seq2, seqLenBeforeTarget2, "+")
        # Genbank file for download
        dump_genbank_file(seq, target, restSites, primerPos, outputDir, geneID, target.start-seqLenBeforeTarget, "+")

    primerFastaFile.close()

    primerResults = runBowtiePrimers(primerFastaFileName, outputDir, genome, bowtieIndexDir, maxOffTargets)
    pairPrimers(primers, primerResults, outputDir)

#Used in main
def make_primers_genome(targets, outputDir, flanks, display_seq_len, genome, limitPrintResults, bowtieIndexDir, twoBitToFaIndexDir,
                        primer3options, guidePadding, enzymeCo, minResSiteLen, strand, geneID, maxOffTargets):
    primers = {}

    primerOpt = get_primer_options(primer3options)

    # RUN PRIMER3 ON TARGET SITES AND CREATE FASTA FILE OF PRIMERS FOR BOWTIE
    primerFastaFileName = '%s/primers.fa' % outputDir
    primerFastaFile = open(primerFastaFileName, 'w')
    for i in range(min(limitPrintResults-1, len(targets))):
        target = targets[i]
        seq, seqLenBeforeTarget = get_primer_query_sequence_2bit(
            target, outputDir, flanks, genome, twoBitToFaIndexDir, strand)
        primer3_output = make_primer_for_target(target, outputDir, seq, seqLenBeforeTarget, primerOpt, guidePadding)
        region = "%s:%s-%s" % (target.chrom, max(0, target.start-flanks), target.end+flanks)
        target_primers, primerPos = parse_primer3_output(target, region, primer3_output, primerFastaFile)
        primers[target.ID] = target_primers

        # Restriction sites
        restSites = dump_restriction_sites(target, seq, flanks, enzymeCo, outputDir, minResSiteLen)
        # Sequence for visualization of locus
        seq2, seqLenBeforeTarget2 = get_primer_query_sequence_2bit(
            target, outputDir, display_seq_len, genome, twoBitToFaIndexDir, strand)
        dump_locus_sequence(target, outputDir, seq2, seqLenBeforeTarget2, strand)
        # Genbank file for download
        dump_genbank_file(seq, target, restSites, primerPos, outputDir, geneID, target.start-seqLenBeforeTarget, strand)

    primerFastaFile.close()

    primerResults = runBowtiePrimers(primerFastaFileName, outputDir, genome, bowtieIndexDir, maxOffTargets)
    pairPrimers(primers, primerResults, outputDir)

#Used in main
def writeIndividualResults(outputDir, maxOffTargets, sortedOutput, guideSize, mode, totalClusters, limitPrintResults, offtargetsTable):
    """ Writes each guide and its offtargets into a file """

    # Initiate list of lists for each cluster
    clusters = [[] for i in range(totalClusters)]

    fileHandler = dict()

    # Limit the number of open files (and results)
    sortedOutput = sortedOutput[0:min(len(sortedOutput), limitPrintResults-1)]

    for i in range(len(sortedOutput)):
        current = sortedOutput[i]
        current.ID = i+1

        # Create new file if not already opened
        if current.ID not in fileHandler:
            resultsFile = '%s/%s.offtargets' % (outputDir, current.ID)
            fileHandler[current.ID] = open(resultsFile, 'w')
        f = fileHandler[current.ID]

        # Add the current TALE pair to the appropriate list in the list of lists, depending on its cluster number
        if mode == TALENS or mode == NICKASE:
            clusterID = current.cluster
            clusters[clusterID-1].append(current)

        offTargets = current.asOffTargetString("", maxOffTargets)
        if not offTargets:
            offTargets = "There are no predicted off-targets."

        f.write(str(current.strandedGuideSeq)+"\n"+offTargets+"\n")

        if mode == CRISPR and not ISOFORMS and current.repStats is not None:
            stats_file = '%s/%s_repStats.json' % (outputDir, current.ID)
            with open(stats_file, 'w') as fp:
                json.dump(current.repStats, fp)
            fp.close()

        if mode == CRISPR and not ISOFORMS and current.repProfile is not None:
            profile_file = '%s/%s_repProfile.csv' % (outputDir, current.ID)
            current.repProfile.to_csv(profile_file, index=False)

        if mode == CRISPR and not ISOFORMS and offtargetsTable:
            off_table = '%s/offtargetsTable.csv' % outputDir
            label = "%s:%s,%s,%s" % (current.chrom, current.start, current.strand, current.strandedGuideSeq)
            off_for_table = map(lambda x: x.asOffTargetString(label, maxOffTargets), current.offTargets)
            with open(off_table, "a") as append_file:
                if len(off_for_table) > 0:
                    append_file.write("\n".join(off_for_table))
                    append_file.write("\n")

    for clust in clusters:
        if len(clust) == 0:
            continue
        bestInCluster = clust[0]

        for member in clust[1:]:
            # Write the other cluster members to file
            fileHandler[bestInCluster.ID].write("%s*%s*%s,%s:%s,%s,%s/%s,%s/%s,%s/%s,%s/%s;" % (
                member.tale1.guideSeq, member.spacerSeq, member.tale2.guideSeq, member.chrom, member.start,
                len(member.offTargetPairs), member.tale1.offTargetsMM[0], member.tale2.offTargetsMM[0],
                member.tale1.offTargetsMM[1], member.tale2.offTargetsMM[1], member.tale1.offTargetsMM[2],
                member.tale2.offTargetsMM[2], member.tale1.offTargetsMM[3], member.tale2.offTargetsMM[3]))

        fileHandler[bestInCluster.ID].write("\n"+current.restrictionSites+"\n")

    for fh in fileHandler.values():
        fh.close()

    return clusters

#Used in main
def getAllowedFivePrime(allowed):
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

#Used in main
def parseTargets(target_string, genome, use_db, data, pad_size, target_region, exon_subset, ups_bp, down_bp,
                 index_dir, output_dir, use_union, make_vis, guideLen):
    targets = []
    vis_coords = []
    target_strand = "+"
    target_size = 0
    gene, isoform, gene_isoforms = (None, None, set())

    pattern = re.compile("(([\.\w]+):)?([\.\,\d]+)\-([\.\,\d]+)")
    is_coordinate = pattern.match(str(target_string))

    if is_coordinate:
        if ISOFORMS:
            sys.stderr.write("--isoforms is not working with coordinate search.\n")
            sys.exit(EXIT['ISOFORMS_ERROR'])

        chrom = is_coordinate.group(2)
        vis_coords.append({"exons": [], "ATG": [], "name": chrom})

        for target in target_string.split(";"):
            m = pattern.match(target)
            if m:
                if m.group(2) is not None and chrom != m.group(2):
                    sys.stderr.write(
                        "Can't target regions on separate chromosomes (%s != %s).\n" % (chrom, m.group(2)))
                    sys.exit(EXIT['GENE_ERROR'])

                start_pos = m.group(3)
                end_pos = m.group(4)
                start_pos = int(start_pos.replace(",", "").replace(".", ""))
                end_pos = int(end_pos.replace(",", "").replace(".", ""))
                target_size += end_pos - start_pos + 1

                if start_pos >= end_pos:
                    sys.stderr.write(
                        "Start position (%s) must be smaller than end position (%s)\n" % (start_pos, end_pos))
                    sys.exit(EXIT['GENE_ERROR'])

                targets.append("%s:%s-%s" % (chrom, max(0, start_pos - pad_size), end_pos + pad_size))
                if make_vis:
                    vis_coords[0]["exons"].append([chrom, start_pos, end_pos, 0, True, "+"])
            else:
                sys.stderr.write("Unknown format: %s\n" % (target))
                sys.exit(EXIT['GENE_ERROR'])

    else:
        if use_db:
            if ISOFORMS:
                sys.stderr.write("--isoforms is not working with database search.\n")
                sys.exit(EXIT['ISOFORMS_ERROR'])
            txInfo = geneToCoord_db(target_string, genome, data)
            txInfo = filterRepeatingNames(txInfo)
        else:
            gene, txInfo = geneToCoord_file(target_string, data)
            txInfo = filterRepeatingNames(txInfo)
            isoform = "union" if use_union else "intersection"
            gene_isoforms = set([str(x[3]) for x in txInfo])
            if target_string in gene_isoforms:
                isoform = target_string
                gene_isoforms = get_isoforms(gene, data)

        target_chr = set([x[0] for x in txInfo])
        target_strand = set([x[6] for x in txInfo])
        isoforms = [str(x[3]) for x in txInfo]
        if len(target_strand) > 1 or len(target_chr) > 1:
            sys.stderr.write(
                "Specify which isoform you want to target as your query " + str(target_string) +
                " returns many isoforms: " + ', '.join(isoforms) +
                " which are from either inconsistent strands or chromosomes.\n")
            sys.exit(EXIT['GENE_ERROR'])
        else:
            target_strand = list(target_strand)[0]
            target_chr = list(target_chr)[0]

        for tx in txInfo:
            tx = list(tx)
            tx[4] = int(tx[4])
            tx[5] = int(tx[5])
            starts = tx[1].split(",")
            ends = tx[2].split(",")
            del starts[-1]
            del ends[-1]
            starts = map(int, starts)
            ends = map(int, ends)
            starts_v = starts[:]
            ends_v = ends[:]
            tx_vis = {"exons": [], "ATG": [], "name": tx[3]}

            if make_vis:
                intron_size = [int(starts_v[x + 1]) - int(ends_v[x]) for x in range(len(starts_v) - 1)]
                intron_size.append(0)
                # tx_vis exons are [chr, start, end, intron_size, isIntron, strand]
                for e in range(len(starts_v)):
                    if ends_v[e] <= tx[4] or starts_v[e] >= tx[5]:
                        tx_vis["exons"].append([tx[0], starts_v[e], ends_v[e], intron_size[e], True, tx[6]])
                    else:
                        if starts_v[e] < tx[4] < ends_v[e]:
                            tx_vis["exons"].append([tx[0], starts_v[e], tx[4], 0, True, tx[6]])
                            starts_v[e] = tx[4]

                        if starts_v[e] < tx[5] < ends_v[e]:
                            tx_vis["exons"].append([tx[0], tx[5], ends_v[e], intron_size[e], True, tx[6]])
                            ends_v[e] = tx[5]
                            intron_size[e] = 0

                        tx_vis["exons"].append([tx[0], starts_v[e], ends_v[e], intron_size[e], False, tx[6]])

                tx_vis["exons"].sort(key=lambda x: x[1]) # sort on starts
                # ATG locations
                prog = Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2> %s/twoBitToFa.err" % (
                    CONFIG["PATH"]["TWOBITTOFA"], tx[0], int(tx[4]) + 1, int(tx[5]) + 1, index_dir,
                    genome, output_dir), stdout=PIPE, shell=True)
                iso_seq = prog.communicate()
                if prog.returncode != 0:
                    sys.stderr.write("Running twoBitToFa when searching isoform sequence failed\n")
                    sys.exit(EXIT['TWOBITTOFA_ERROR'])

                iso_seq = iso_seq[0]
                iso_seq = iso_seq.split("\n")
                iso_seq = Seq(''.join(iso_seq[1:]).upper())
                # splicing
                iso_seq_spl = ""
                for e in tx_vis["exons"]:
                    if not e[4]:
                        iso_seq_spl += iso_seq[(e[1] - tx[4]):(e[2] - tx[4])]
                atg = "ATG" if tx[6] != "-" else "CAT"
                tx_atg = [m.start() for m in re.finditer(atg, str(iso_seq_spl)) if m.start() % 3 == 0]
                tx_atg.sort()
                for atg1 in tx_atg: # every ATG as 3 x 1bp as they can span across two exons...
                    atg2 = atg1 + 1
                    atg3 = atg1 + 2
                    shift_atg1, shift_atg2, shift_atg3, exon_len = 0, 0, 0, 0
                    for e in tx_vis["exons"]: # exons are sorted
                        if not e[4]:
                            exon_len += (e[2] - e[1])
                            if atg1 > exon_len:
                                shift_atg1 += e[3]
                            if atg2 > exon_len:
                                shift_atg2 += e[3]
                            if atg3 > exon_len:
                                shift_atg3 += e[3]
                    tx_vis["ATG"].extend([atg1 + shift_atg1 + tx[4], atg2 + shift_atg2 + tx[4],
                                          atg3 + shift_atg3 + tx[4]])

                vis_coords.append(tx_vis)

            # restrict isoforms
            coords = map(lambda x: [tx[0], x[0], x[1]], zip(starts, ends))
            if tx[6] == "-":
                coords.reverse()
            coords = subsetExons(exon_subset, coords)
            if tx[6] == "-":
                coords.reverse()

            # Truncate to region
            if target_region == "CODING":
                coords = truncateToCoding(tx[4], tx[5], coords)
            elif target_region == "UTR5":
                if tx[6] == "+":
                    coords = truncateToUTR5(tx[4], coords)
                else:
                    coords = truncateToUTR3(tx[5], coords)
            elif target_region == "PROMOTER":
                coords = truncateToPROMOTER(tx[6], coords, ups_bp, down_bp)
            elif target_region == "UTR3":
                if tx[6] == "+":
                    coords = truncateToUTR3(tx[5], coords)
                else:
                    coords = truncateToUTR5(tx[4], coords)
            elif target_region == "SPLICE":
                coords = truncateToSplice(coords)
            elif target_region != "WHOLE":
                sys.stderr.write("Unknown region: %s\n" % target_region)
                sys.exit(EXIT['PYTHON_ERROR'])

            # filter exons that are too truncated
            coords = [x for x in coords if x[1] < x[2]]
            if not coords:
                if gene_isoforms:
                    gene_isoforms.remove(tx[3])
                if vis_coords:
                    del vis_coords[-1]

            # compute intersection/union on all exons
            if txInfo[0][3] == tx[3]:  # if this is first of the isoforms
                for x in coords:
                    targets.extend(range(x[1], x[2] + 1))
                targets = set(targets)
            else:
                if not use_union:
                    targets_ = []
                    for x in coords:
                        targets_.extend(range(x[1], x[2] + 1))

                    if len(targets_) >= guideLen: # cover cases where some transcripts provide short or none bp
                        targets &= set(targets_)

                    if len(targets) < guideLen:
                        sys.stderr.write(
                            "Computing intersection over specified isoforms resulted in lack of targets." +
                            " Consider either using specific isoform as input: " + ', '.join(isoforms) +
                            " or using --consensusUnion to compute union instead of intersection " +
                            "of your isoforms (on the website you can find it in " +
                            "Options -> General -> Isoform consensus determined by -> Union.")
                        sys.exit(EXIT['GENE_ERROR'])
                else:
                    targets_ = []
                    for x in coords:
                        targets_.extend(range(x[1], x[2] + 1))
                    targets |= set(targets_)

        target_size = len(targets)
        if target_size < guideLen:
            sys.stderr.write("Search region is too small. You probably want to specify -t option as WHOLE")
            sys.exit(EXIT['GENE_ERROR'])

        starts, ends = bins(targets)
        if ISOFORMS:
            targets = map(lambda x: "%s:%s-%s" % (target_chr, x[0], x[1]), zip(starts, ends))
        else:
            targets = map(lambda x: "%s:%s-%s" % (target_chr, x[0] - pad_size, x[1] + pad_size), zip(starts, ends))

    if target_size > TARGET_MAX:
        sys.stderr.write("Search region is too large (%s nt). Maximum search region is %s nt.\n" % (
            target_size, TARGET_MAX))
        sys.exit(EXIT['GENE_ERROR'])

    return targets, vis_coords, target_strand, gene, isoform, gene_isoforms

#Used in main
def parseFastaTarget(fasta_file, candidate_fasta_file, target_size, eval_and_print):
    """ Parse a FASTA file as input for targeting """

    fasta_file = list(SeqIO.parse(fasta_file, 'fasta'))
    seq_name, sequence = fasta_file[0].id, str(fasta_file[0].seq)

    name = "%s:0-%s" % (seq_name, len(sequence))
    id_name = "C:" + name
    sequence = sequence.upper()
    sequence = "".join(sequence.split())

    dna_pattern = re.compile(r'([^ACGTNacgtn])')
    if dna_pattern.search(sequence):
        sys.stderr.write("Input sequence contains illegal characters.\n")
        sys.exit(EXIT['GENE_ERROR'])

    sequences = {}
    candidate_fasta_file = open(candidate_fasta_file, 'w')

    # Loop over sequence, write every k-mer into file in which k-mer ends in as PAM in fasta format
    for num in range(0, len(sequence)-(target_size-1)):

        if (num - DOWNSTREAM_NUC) > 0:
                start5prim = num - DOWNSTREAM_NUC
        else:
                start5prim = 0

        if (num + target_size + DOWNSTREAM_NUC) > len(sequence):
                end3prim = len(sequence)
        else:
                end3prim = num + target_size + DOWNSTREAM_NUC

        downstream_5prim = sequence[start5prim:num]
        downstream_3prim = sequence[(num + target_size):end3prim]

        if eval_and_print(id_name, target_size, sequence[num:(num + target_size)], num,
                          candidate_fasta_file, downstream_5prim, downstream_3prim):
            sequences[id_name] = sequence

    return sequences, [name], [{"exons": [[seq_name, 1, len(sequence), 0, 20, "+"]],
                                "ATG": [], "name": seq_name}], sequence, "+"

#Used in main
def connect_db(database_string):
    import MySQLdb

    m = re.compile("(.+):(.+)@(.+)/(.+)").search(database_string)
    if not m:
        sys.stderr.write("Wrong syntax for connection string: username:pass@localhost/db_name")
        sys.exit(EXIT["DB_ERROR"])

    try:
        db = MySQLdb.connect(user = m.group(1), passwd = m.group(2), host = m.group(3), db = m.group(4))
    except:
        sys.stderr.write("Could not connect to database\n")
        sys.exit(EXIT['DB_ERROR'])

    return db

#Used in main and tests
def getMismatchVectors(pam, gLength, cong):

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

#Used in main
def getCpf1MismatchVectors(pam, gLength):

    allowed = [True] * (gLength -len(pam))
    count = [True] * (gLength -len(pam))

    for char in pam[::-1]:
        count.insert(0, False)
        if char == "N":
            allowed.insert(0,True)
        else:
            allowed.insert(0,False)

    return allowed, count

#Used in main
def mode_select(var, index, MODE):
    """ Selects a default depending on mode for options that have not been set """

    if var is not None:
        return var

    if MODE == CRISPR:
        return CRISPR_DEFAULT[index]

    elif MODE == TALENS:
        return TALEN_DEFAULT[index]

    elif MODE == CPF1:
        return CPF1_DEFAULT[index]

    elif MODE == NICKASE:
        return NICKASE_DEFAULT[index]

    sys.stderr.write("Unknown model %s\n" % MODE)
    sys.exit(EXIT['PYTHON_ERROR'])

#Used in main
def print_bed(mode, vis_cords, targets, output_file, description): # bed is 0-based
    bed_file = open(output_file, 'w')

    if mode == CRISPR:
        thresholds = [0, 1000]
    elif mode == CPF1:
        thresholds = [300, 1000]
    elif mode == NICKASE:
        thresholds = [3000, 6000]
    else:
        thresholds = [10000, 15000]

    if targets is not None:

        chromosome = vis_cords[0]["exons"][0][0]
        min_loc = min([x["exons"][0][1] for x in vis_cords])
        max_loc = max([x["exons"][-1][2] for x in vis_cords])

        header = """track name=CHOPCHOP description=""" + description + """ visibility="pack" itemRgb="On"\n"""
        bed_file.write("browser position {0}:{1}-{2}\n".format(chromosome, min_loc, max_loc))
        bed_file.write(header)

        for target in targets:

            color = "0,128,0"  # green
            if target[2] >= thresholds[0]:
                color = "255,255,0"  # yellow
            if target[2] >= thresholds[1]:
                color = "255,0,0"  # red

            if mode == CRISPR or mode == CPF1:
                start = target[1] - 1
                stop = target[1] + target[3] - 1
            else:
                start = target[6] - 1
                stop = target[7] - 1

            bed_line = "{0}\t{1}\t{2}\tRanked:{3}\t{4}\t{5}\t{1}\t{2}\t{6}\n".format(chromosome, start, stop,
                                                                                     target[0], 0, target[4], color)
            bed_file.write(bed_line)

    bed_file.close()

#Used in main
def print_genbank(mode, name, seq, exons, targets, chrom, seq_start, seq_end, strand, output_file, description): # different than other dump_gb
    genbank_file = open(output_file, 'w')
    loci = chrom + ":" + str(seq_start) + "-" + str(seq_end)
    if len(name) > 10: # almost always... Genbank seems a bit outdated as format
        name = name[-10:]
    if len(loci) > 10: # almost always...
        loci = name[-10:]
    record = SeqRecord(Seq(seq, IUPACAmbiguousDNA()), description=description,
                       name=name, id=loci)
    gene_strand = 1 if strand == "+" else -1
    # genbank is 0-based
    if len(targets) > 0:
        for target in targets:
            ts = 1 if target[4] == "+" else -1
            if ISOFORMS:
                ts = gene_strand

            if mode == CRISPR or mode == CPF1:
                start = target[1] - 1
                stop = target[1] + target[3] - 1
            else:
                start = target[6] - 1
                stop = target[7] - 1

            record.features.append(SeqFeature(FeatureLocation(start-seq_start, stop-seq_start,
                                                              strand=ts), type="Target_%s" % target[0]))

    if len(exons) > 0:
        for exon in exons:
            record.features.append(SeqFeature(FeatureLocation(exon[1]-seq_start, exon[2]-seq_start,
                                                              strand=gene_strand), type="gene_loci"))

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("ignore")
        SeqIO.write(record, genbank_file, "genbank")
    genbank_file.close()

#Used in main
def rna_folding_metric(specie, tx_id, tx_start, tx_end):
    mean_bpp = 0
    file_path = CONFIG["PATH"]["ISOFORMS_MT_DIR"] + "/" + specie + "/" + tx_id + ".mt"
    if os.path.isfile(file_path):
        mt = pandas.read_csv(file_path, sep="\t", header=None, skiprows=tx_start, nrows=tx_end - tx_start)
        mean_bpp = numpy.mean(mt[1].tolist())

    return mean_bpp

#used in main
def tx_relative_coordinates(visCoords, tx_id, start, end):
    tx_start, tx_end = -1, -1
    exons = [e["exons"] for e in visCoords if e["name"] == tx_id][0]
    e_id = -1
    for i, e in enumerate(exons):
        if e[1] <= (start - 1) and e[2] >= (end - 1):
            e_id = i
            break

    if e_id is not -1:
        for i in range(0, e_id) if exons[0][5] == "+" else range(e_id + 1, len(exons)):
            tx_start += exons[i][2] - exons[i][1]

        tx_start += (exons[e_id][1] - start - 1) if exons[0][5] == "+" else (exons[e_id][2] - end - 1)
        tx_end = tx_start + end - start

    return tx_start, tx_end

#
def FastaToViscoords(sequences, strand):
    """ Makes the exons in 'sequences' array generated in coordToFasta json readable for visualization"""
    exonstart = []
    exonend = []
    exonsequence = []

    for exon in sequences:
        # sys.stderr.write("%s\n" % exon)
        exonlist = exon.split(':')
        exoncoord = exonlist[2].split('-')
        exonstart.append(exoncoord[0])
        exonend.append(exoncoord[1])
        seq = sequences[exon]
        if strand == "-":
            seq = complement(seq)

        exonsequence.append(seq)

    return zip(exonstart, exonend, exonsequence)

#Used in main
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
