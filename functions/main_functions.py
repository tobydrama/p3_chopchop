#####################
##
## FUNCTIONS ONLY USED IN MAIN
##
from typing import Union, List

import scipy.stats as ss
import logging
import warnings
import re
import sys
import json

import config
from classes.ProgramMode import ProgramMode
from Vars import CONFIG, EXIT, NICKASE_DEFAULT, PRIMER_OFF_TARGET_MIN
from Vars import DOWNSTREAM_NUC, CPF1_DEFAULT, TALEN_DEFAULT, CRISPR_DEFAULT
from subprocess import Popen, PIPE
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


# Here lies concatenate_feature_sets

# Used in main
def coord_to_fasta(regions, fasta_file, outputDir, targetSize, evalAndPrintFunc, nonOver, indexDir, genome, strand, ext):
    """ Extracts the sequence corresponding to genomic coordinates from a FASTA file """

    ext = 0 if config.isoforms else ext # for genomic context for some models
    sequences = {}
    fasta_file = open(fasta_file, 'w')
    fasta_seq = ""

    if config.isoforms and strand == "-":
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

        output = output[0].decode()
        exons = output.split("\n")
        dna = ''.join(exons[1:]).upper()
        ext_dna = dna
        dna = dna[ext:(len(dna)-ext)]
        if len(dna) != (finish - start):  # something is wrong with what was fetched by twoBitToFa
            continue

        if config.isoforms and strand == "-":
            dna = str(Seq(dna).reverse_complement())

        # Write exon sequences to text file user can open in ApE. exon-intron junctions in lowercase.
        fasta_seq += dna[0].lower()+dna[1:-1]+dna[-1].lower()

        # Add 1 due to BED 0-indexing
        name = "C:%s:%d-%d" % (chrom, start, finish)

        # Loop over exon sequence, write every g-mer into file in which g-mer ends in PAM in fasta format
        positions = list(range(0, len(dna)-(targetSize-1)))
        while len(positions) != 0:
            num = positions.pop(0)
            downstream_5prim = ext_dna[num:(num + ext)]
            g_end = num + ext + targetSize
            downstream_3prim = ext_dna[g_end:(g_end + ext)]
            if evalAndPrintFunc(name, targetSize, dna[num:(num + targetSize)],
                                len(dna) - num - targetSize if config.isoforms and strand == "-" else num, fasta_file,
                                downstream_5prim, downstream_3prim):
                if nonOver:  # positions overlapping those of this guide
                    for p in range(num, num + targetSize):
                        if p in positions:
                            positions.remove(p)

                if name not in sequences:
                    sequences[name] = dna

    fasta_file.close()

    if config.isoforms and strand == "-":
        fasta_seq = str(Seq(fasta_seq).reverse_complement())

    return sequences, fasta_seq


# Used in main
def run_bowtie(PAMlength, unique_method_cong, fasta_file, output_dir,
               max_off_targets, index_dir, genome, max_mismatches):
    logging.info("Running bowtie.")

    bwt_results_file = '%s/output.sam' % output_dir
    if unique_method_cong and not config.isoforms:
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

    if config.isoforms: # When ISFORMS we don't check reverse complement
        command += "--norc "

    command += "2> %s/bowtie.err" % output_dir

    prog = Popen(command, shell=True)
    prog.wait()

    if prog.returncode != 0:
        sys.stderr.write("Running bowtie failed\n")
        sys.exit(EXIT['BOWTIE_ERROR'])

    logging.debug("Finished running bowtie.")

    return bwt_results_file


# Used in main
def write_individual_results(outputDir, maxOffTargets, sortedOutput, guideSize, mode, totalClusters, limitPrintResults, offtargetsTable):
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
        if mode == ProgramMode.TALENS or mode == ProgramMode.NICKASE:
            clusterID = current.cluster
            clusters[clusterID-1].append(current)

        offTargets = current.as_off_target_string("", maxOffTargets)
        if not offTargets:
            offTargets = "There are no predicted off-targets."

        f.write(str(current.strandedGuideSeq)+"\n"+offTargets+"\n")

        if mode == ProgramMode.CRISPR and not config.isoforms and current.repStats is not None:
            stats_file = '%s/%s_repStats.json' % (outputDir, current.ID)
            with open(stats_file, 'w') as fp:
                json.dump(current.repStats, fp)
            fp.close()

        if mode == ProgramMode.CRISPR and not config.isoforms and current.repProfile is not None:
            profile_file = '%s/%s_repProfile.csv' % (outputDir, current.ID)
            current.repProfile.to_csv(profile_file, index=False)

        if mode == ProgramMode.CRISPR and not config.isoforms and offtargetsTable:
            off_table = '%s/offtargetsTable.csv' % outputDir
            label = "%s:%s,%s,%s" % (current.chrom, current.start, current.strand, current.strandedGuideSeq)
            off_for_table = map(lambda x: x.as_off_target_string(label, maxOffTargets), current.offTargets)
            with open(off_table, "a") as append_file:
                if len(list(off_for_table)) > 0:
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


# Used in main
def parse_fasta_target(fasta_file, candidate_fasta_file, target_size, eval_and_print):
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


# Used in main
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


# Used in main
def mode_select(var: any, index: str, mode: ProgramMode):
    """ Selects a default depending on mode for options that have not been set """
    if var is not None:
        return var

    if mode == ProgramMode.CRISPR:
        return CRISPR_DEFAULT[index]

    elif mode == ProgramMode.TALENS:
        return TALEN_DEFAULT[index]

    elif mode == ProgramMode.CPF1:
        return CPF1_DEFAULT[index]

    elif mode == ProgramMode.NICKASE:
        return NICKASE_DEFAULT[index]

    sys.stderr.write("Unknown model %s\n" % mode)
    sys.exit(EXIT['PYTHON_ERROR'])


# Used in main
def print_bed(mode, vis_cords, targets, output_file, description): # bed is 0-based
    bed_file = open(output_file, 'w')

    if mode == ProgramMode.CRISPR:
        thresholds = [0, 1000]
    elif mode == ProgramMode.CPF1:
        thresholds = [300, 1000]
    elif mode == ProgramMode.NICKASE:
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

            if mode == ProgramMode.CRISPR or mode == ProgramMode.CPF1:
                start = target[1] - 1
                stop = target[1] + target[3] - 1
            else:
                start = target[6] - 1
                stop = target[7] - 1

            bed_line = "{0}\t{1}\t{2}\tRanked:{3}\t{4}\t{5}\t{1}\t{2}\t{6}\n".format(chromosome, start, stop,
                                                                                     target[0], 0, target[4], color)
            bed_file.write(bed_line)

    bed_file.close()


# Used in main
def print_genbank(mode, name, seq, exons, targets, chrom, seq_start, seq_end, strand, output_file, description): # different than other dump_gb
    genbank_file = open(output_file, 'w')
    loci = chrom + ":" + str(seq_start) + "-" + str(seq_end)
    if len(name) > 10: # almost always... Genbank seems a bit outdated as format
        name = name[-10:]
    if len(loci) > 10: # almost always...
        loci = name[-10:]
    record = SeqRecord(Seq(seq), description=description,
                       name=name, id=loci)
    gene_strand = 1 if strand == "+" else -1
    # genbank is 0-based
    if len(targets) > 0:
        for target in targets:
            ts = 1 if target[4] == "+" else -1
            if config.isoforms:
                ts = gene_strand

            if mode == ProgramMode.CRISPR or mode == ProgramMode.CPF1:
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

        # No clue what this does, it might fix something to do with Bio.Alphabet deprecation stuff.
        # TODO: Have an actual bioinformatics person look at this.
        record.annotations["molecule_type"] = "DNA"

        SeqIO.write(record, genbank_file, "genbank")
    genbank_file.close()

# Used in FastaToViscord
def complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))

# Used in main
def fasta_to_viscoords(sequences, strand):
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


# Used in main
def cluster_pairs(pairs):
    """ Clusters paired sequences according to overlap, so user knows which TALE pairs are redundant """

    # Sets the starting pair of TALEs to be compared to
    first = pairs[0]
    cluster = 1
    first.cluster = cluster
    in_cluster = 0

    # Compares each TALE pair to previous pair in list to see whether redundant. Assigns cluster number accordingly
    for i in range(1,len(pairs)):
        cur = pairs[i]
        prev = pairs[i-1]

        # Specifically, compares location of spacer (by comparing location of tales) to see whether there is overlap,
        # and therefore TALE pairs are redundant
        if ((cur.spacerStart <= prev.spacerEnd) and (cur.spacerEnd >= prev.spacerStart) and
                    in_cluster < PRIMER_OFF_TARGET_MIN):

            cur.cluster = cluster
            in_cluster += 1
        else:
            # If not redundant, increase cluster number
            cluster += 1
            cur.cluster = cluster
            in_cluster = 0

    return cluster, pairs
