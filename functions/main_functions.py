import json
import logging
import re
import sys
import warnings
from subprocess import Popen, PIPE

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

import config
from classes.ProgramMode import ProgramMode
from constants import DOWNSTREAM_NUC, CPF1_DEFAULT, TALEN_DEFAULT, CRISPR_DEFAULT, EXIT, NICKASE_DEFAULT


# Here lies concatenate_feature_sets

# Used in main
def coord_to_fasta(regions, fasta_file, output_dir, target_size, eval_and_print_func, non_over, index_dir, genome,
                   strand, ext):
    """ Extracts the sequence corresponding to genomic coordinates from a FASTA file """

    ext = 0 if config.isoforms else ext  # for genomic context for some models
    sequences = {}
    fasta_file = open(fasta_file, 'w')
    fasta_seq = ""

    if config.isoforms and strand == "-":
        regions = regions[::-1]

    for region in regions:
        # Extracts chromosome number and region start and end
        chrom = region[0:region.rfind(':')]
        start = int(region[region.rfind(':') + 1:region.rfind('-')])
        finish = int(region[region.rfind('-') + 1:])
        start = max(start, 0)

        if ext == 0 and finish == start:
            continue

        # Run twoBitToFa program to get actual dna sequence corresponding to input genomic coordinates
        # Popen runs twoBitToFa program. PIPE pipes stdout.
        prog = Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2> %s/twoBitToFa.err" % (
            config.path("TWOBITTOFA"), chrom, start - ext, finish + ext, index_dir, genome, output_dir), stdout=PIPE,
                     shell=True)

        # Communicate converts stdout to a string
        output = prog.communicate()
        if prog.returncode != 0:
            sys.stderr.write("Running twoBitToFa failed\n")
            sys.exit(EXIT['TWOBITTOFA_ERROR'])

        output = output[0].decode()
        exons = output.split("\n")
        dna = ''.join(exons[1:]).upper()
        ext_dna = dna
        dna = dna[ext:(len(dna) - ext)]
        if len(dna) != (finish - start):  # something is wrong with what was fetched by twoBitToFa
            continue

        if config.isoforms and strand == "-":
            dna = str(Seq(dna).reverse_complement())

        # Write exon sequences to text file user can open in ApE. exon-intron junctions in lowercase.
        fasta_seq += dna[0].lower() + dna[1:-1] + dna[-1].lower()

        # Add 1 due to BED 0-indexing
        name = "C:%s:%d-%d" % (chrom, start, finish)

        # Loop over exon sequence, write every g-mer into file in which g-mer ends in PAM in fasta format
        positions = list(range(0, len(dna) - (target_size - 1)))
        while len(positions) != 0:
            num = positions.pop(0)
            downstream_5prim = ext_dna[num:(num + ext)]
            g_end = num + ext + target_size
            downstream_3prim = ext_dna[g_end:(g_end + ext)]
            if eval_and_print_func(name, target_size, dna[num:(num + target_size)],
                                   len(dna) - num - target_size if config.isoforms and strand == "-" else num,
                                   fasta_file, downstream_5prim, downstream_3prim):
                if non_over:  # positions overlapping those of this guide
                    for p in range(num, num + target_size):
                        if p in positions:
                            positions.remove(p)

                if name not in sequences:
                    sequences[name] = dna

    fasta_file.close()

    if config.isoforms and strand == "-":
        fasta_seq = str(Seq(fasta_seq).reverse_complement())

    return sequences, fasta_seq


# Used in main
def run_bowtie(pam_length, unique_method_cong, fasta_file, output_dir,
               max_off_targets, index_dir, genome, max_mismatches):
    logging.info("Running bowtie.")

    bwt_results_file = '%s/output.sam' % output_dir
    if unique_method_cong and not config.isoforms:
        # When ISOFORMS dna string is not reverse complemented and Cong can't be used
        # the -l alignment mode specifies a seed region to search for the number of mismatches specified with the
        # -n option. Outside of that seed, up to 2 mismatches are searched.
        # E.g. -l 15 -n 0 will search the first 15 bases with no mismatches, and the rest with up to 3 mismatches
        command = "%s -p %s -l %d -n %d -m %d --sam-nohead -k %d %s/%s -f %s -S %s " % (
            config.path("BOWTIE"), config.threads(), (pam_length + 11), max_mismatches, max_off_targets,
            max_off_targets, index_dir, genome, fasta_file, bwt_results_file)
    else:
        command = "%s -p %s -v %d --sam-nohead -k %d %s/%s -f %s -S %s " % (
            config.path("BOWTIE"), config.threads(), max_mismatches, max_off_targets, index_dir, genome, fasta_file,
            bwt_results_file)

    if config.isoforms:  # When ISFORMS we don't check reverse complement
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
def write_individual_results(output_dir: str, max_off_targets: int, sorted_output: list, program_mode: ProgramMode,
                             total_clusters: int, off_targets_table: bool) -> list:
    """ Writes each guide and its offtargets into a file """

    # Initiate list of lists for each cluster
    clusters = [[] for _ in range(total_clusters)]

    # Info to write to 'offtargets.json'
    off_targets_info = dict()

    for i, guide in enumerate(sorted_output):
        guide.id = i + 1

        if guide.id not in off_targets_info:
            off_targets = guide.as_off_target_string("", max_off_targets)

            if not off_targets:
                off_targets = "There are no predicted off-targets."

            off_targets_info[guide.id] = {'stranded guide seq': guide.stranded_guide_seq,
                                          'off-targets': off_targets}

        # Add the current TALE pair to the appropriate list in the list of lists, depending on its cluster number
        if program_mode == ProgramMode.TALENS or program_mode == ProgramMode.NICKASE:
            cluster_id = guide.cluster
            clusters[cluster_id - 1].append(guide)

        if program_mode == ProgramMode.CRISPR and not config.isoforms:
            if guide.repair_stats is not None:
                stats_file = f"{output_dir}/{guide.id}_repStats.json"

                with open(stats_file, 'w') as fp:
                    json.dump(guide.repair_stats, fp)

            if guide.repair_profile is not None:
                profile_file = f"{output_dir}/{guide.id}_repProfile.csv"
                guide.repair_profile.to_csv(profile_file, index=False)

            if off_targets_table:
                off_table = f"{output_dir}/offtargetsTable.csv"
                label = f"{guide.chrom}:{guide.start},{guide.strand},{guide.stranded_guide_seq}"
                off_for_table = map(lambda x: x.as_off_target_string(label, max_off_targets), guide.off_targets)

                with open(off_table, 'a') as append_file:
                    if len(list(off_for_table)) > 0:
                        append_file.write("\n".join(off_for_table))
                        append_file.write("\n")

            for cluster in clusters:
                if len(cluster) == 0:
                    continue

                best_in_cluster = cluster[0]

                for member in cluster[1:]:
                    # Write the other cluster members to file
                    off_targets_info[best_in_cluster.id]['cluster'] = "%s*%s*%s,%s:%s,%s,%s/%s,%s/%s,%s/%s,%s/%s;" % (
                        member.tale1.guide_seq, member.spacer_seq, member.tale2.guide_seq, member.chrom, member.start,
                        len(member.off_target_pairs), member.tale1.off_targets_mm[0], member.tale2.off_targets_mm[0],
                        member.tale1.off_targets_mm[1], member.tale2.off_targets_mm[1], member.tale1.off_targets_mm[2],
                        member.tale2.off_targets_mm[2], member.tale1.off_targets_mm[3], member.tale2.off_targets_mm[3]
                    )

                off_targets_info[best_in_cluster.id]['restriction sites'] = guide.restriction_sites

    for key, val in off_targets_info.items():
        with open(f"{output_dir}/{key}.offtargets.json", 'w') as f:
            json.dump(val, f)

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
    for num in range(0, len(sequence) - (target_size - 1)):

        if (num - DOWNSTREAM_NUC) > 0:
            start_5prim = num - DOWNSTREAM_NUC
        else:
            start_5prim = 0

        if (num + target_size + DOWNSTREAM_NUC) > len(sequence):
            end_3prim = len(sequence)
        else:
            end_3prim = num + target_size + DOWNSTREAM_NUC

        downstream_5prim = sequence[start_5prim:num]
        downstream_3prim = sequence[(num + target_size):end_3prim]

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
        db = MySQLdb.connect(user=m.group(1), passwd=m.group(2), host=m.group(3), db=m.group(4))
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
def print_bed(mode, vis_cords, targets, output_file, description):  # bed is 0-based
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
def print_genbank(mode, name, seq, exons, targets, chrom, seq_start, seq_end, strand, output_file,
                  description):  # different than other dump_gb
    genbank_file = open(output_file, 'w')
    loci = chrom + ":" + str(seq_start) + "-" + str(seq_end)
    if len(name) > 10:  # almost always... Genbank seems a bit outdated as format
        name = name[-10:]
    if len(loci) > 10:  # almost always...
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

            record.features.append(SeqFeature(FeatureLocation(start - seq_start, stop - seq_start,
                                                              strand=ts), type="Target_%s" % target[0]))

    if len(exons) > 0:
        for exon in exons:
            record.features.append(SeqFeature(FeatureLocation(exon[1] - seq_start, exon[2] - seq_start,
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
    exon_start = []
    exon_end = []
    exon_sequence = []

    for exon in sequences:
        # sys.stderr.write("%s\n" % exon)
        exon_list = exon.split(':')
        exon_coord = exon_list[2].split('-')
        exon_start.append(exon_coord[0])
        exon_end.append(exon_coord[1])
        seq = sequences[exon]
        if strand == "-":
            seq = complement(seq)

        exon_sequence.append(seq)

    return zip(exon_start, exon_end, exon_sequence)
