import os.path
import gzip
from functools import partial
import sys
import io
from Bio import SeqIO
from itertools import chain
import json


def store_gff_as_dict(gff_path, evalue_threshold=1e-06):
    """
    :param gff_path: path to gff file from barrnap
    :param evalue_threshold: max similarity e-value of predicted subunit to add to the dictionary
    :return: dictionary: {read_name: [(subunit_name, start_pos, end_pos, strand, e-value), ...]}
    """
    gff_dict = {}
    gff_set_dict = {}
    with open(gff_path) as gff:
        gff.readline()  # skip header
        for line in gff:
            read_name = line.split('\t')[0]
            start_pos = int(line.split('\t')[3]) - 1
            end_pos = int(line.split('\t')[4]) - 1
            strand = line.split('\t')[6]
            evalue = float(line.split('\t')[5])
            if evalue > float(evalue_threshold):
                continue
            subunit_name = line.split('\t')[8].split(';')[0].split('=')[1]
            gff_reduced_entry = (subunit_name, start_pos, end_pos, strand)
            if read_name in gff_dict:
                if gff_reduced_entry not in gff_set_dict[read_name]:
                    gff_dict[read_name].append(gff_reduced_entry)
                    gff_set_dict[read_name].add(gff_reduced_entry)
            else:
                gff_dict[read_name] = [gff_reduced_entry]
                gff_set_dict[read_name] = {gff_reduced_entry}
    return gff_dict


def pick_rdna_reads_to_assembly(gff_dict, out_prefix):
    """
    :param gff_dict: dictionary output of store_gff_as_dict() function
    :param out_prefix: files prefix to store stats and read ids
    :return: list with id of reads which contain 18S, 5.8S, 28S subunits (in this or reverse order)
    """
    grouped_rdna_reads = {}
    grouped_rdna_reads_w_diff_strands = {}

    for read_name in gff_dict:
        found_subunits = None
        # all rDNA subunits in one read are expected to have same strand orientation, do analysis separately otherwise
        strands = set([gff_entry[3] for gff_entry in gff_dict[read_name]])
        if len(strands) == 1:
            for gff_entry in gff_dict[read_name]:
                if found_subunits is None:
                    found_subunits = [gff_entry[0]]
                elif found_subunits[-1] == gff_entry[0] and gff_entry[0] != '5S_rRNA':
                    continue
                else:
                    found_subunits.append(gff_entry[0])

            if '+' in strands:
                read_rdna_composition = '-'.join(found_subunits).replace('_rRNA', '').replace('_', '.')
            else:
                read_rdna_composition = '-'.join(reversed(found_subunits)).replace('_rRNA', '').replace('_', '.')

            if read_rdna_composition in grouped_rdna_reads:
                grouped_rdna_reads[read_rdna_composition].append(read_name)
            else:
                grouped_rdna_reads[read_rdna_composition] = [read_name]
        else:
            found_subunits_strand = []
            for gff_entry in gff_dict[read_name]:
                if found_subunits is None:
                    found_subunits = [gff_entry[0]]
                    found_subunits_strand.append(gff_entry[3])
                elif found_subunits[-1] == gff_entry[0] and found_subunits_strand[-1] == gff_entry[3] and gff_entry[0] != '5S_rRNA':
                    continue
                else:
                    found_subunits.append(gff_entry[0])
                    found_subunits_strand.append(gff_entry[3])

            read_rdna_composition = '_'.join([subunit.replace('_rRNA', '').replace('_', '.') + '(' + strand + ')'
                                              for subunit, strand in zip(found_subunits, found_subunits_strand)])

            if read_rdna_composition in grouped_rdna_reads_w_diff_strands:
                grouped_rdna_reads_w_diff_strands[read_rdna_composition].append(read_name)
            else:
                grouped_rdna_reads_w_diff_strands[read_rdna_composition] = [read_name]

    # Display and save rDNA read quantities
    # NOR rDNA (18S, 5.8S, 28S)
    nor_rdna_compositions = [composition for composition in grouped_rdna_reads.keys()
                             if '5S' not in composition and '5.8S' in composition
                             and len(composition.split('-')) >= 2 and '28S-5.8S' not in composition
                             and '18S-28S' not in composition]
    nor_rdna_compositions_for_assembly = [composition for composition in nor_rdna_compositions
                                          if '5S' not in composition and '18S-5.8S-28S' in composition]
    nor_rdna_compositions = sorted(nor_rdna_compositions, key=lambda x: len(grouped_rdna_reads[x]), reverse=True)
    stats_strings = []
    stats_strings.append('\n====== NOR rDNA reads ======\n'
                         'Nucleolus organiser region (NOR) rDNA reads have no 5S subunits,\n'
                         'cover at least one 5.8S and either 18S or 28S\n'
                         'The reads choosen for assembly (*) cover at least one\n'
                         'full rDNA operon (18S-5.8S-28S) in correct order and no 5S')
    for rdna_composition in nor_rdna_compositions:
        choosen_one = '(*)' if rdna_composition in nor_rdna_compositions_for_assembly else '( )'
        stats_strings.append('\t' + choosen_one + str(len(grouped_rdna_reads[rdna_composition])) + 'x(' + rdna_composition + ')')
    stats_strings.append('Total NOR long reads: ' + str(sum([len(grouped_rdna_reads[composition])
                                              for composition in nor_rdna_compositions])))
    stats_strings.append('Total NOR long reads for assembly: ' + str(sum([len(grouped_rdna_reads[composition])
                                                           for composition in nor_rdna_compositions_for_assembly])))

    # Other rDNA (5S arrays)
    other_rdna_compositions = [composition for composition in grouped_rdna_reads.keys()
                               if composition.count('5S') >= 2 and len(set(composition.split('-'))) == 1]
    other_rdna_compositions = sorted(other_rdna_compositions, key=lambda x: len(grouped_rdna_reads[x]), reverse=True)
    stats_strings.append('\n====== 5S rDNA reads ======\n'
                         'The reads have at least two 5S rDNA and no other rDNA')
    for rdna_composition in other_rdna_compositions:
        stats_strings.append('\t' + str(len(grouped_rdna_reads[rdna_composition])) + 'x(' + rdna_composition + ')')
    stats_strings.append('Total 5S long reads: ' + str(sum([len(grouped_rdna_reads[composition])
                                              for composition in other_rdna_compositions])))

    # Ambiguous reads (e.g., with 5S and NOR rDNA in one read)
    amb_rdna_compositions = [composition for composition in grouped_rdna_reads.keys()
                             if composition not in nor_rdna_compositions
                             and composition not in other_rdna_compositions]
    amb_rdna_compositions = sorted(amb_rdna_compositions, key=lambda x: len(grouped_rdna_reads[x]), reverse=True)
    stats_strings.append('\n====== Dubious rDNA reads ======\n'
                         'The reads have non-canonical composition (subunits order, or combination),\n'
                         'they likely originate from mitochondrial rRNA genes / rRNA pseudogenes,\n'
                         'cover rare rRNA subunit variation (to the point it is not predicted by HMM-model),\n'
                         'or the reads are too short to contain enough rDNA subunits')
    for rdna_composition in amb_rdna_compositions:
        stats_strings.append('\t' + str(len(grouped_rdna_reads[rdna_composition])) + 'x(' + rdna_composition + ')')
    stats_strings.append('Total dubious reads: ' + str(sum([len(grouped_rdna_reads[composition])
                                              for composition in amb_rdna_compositions])))

    # Reads with differently oriented subunits
    rdna_reads_w_diff_strands_compositions = [composition for composition in grouped_rdna_reads_w_diff_strands.keys()]
    rdna_reads_w_diff_strands_compositions = sorted(rdna_reads_w_diff_strands_compositions,
                                                    key=lambda x: len(grouped_rdna_reads_w_diff_strands[x]), reverse=True)
    stats_strings.append('\n====== Reads with differently orientated rDNA subunits ======\n'
                         'The predicted rDNA subunits in these reads have different strand orientations.\n'
                         'While it is expected that all subunits in one rDNA operon (5S array) have same strand orientation,\n'
                         'cases with rRNA gene variants with reversed subunits have been reported.\n'
                         'Another explanation could be erroneous HMM model prediction, or bad PacBio adapter detection.\n'
                         'It worth to consider the quantity of the reads and their composition.')
    for rdna_composition in rdna_reads_w_diff_strands_compositions:
        stats_strings.append('\t' + str(len(grouped_rdna_reads_w_diff_strands[rdna_composition]))
                             + 'x(' + rdna_composition + ')')
    stats_strings.append('Total reads having rDNA subunits with different strand orientation: '
                         + str(sum([len(grouped_rdna_reads_w_diff_strands[composition])
                                    for composition in rdna_reads_w_diff_strands_compositions]))
                         + '\n====================================\n')

    with open(out_prefix + '.stats', 'w') as out_stats:
        out_stats.write('\n'.join(stats_strings))

    with open(out_prefix + '_for_assembly.list', 'w') as out_id_list:
        out_id_list.write('\n'.join([read_id for read_ids in [grouped_rdna_reads[composition]
                                                              for composition in nor_rdna_compositions_for_assembly]
                                     for read_id in read_ids]))

    all_rdna_reads_dict = dict(chain.from_iterable(d.iteritems() for d in (grouped_rdna_reads, grouped_rdna_reads_w_diff_strands)))
    with open(out_prefix + '.json', 'w') as out_json:
        json.dump(all_rdna_reads_dict, out_json, sort_keys=True, indent=4)


def extract_rdna_reads(reads_path, reads_id_list_path, output_path):
    """
    :param reads_path: path to long reads
    :param reads_id_list_path: path to list containing selected read ids
    :param output_path: path to output reads
    """
    reads_extention = os.path.splitext(reads_path)[1].lower()
    _open = open
    if reads_extention == '.gz':
        _open = partial(gzip.open, mode='rt')
        reads_extention = os.path.splitext(os.path.splitext(reads_path)[0])[1].lower()
    if reads_extention in ['.fasta', '.fa', '.fna']:
        reads_encoding = 'fasta'
    elif reads_extention in ['.fastq', '.fq']:
        reads_encoding = 'fastq'
    else:
        sys.exit('WARNING: long reads file '
                 + reads_path +
                 ' has unknown format. Allowed formats: fasta, fa, fna, fastq, fq '
                 '(not case-sensitive, can be compressed with gzip). Skipping file...')

    with _open(reads_path) as reads_opener, open(reads_id_list_path) as reads_id_list, open(output_path, 'w') as output_fasta:
        reads_id_list = {read_id.strip() for read_id in reads_id_list.readlines()}
        reads = io.BufferedReader(reads_opener)
        for read in SeqIO.parse(reads, reads_encoding):
            if read.id in reads_id_list:
                SeqIO.write(read, output_fasta, 'fasta')

