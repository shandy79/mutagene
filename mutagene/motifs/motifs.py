import re
import math
import scipy.stats as stats
import numpy as np
import pandas as pd

import pprint
import logging
logger = logging.getLogger(__name__)


nucleotides = "ACGT"  # this order of nucleotides is important for reversing
complementary_nucleotide = dict(zip(nucleotides, reversed(nucleotides)))
complementary_nucleotide['N'] = 'N'

bases_dict = {
    "A": "A", "G": "G", "T": "T", "C": "C",
    "W": "AT", "S": "CG", "M": "AC", "K": "GT", "R": "AG", "Y": "CT",
    "B": "TCG", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ATGC"}

extended_nucleotides = "ACTGWSMKRYBDHVN"

comp_dict = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "W": "AT", "S": "CG", "K": "AC", "M": "GT", "Y": "AG", "R": "CT",
    "V": "TCG", "H": "AGT", "D": "ACT", "B": "ACG", "N": "ATGC"}


motifs = [
    {
        'name': 'APOBEC1 and APOBEC3A/B',
        'logo': 'T[C>K]W',
        'motif': 'TCW',
        'position': 1,
        'ref': 'C',
        'alt': 'K',
        'references': ' Biochemistry  2011;76:131–46. Nat Immunol  2001;2:530–6. Biochim Biophys Acta  1992;1171:11–18'
    },
    {
        'name': 'APOBEC3G',
        'logo': 'C[C>K]R',
        'motif': 'CCR',
        'position': 1,
        'ref': 'C',
        'alt': 'K',
        'references': ' Biochemistry  2011;76:131–46'
    },
    {
        'name': 'Spontaneous G:C>A:T mutations',
        'logo': '[C>T]G',
        'motif': 'CG',
        'position': 0,
        'ref': 'C',
        'alt': 'T',
        'references': 'Hum Genet  1988;78:151–5'
    },
    {
        'name': 'UV Light',
        'logo': 'Y[C>T]',
        'motif': 'YC',
        'position': 1,
        'ref': 'C',
        'alt': 'T',
        'references': 'JNCL Natl Cancer Inst (2018)'
    },
    {
        'name': 'Pol Eta',
        'logo': 'W[A>T]',
        'motif': 'WA',
        'position': 1,
        'ref': 'A',
        'alt': 'T',
        'references': 'Nat Genet  2013;45:970–6'
    },
    {
        'name': 'AID',
        'logo': 'W[R>S]C',
        'motif': 'WRC',
        'position': 1,
        'ref': 'R',
        'alt': 'S',
        'references': 'Nature  2003;424:103–7'
    },
]


def identify_motifs(samples_mutations, custom_motif=None, strand=None):
    motif_matches = []

    if strand is None:
        strand = '*'
    if strand == '*':
        strand = '+-'

    if custom_motif:
        search_motifs = scanf_motif(custom_motif)
    else:
        search_motifs = motifs.copy()
    # search_motifs.extend(scanf_motif(custom_motif))

    for sample, mutations in samples_mutations.items():
        # print(sample, len(mutations))
        if mutations is not None and len(mutations) > 0:
            first_mut_seq_with_coords = mutations[0][-1]
            window_size = (len(first_mut_seq_with_coords) - 1) // 2

            for m in search_motifs:
                for s in strand:
                    # print("IDENTIFYING MOTIF: ", m['name'])
                    result = get_enrichment(mutations, m['motif'], m['position'], m['ref'], m['alt'], window_size, s)

                    debug_data = {'sample': sample, 'motif': m['logo'], 'strand': s}
                    debug_data.update(result)
                    debug_string = pprint.pformat(debug_data, indent=4)
                    logger.debug(debug_string)

                    if result['mutation_load'] == 0:
                        continue

                    motif_matches.append({
                        'sample': sample,
                        'name': m['name'],
                        'motif': m['logo'],
                        'strand': s,
                        'enrichment': result['enrichment'],
                        'pvalue': result['pvalue_fisher'],
                        'mutations_low_est': result['mutation_load'],
                        'mutations_high_est': result['bases_mutated_in_motif'],
                    })
    return motif_matches


def scanf_motif(custom_motif):
    """ recognize motif syntax like A[C>T]G and create a motif entry """
    m = re.search(
        r'([' + extended_nucleotides + ']*)\[([' + nucleotides + '])\>([' + extended_nucleotides + '])\]([' + extended_nucleotides + ']*)',
        custom_motif.upper())
    if m:
        g = m.groups('')
        # print("GROUPS", m.group(1), m.group(2), m.group(3), m.group(4))
        entry = {}
        entry['logo'] = m.group(0)
        entry['motif'] = g[0] + g[1] + g[3]
        entry['position'] = len(g[0])
        entry['ref'] = g[1]
        entry['alt'] = g[2]
        if entry['ref'] == entry['alt']:
            return []
        entry['name'] = 'Custom motif'
        entry['references'] = ''
        return [entry, ]
    return []


def get_stats(motif_mutation_count, mutation_count, motif_count, ref_count):
    contingency_table = np.array(
        [[mutation_count, motif_mutation_count],
         [ref_count, motif_count]])

    p_val_fisher = stats.fisher_exact(contingency_table, alternative="less")[1]
    try:
        p_val_chi2 = stats.chi2_contingency(contingency_table)[1]
    except ValueError:
        p_val_chi2 = 1.0
    # if p_value <= 0.05:
    #  qvalues = multipletests(pvals=p_value, method='fdr_bh')
    #  print(qvalues)
    #  if qvalues[3] <= 0.05:
    #     print("significant")
    # print("odds_ratio: ", "p-value")
    return p_val_fisher, p_val_chi2


def get_rev_comp_seq(sequence):
    # rev_comp_seq = "".join([complementary_nucleotide[i] for i in reversed(sequence)])
    rev_comp_seq = [(i[0], i[1], complementary_nucleotide[i[2]], "-") for i in reversed(sequence)]
    return rev_comp_seq


def mutated_base(mutation, ref, alt):
    """
    :param mutation: [(record.CHROM, record.POS, record.REF, record.ALT)]
    :param ref: the nucleotide base pre-mutation
    :param alt: the nucleotide base post-mutation
    :return: how many mutations match specified ref and alt
    """
    assert ref != alt, "mutation should have different ref and alt nucleotides"

    # makes sure single base substitution
    if mutation[3] and mutation[2] and len(mutation[2]) == 1 \
            and len(mutation[3]) == 1 and len(mutation[3]) == 1 \
            and len(mutation[2]) == 1 and mutation[2] != mutation[3]:
        # mutation matches the substitution
        if mutation[2] in ref and mutation[3] in alt:
            return True


def find_matching_motifs(seq, motif, motif_position):
    # print("Looking for motif {} in {}, {}".format(motif, sequence, len(sequence) - len(motif)))
    for i in range(len(seq) - len(motif) + 1):
        s = seq[i: i + len(motif)]
        for j, char in enumerate(motif):
            if s[j][2] not in bases_dict[char]:
                break
        else:
            yield seq[i + motif_position]


def find_matching_bases(seq, ref, motif, motif_position):
    # print("Looking for motif {} in {}, {}".format(motif, seq, len(seq) - len(motif)))
    ss = ""
    for i in range(motif_position, len(seq) - (len(motif) - motif_position)):
        ss += seq[i][2]
    # print(ref, ss)
    for i in range(motif_position, len(seq) - (len(motif) - motif_position)):
        s = seq[i][2]
        if s in bases_dict[ref]:
            yield seq[i]


def get_enrichment(mutations, motif, motif_position, ref, alt, range_size, strand):
    assert range_size >= 0
    assert len(ref) == 1
    assert len(alt) == 1
    assert 0 <= motif_position < len(motif)

    matching_bases = set()
    matching_motifs = set()
    matching_mutated_motifs = set()
    matching_mutated_bases = set()

    # extra loop for sample in sample list
    for chrom, pos, transcript_strand, x, y, seq in mutations:
        # extract the longest sequence we would ever need (motif + range_size)
        mutation = chrom, pos, x, y
        rev_seq = get_rev_comp_seq(seq)

        # if strand == '+':
        if transcript_strand == strand:
            # not mutated:
            for ref_match in find_matching_bases(seq, ref, motif, motif_position):
                matching_bases.add(ref_match[0:2])
            for motif_match in find_matching_motifs(seq, motif, motif_position):
                matching_motifs.add(motif_match[0:2])

            # mutated:
            if mutated_base(mutation, bases_dict[ref], bases_dict[alt]):
                # m = (mutation[0], mutation[1], mutation[2], "+")
                matching_mutated_bases.add(mutation[0:2])

                context_of_mutation = seq[range_size - motif_position: range_size - motif_position + len(motif)]
                for motif_match in find_matching_motifs(context_of_mutation, motif, motif_position):
                    matching_mutated_motifs.add(motif_match[0:2])

        # elif strand == '-':
        elif transcript_strand != strand:
            # rev compl: not mutated:
            for ref_match in find_matching_bases(rev_seq, ref, motif, len(motif) - motif_position - 1):
                matching_bases.add(ref_match[0:2])

            for motif_match in find_matching_motifs(rev_seq, motif, len(motif) - motif_position - 1):
                matching_motifs.add(motif_match[0:2])

            # rev compl: mutated:
            if mutated_base(mutation, comp_dict[ref], comp_dict[alt]):
                # m = (mutation[0], mutation[1], mutation[2], "-")
                matching_mutated_bases.add(mutation[0:2])

                # rev comp:
                context_of_mutation = rev_seq[range_size - motif_position: range_size - motif_position + len(motif)]
                for motif_match in find_matching_motifs(context_of_mutation, motif, len(motif) - motif_position - 1):
                    matching_mutated_motifs.add(motif_match[0:2])

        # if seq[0][0] == '19':
        #     print(seq)

        # if seq[0][0] == '19' and seq[0][1] == 51022172:
        #     for s in seq:
        #         print(s[2], end='')
        #     print()

    motif_mutation_count = len(matching_mutated_motifs)  # bases mutated in motif
    mutation_count = len(matching_mutated_bases - matching_mutated_motifs)  # bases mutated not in motif
    motif_count = len(matching_motifs)  # bases in motif
    ref_count = len(matching_bases - matching_motifs)  # bases not in motif

    stat_motif_count = len(matching_motifs - matching_mutated_motifs)  # bases not mutated in motif
    stat_ref_count = len(matching_bases - matching_motifs - matching_mutated_bases)  # bases not mutated not in motif

    # try:
    #     enrichment = (motif_mutation_count / mutation_count) / (motif_count / ref_count)
    # except ZeroDivisionError:
    #     enrichment = 0.0

    # haldane correction:
    motif_mutation_count += 0.5
    mutation_count += 0.5
    stat_motif_count += 0.5
    stat_ref_count += 0.5

    try:
        enrichment = (motif_mutation_count / mutation_count) / (stat_motif_count / stat_ref_count)
    except ZeroDivisionError:
        enrichment = 0.0

    p_val_fisher, p_val_chi2 = get_stats(motif_mutation_count, mutation_count, stat_motif_count, stat_ref_count)

    if enrichment > 1 and p_val_fisher < 0.05 and p_val_chi2 < 0.05:
        mut_load = (motif_mutation_count * (enrichment - 1)) / enrichment
    else:
        mut_load = 0.0

    table = pd.DataFrame(data={
        "'{}>{}' mutation".format(ref, alt): [motif_mutation_count, mutation_count],
        "no '{}>{}' mutation".format(ref, alt): [stat_motif_count, stat_ref_count]},
        index=("'{}' motif".format(motif), "no '{}' motif".format(motif)))
    logger.debug("\n" + table.to_string() + "\n")

    result = {
        'enrichment': enrichment,
        'mutation_load': math.ceil(mut_load),
        'pvalue_fisher': p_val_fisher,
        'pvalue_chi2': p_val_chi2,
        'bases_mutated_in_motif': motif_mutation_count,
        'bases_mutated_not_in_motif': mutation_count,
        'bases_in_motif': motif_count,
        'bases_not_in_motif': ref_count,
        'bases_not_mutated_in_motif': stat_motif_count,
        'bases_not_mutated_not_in_motif': stat_ref_count,
        'total_mutations': len(mutations)
    }

    return result