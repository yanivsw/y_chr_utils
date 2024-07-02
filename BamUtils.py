from dataclasses import dataclass, fields
from pydantic import validate_arguments
import numpy as np
import pandas as pd
import subprocess
import re

@validate_arguments
@dataclass
class BamSequence:
    QNAME: str  # Query template NAME
    FLAG: int   # bitwise FLAG
    RNAME: str  # Reference sequence NAME
    POS: int    # 1-based leftmost mapping POSition
    MAPQ: int   # MAPping Quality
    CIGAR: str  # CIGAR string
    RNEXT: str  # Reference name of the mate/next read
    PNEXT: str  # Position of the mate/next read
    TLEN: str   # observed Template LENgth
    SEQ: str    # segment SEQuence
    QUAL: str   # ASCII of Phred-scaled base QUALity+33
    RG: str = ''     # Read group (use to check sequencing library)
    MD: str = ''     # String encoding mismatched and deleted reference bases

@dataclass
class DecodedBamFlag:
    read_paired: int = 0
    read_mapped_in_proper_pair: int = 0
    read_unmapped: int = 0
    mate_unmapped: int = 0
    read_reverse_strand: int = 0
    mate_reverse_strand: int = 0
    first_in_pair: int = 0
    second_in_pair: int = 0
    not_primary_alignment: int = 0
    read_fails_quality_checks: int = 0
    read_is_duplicate: int = 0
    supplementary_alignment: int = 0

    SamFlags = [["read paired", 0x1],
                ["read mapped in proper pair", 0x2],
                ["read unmapped", 0x4],
                ["mate unmapped", 0x8],
                ["read reverse strand", 0x10],
                ["mate reverse strand", 0x20],
                ["first in pair", 0x40],
                ["second in pair", 0x80],
                ["not primary alignment", 0x100],
                ["read fails platform/vendor quality checks", 0x200],
                ["read is PCR or optical duplicate", 0x400],
                ["supplementary alignment", 0x800]]

    def decode_flag(self, input_flag):
        for decoded_flag_field, flag_info in zip(reversed(fields(self)), reversed(self.SamFlags)):
            if input_flag < flag_info[1]:
                continue
            else:
                bitshift = len(bin(flag_info[1])[2:]) - 1
                decoded_flag = (input_flag & flag_info[1]) >> bitshift
                setattr(self, decoded_flag_field.name, decoded_flag)
        return self

def write_bam(bam_file_out, seq_str_array, header_str):
    out_seq_str = '\n'.join(seq_str_array) + '\n'
    bam_outstr = header_str + out_seq_str
    samtools_str = "samtools view -bh -o {} -".format(bam_file_out)
    subprocess.run(samtools_str, input=bam_outstr, shell=True, capture_output=True, text=True)

def read_bam(bam_file):
    samtools_str = "samtools view {}".format(bam_file)
    samtools_output = subprocess.run(samtools_str, shell=True, capture_output=True, text=True)

    samtools_str_header = "samtools view -H {}".format(bam_file)
    samtools_output_header = subprocess.run(samtools_str_header, shell=True, capture_output=True, text=True)

    header_str = samtools_output_header.stdout
    sequence_strings = samtools_output.stdout.split("\n")[:-1]

    return {
        'header' : header_str,
        'sequence_strings' : sequence_strings
    }

def parse_md(md, reconstructed_alignment):
    reconstructed_reference = ''
    md = list(filter(None, re.split('(\d+)', md)))
    ref_pos = 0
    for x in md:
        if (x != '0'):
            if x.isdigit():
                reconstructed_reference += reconstructed_alignment[ref_pos:(ref_pos + int(x))]
                ref_pos += int(x)
            elif (x[0] != '^'):
                reconstructed_reference += x
                ref_pos += len(x)
            elif x[0] == '^':
                reconstructed_reference += x[1:]
                ref_pos += len(x[1:])

    return reconstructed_reference

def parse_cigar(seq_str, cigar_str, qual_str=None):
    cigar = tuple(zip(*[iter(re.findall(r'[^\W\d_]+|\d+', cigar_str))] * 2))

    # reconstructed_base_qual = ''
    # full_reconstructed_alignment = ''
    # full_reconstructed_base_qual = ''
    # insertions = []

    reconstructed_alignment = ''
    alignment_pos = 0
    for x in cigar:
        temp = int(x[0])
        if ('M' in x) | ('=' in x) | ('X' in x):
            reconstructed_alignment += seq_str[alignment_pos:(alignment_pos + temp)]
            # reconstructed_base_qual += qual_str[alignment_pos:(alignment_pos + temp)]

            # full_reconstructed_alignment += seq_str[alignment_pos:(alignment_pos + temp)]
            # full_reconstructed_base_qual += qual_str[alignment_pos:(alignment_pos + temp)]

            alignment_pos += temp
        elif ('I' in x) | ('S' in x):
            # full_reconstructed_alignment += seq_str[alignment_pos:(alignment_pos + temp)]
            # full_reconstructed_base_qual += qual_str[alignment_pos:(alignment_pos + temp)]

            # insertions.append([alignment_pos, temp, seq_str[alignment_pos:(alignment_pos + temp)]])

            alignment_pos += temp
        elif ('D' in x) | ('N' in x):
            reconstructed_alignment += '-' * temp
            # reconstructed_base_qual += '!' * temp

            # full_reconstructed_alignment += '-' * temp
            # full_reconstructed_base_qual += '!' * temp
        elif ('H' in x) | ('P' in x):
            continue

    return reconstructed_alignment

def is_seq_deaminated(read_reverse_strand, reconstructed_reference, reconstructed_alignment):
    ref_ends = reconstructed_reference[0:3] + reconstructed_reference[-3:]
    seq_ends = reconstructed_alignment[0:3] + reconstructed_alignment[-3:]

    if ref_ends != seq_ends:
        for ref_pos, seq_pos in zip(ref_ends, seq_ends):
            if (read_reverse_strand == 0) & (ref_pos == 'C') & (seq_pos == 'T'):
                return True
            if (read_reverse_strand == 1) & (ref_pos == 'G') & (seq_pos == 'A'):
                return True

    return False

def get_read_lengths(bam_file):
    samtools_str = "samtools view {}".format(bam_file)
    samtools_output = subprocess.run(samtools_str, shell=True, capture_output=True, text=True)

    output = samtools_output.stdout.split("\n")[:-1]
    sequences = [i.split('\t') for i in output]

    seq_length = []
    for seq in sequences:
        sequence = {'SEQ' : seq[9],
                    'RG'  : seq[[idx for idx, s in enumerate(seq) if 'RG:' in s][0]].split(':')[-1]}

        seq_length.append([sequence['SEQ'], len(sequence['SEQ']), sequence['RG']])

    seq_length_df = pd.DataFrame(seq_length, columns=['sequence', 'sequence_length', 'RG'])

    seq_length_stats = seq_length_df[['sequence_length', 'RG']].groupby(['RG']).agg(
                                        count = ('sequence_length', 'count'),
                                        median = ('sequence_length', np.median),
                                        mean = ('sequence_length', np.mean),
                                        max = ('sequence_length', np.max),
                                        min = ('sequence_length', np.min))

    seq_length_stats.loc['Combined'] = [len(seq_length_df['sequence_length']),
                                        seq_length_df['sequence_length'].median(),
                                        seq_length_df['sequence_length'].mean(),
                                        seq_length_df['sequence_length'].max(),
                                        seq_length_df['sequence_length'].min()]

    return seq_length_df, seq_length_stats

def get_read_lengths_deaminated(bam_file):
    samtools_str = "samtools view {}".format(bam_file)
    samtools_output = subprocess.run(samtools_str, shell=True, capture_output=True, text=True)

    output = samtools_output.stdout.split("\n")[:-1]
    sequences = [i.split('\t') for i in output]

    clipped_seqs = []
    out_seqs = []
    for seq in sequences:
        sequence = {'CIGAR' : seq[5],
                    'SEQ'   : seq[9],
                    'FLAG'  : int(seq[1]),
                    'MD'    : seq[[idx for idx, s in enumerate(seq) if 'MD:' in s][0]].split(':')[2],
                    'RG'    : seq[[idx for idx, s in enumerate(seq) if 'RG:' in s][0]].split(':')[2]}

        if 'S' in sequence['CIGAR']:
            clipped_seqs.append(sequence)
            continue

        reconstructed_alignment = parse_cigar(sequence['SEQ'], sequence['CIGAR'])
        reconstructed_reference = parse_md(sequence['MD'], reconstructed_alignment)

        decoded_flag = DecodedBamFlag().decode_flag(sequence['FLAG'])

        if is_seq_deaminated(decoded_flag.read_reverse_strand, reconstructed_reference, reconstructed_alignment):
            out_seqs.append([sequence['SEQ'], len(sequence['SEQ']), sequence['RG']])

    seq_length_df = pd.DataFrame(out_seqs, columns=['sequence', 'sequence_length', 'RG'])

    seq_length_stats = seq_length_df[['sequence_length', 'RG']].groupby(['RG']).agg(
                                        count = ('sequence_length', 'count'),
                                        median = ('sequence_length', np.median),
                                        mean = ('sequence_length', np.mean),
                                        max = ('sequence_length', np.max),
                                        min = ('sequence_length', np.min))

    seq_length_stats.loc['Combined'] = [len(seq_length_df['sequence_length']),
                                        seq_length_df['sequence_length'].median(),
                                        seq_length_df['sequence_length'].mean(),
                                        seq_length_df['sequence_length'].max(),
                                        seq_length_df['sequence_length'].min()]

    return seq_length_df, seq_length_stats

def print_informative_positions_results(input_df, state1, state2, bam_file, deaminated_reads, output_per_rg=False):
    from statsmodels.stats import proportion
    pd.options.mode.chained_assignment = None

    # Get duplicated sequences and generate a DataFrame of unique sequences
    duplicate_seqs = input_df[input_df.duplicated(subset=['QNAME'], keep=False)].groupby(['QNAME']).agg(tuple).reset_index()
    removed_duplicate_seqs = input_df.drop_duplicates(subset=['QNAME'], keep=False)

    state1_count = len(removed_duplicate_seqs[removed_duplicate_seqs['state'] == state1])
    state2_count = len(removed_duplicate_seqs[removed_duplicate_seqs['state'] == state2])
    conflicting_count = 0

    # Group the unique sequences by read group and count occurrences of each state
    count_per_rg = pd.DataFrame(removed_duplicate_seqs.groupby('RG')['state'].value_counts()).rename(columns={'state':'count'})

    # Process duplicate sequences
    if len(duplicate_seqs) > 0:
        # Identify conflicting sequences (those containing both state1 and state2) and get unique sequences
        conflicting = duplicate_seqs['state'].apply(lambda x: (state2 in x) & (state1 in x))
        non_conflicting_seqs = duplicate_seqs[~conflicting]
        conflicting_seqs = duplicate_seqs[conflicting]

        state1_count += len(non_conflicting_seqs[non_conflicting_seqs['state'].apply(lambda x: state1 in x)])
        state2_count += len(non_conflicting_seqs[non_conflicting_seqs['state'].apply(lambda x: state2 in x)])

        # Extract the read group and state for non-conflicting sequences
        non_conflicting_seqs['RG'] = [x[0] for x in non_conflicting_seqs['RG'].to_list()]
        non_conflicting_seqs['state'] = [x[0] for x in non_conflicting_seqs['state'].to_list()]

        # Count occurrences of each state in non-conflicting sequences and add to
        non_conflicting_seqs_count_per_rg = pd.DataFrame(non_conflicting_seqs.groupby('RG')['state'].value_counts()).rename(columns={'state':'count'})
        count_per_rg = count_per_rg.add(non_conflicting_seqs_count_per_rg, fill_value=0)

        # Process conflicting sequences
        if len(conflicting_seqs) > 0:
            conflicting_count += len(conflicting_seqs)

            count_per_rg.index = count_per_rg.index.set_levels([state2, state1, 'conflicting'], level=1)

            conflicting_seqs['RG'] = [sublist[0] for sublist in conflicting_seqs['RG'].to_list()]
            conflicting_seqs['state'] = 'conflicting'

            conflicting_seqs_count_per_rg = pd.DataFrame(conflicting_seqs.groupby('RG')['state'].value_counts()).rename(columns={'state':'count'})
            count_per_rg = count_per_rg.add(conflicting_seqs_count_per_rg, fill_value=0).astype(int)

    total_count = state1_count + state2_count + conflicting_count

    if total_count > 0:
        state1_proportion = (state1_count / total_count) * 100
        state2_proportion = (state2_count / total_count) * 100
        conflicting_proportion = (conflicting_count / total_count) * 100

        # df = pd.DataFrame(processed_sequences, columns=['pos', 'QNAME', 'RG', 'state', 'state1_allele', 'state2_allele', 'base', 'read_reverse_strand', 'base_index', 'reconstructed_alignment'])
        # df.to_csv('res.csv', sep='\t')

        state1_CI = proportion.proportion_confint(count = state1_count, nobs = total_count, alpha = 0.05)
        state2_CI = proportion.proportion_confint(count = state2_count, nobs = total_count, alpha = 0.05)
        conflicting_CI = proportion.proportion_confint(count = conflicting_count, nobs = total_count, alpha = 0.05)

        out_str = """BAM file: {}
Deaminated reads: {}
Num {} = {} \t % {} = {:.4f} \t 95% CI = ({:.2f}, {:.2f})
Num {} = {} \t % {} = {:.4f} \t 95% CI = ({:.2f}, {:.2f})
Num conflicting = {} \t % conf = {:.4f} \t 95% CI = ({:.2f}, {:.2f})""".format(bam_file, deaminated_reads,
                                                                               state1, state1_count, state1, state1_proportion, state1_CI[0]*100, state1_CI[1]*100,
                                                                               state2, state2_count, state2, state2_proportion, state2_CI[0]*100, state2_CI[1]*100,
                                                                               conflicting_count, conflicting_proportion, conflicting_CI[0]*100, conflicting_CI[1]*100)
        print(out_str)
    else:
        out_str = """BAM file: {}
Deaminated reads: {}
****No overlapping reads****
"""
        print(out_str)

    if output_per_rg:
        for rg in count_per_rg.index.get_level_values('RG').unique():
            rg_data = count_per_rg.loc[rg].reset_index()
            total_count = rg_data['count'].sum()

            if total_count > 0:
                state1_count = 0
                state2_count = 0
                conflicting_count = 0

                if state1 in rg_data['state'].values:
                    state1_count = rg_data[rg_data['state'] == state1]['count'].values[0]
                if state2 in rg_data['state'].values:
                    state2_count = rg_data[rg_data['state'] == state2]['count'].values[0]
                if 'conflicting' in rg_data['state'].values:
                    conflicting_count = rg_data[rg_data['state'] == 'conflicting']['count'].values[0]

                state1_proportion = (state1_count / total_count) * 100
                state2_proportion = (state2_count / total_count) * 100
                conflicting_proportion = (conflicting_count / total_count) * 100

                state1_CI = proportion.proportion_confint(count = state1_count, nobs = total_count, alpha = 0.05)
                state2_CI = proportion.proportion_confint(count = state2_count, nobs = total_count, alpha = 0.05)
                conflicting_CI = proportion.proportion_confint(count = conflicting_count, nobs = total_count, alpha = 0.05)

                out_str = '''
    Library - {}
    Num {} = {} \t % {} = {:.4f} \t 95% CI = ({:.2f}, {:.2f})
    Num {} = {} \t % {} = {:.4f} \t 95% CI = ({:.2f}, {:.2f})
    Num conflicting = {} \t % conf = {:.4f} \t 95% CI = ({:.2f}, {:.2f})'''.format(rg,
                                                                                   state1, int(state1_count), state1, state1_proportion, state1_CI[0]*100, state1_CI[1]*100,
                                                                                   state2, int(state2_count), state2, state2_proportion, state2_CI[0]*100, state2_CI[1]*100,
                                                                                   conflicting_count, int(conflicting_proportion), conflicting_CI[0]*100, conflicting_CI[1]*100)
                print(out_str)
            else:
                out_str = '''
    Library - {}
    ****No overlapping reads****'''
                print(out_str)

def count_informative_positions(positions, bam_file, deaminated_reads, mask, stranded, output_per_rg=False):
    states = list(positions.columns.values)[1:]
    state1 = states[0]
    state2 = states[1]

    processed_sequences = []
    clipped_seqs = []
    for index, row in positions.iterrows():
        pos = row['position']
        state1_allele = row[state1]
        state2_allele = row[state2]

        samtools_str = "samtools view {} Y:{}-{}".format(bam_file, pos, pos)
        samtools_output = subprocess.run(samtools_str, shell=True, capture_output=True, text=True)

        output = samtools_output.stdout.split("\n")[:-1]
        sequences = [i.split('\t') for i in output]

        for seq in sequences:
            sequence = BamSequence(*(seq[0:11] +
                                    [seq[[idx for idx, s in enumerate(seq) if 'RG:' in s][0]].split(':')[-1]] +
                                    [seq[[idx for idx, s in enumerate(seq) if 'MD:' in s][0]].split(':')[-1]]))

            if 'S' in sequence.CIGAR:
                clipped_seqs.append(sequence)
                continue

            reconstructed_alignment = parse_cigar(sequence.SEQ, sequence.CIGAR)

            base_index = abs(pos - sequence.POS)
            base = reconstructed_alignment[base_index]

            # sequence_base_qual = [(ord(i) - 33) for i in list(reconstructed_base_qual)]
            # full_sequence_base_qual = [(ord(i) - 33) for i in list(full_reconstructed_base_qual)]
            # base_qual = sequence_base_qual[base_index]

            decoded_flag = DecodedBamFlag().decode_flag(sequence.FLAG)

            # Only count a position if the read is deaminated
            if deaminated_reads:
                reconstructed_reference = parse_md(sequence.MD, reconstructed_alignment)

                if not is_seq_deaminated(decoded_flag.read_reverse_strand, reconstructed_reference, reconstructed_alignment):
                    continue

            # For stranded option as in count_informative_positions.pl
            # Don't count bases:
                # On the -ve strand (read_reverse_strand = 1) if one of the informative states is G
                # On the +ve strand (read_reverse_strand = 0) if one of the informative states is C
            if stranded:
                if (decoded_flag.read_reverse_strand == 0) & ( (state1_allele == 'C') | (state2_allele == 'C') ):
                    # Don't count position if either position = C
                    continue
                if (decoded_flag.read_reverse_strand == 1) & ( (state1_allele == 'G') | (state2_allele == 'G') ):
                    # Don't count position if either position = G
                    continue

            if mask > 0:
                base_masked = (base_index < mask) | (base_index >= (len(reconstructed_alignment) - mask))
                if base_masked:
                    if (decoded_flag.read_reverse_strand == 0) & (base == 'T'):
                        continue
                    if (decoded_flag.read_reverse_strand == 1) & (base == 'A'):
                        continue

            state = ''
            if base == state1_allele:
                state = state1
            elif base == state2_allele:
                state = state2

            # extra_info = [state1_allele, state2_allele, base, decoded_flag.read_reverse_strand, base_index, sequence.MD, reconstructed_reference, sequence.CIGAR, reconstructed_alignment, full_reconstructed_alignment]
            extra_info = []
            if state:
                processed_sequences.append([pos, sequence.QNAME, sequence.RG, state] + extra_info)

    # df_out = pd.DataFrame(processed_sequences, columns=['pos', 'QNAME', 'RG', 'state', 'state1_allele', 'state2_allele', 'base', 'read_reverse_strand', 'base_index', 'MD', 'reconstructed_reference', 'CIGAR', 'reconstructed_alignment', 'full_reconstructed_alignment'])
    df_out = pd.DataFrame(processed_sequences, columns=['pos', 'QNAME', 'RG', 'state'])

    print_informative_positions_results(df_out, state1, state2, bam_file, deaminated_reads, output_per_rg)

def get_substitution_patterns(bam_file, end_length):
    samtools_str = "samtools view {}".format(bam_file)
    samtools_output = subprocess.run(samtools_str, shell=True, capture_output=True, text=True)

    output = samtools_output.stdout.split("\n")[:-1]
    sequences = [i.split('\t') for i in output]

    translation_table = str.maketrans('ATCG','TAGC')

    clipped_seqs = []
    counts = []
    for idx, seq in enumerate(sequences):
        sequence = {'CIGAR' : seq[5],
                    'SEQ'   : seq[9],
                    'FLAG'  : int(seq[1]),
                    'MD'    : seq[[idx for idx, s in enumerate(seq) if 'MD:' in s][0]].split(':')[-1],
                    'RG'    : seq[[idx for idx, s in enumerate(seq) if 'RG:' in s][0]].split(':')[-1]}

        if 'S' in sequence['CIGAR']:
            clipped_seqs.append(sequence)
            continue

        reconstructed_alignment = parse_cigar(sequence['SEQ'], sequence['CIGAR'])
        reconstructed_reference = parse_md(sequence['MD'], reconstructed_alignment)

        # insertions = []
        # if insertions:
        #     for insertion in insertions:
        #         reconstructed_reference = ''.join((reconstructed_reference[:insertion[0]], '-' * insertion[1], reconstructed_reference[insertion[0]:]))
        #         reconstructed_alignment = ''.join((reconstructed_alignment[:insertion[0]], insertion[2], reconstructed_alignment[insertion[0]:]))

        #         if len(insertions) > 1:
        #             insertions = [[x[0] + insertion[1], x[1], x[2]] for x in insertions]

        if len(reconstructed_reference) != len(reconstructed_alignment):
            print('Diff length ' + reconstructed_reference + ' ' + reconstructed_alignment)

        decoded_flag = DecodedBamFlag().decode_flag(sequence['FLAG'])
        if (decoded_flag.read_reverse_strand == 1):
            reconstructed_reference = reconstructed_reference[::-1].translate(translation_table)
            reconstructed_alignment = reconstructed_alignment[::-1].translate(translation_table)

        ref_ends = reconstructed_reference[0:end_length] + reconstructed_reference[-end_length:]
        seq_ends = reconstructed_alignment[0:end_length] + reconstructed_alignment[-end_length:]

        counts.append([[ref_base + seq_base for ref_base, seq_base in zip(ref_ends, seq_ends)], sequence['RG']])

    counts_df = pd.DataFrame(counts, columns=['substitutions', 'RG'])

    return counts_df

def get_substitution_proportions(input):
    a_cols = ['A-', 'AA', 'AC', 'AG', 'AT']
    c_cols = ['C-', 'CA', 'CC', 'CG', 'CT']
    g_cols = ['G-', 'GA', 'GC', 'GG', 'GT']
    t_cols = ['T-', 'TA', 'TC', 'TG', 'TT']

    plot_cols = ['AC', 'AG', 'AT', 'CA', 'CG', 'GC', 'GT', 'TA', 'TC', 'TG', 'CT', 'GA']

    counts = pd.DataFrame(input).apply(pd.value_counts).T
    counts.columns = [col.replace('N', '-') for col in counts.columns.values]

    out = pd.DataFrame()
    out[a_cols] = counts[a_cols].div(counts[a_cols].sum(axis=1), axis=0)
    out[c_cols] = counts[c_cols].div(counts[c_cols].sum(axis=1), axis=0)
    out[g_cols] = counts[g_cols].div(counts[g_cols].sum(axis=1), axis=0)
    out[t_cols] = counts[t_cols].div(counts[t_cols].sum(axis=1), axis=0)

    return out[plot_cols]

def plot_deamination(ax1, ax2, count, end_length, y_maximum):
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    df_5 = count.T.iloc[:, 0:end_length].T
    df_3 = count.T.iloc[:, -end_length:].T

    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for col in df_5.columns:
        test_5 = df_5[col]
        test_3 = df_3[col].reset_index(drop=True)

        if col == 'CT':
            ax1.plot(test_5.index.values, test_5.values, linewidth=2, color=colours[0], label='C-T')
            ax2.plot(test_3.index.values, test_3.values, linewidth=2, color=colours[0], label='C-T')
        if col == 'GA':
            ax1.plot(test_5.index.values, test_5.values, linewidth=2, color=colours[1], label='G-A')
            ax2.plot(test_3.index.values, test_3.values, linewidth=2, color=colours[1], label='G-A')
        else:
            ax1.plot(test_5.index.values, test_5.values, linewidth=0.5, color='grey', label='other')
            ax2.plot(test_3.index.values, test_3.values, linewidth=0.5, color='grey', label='other')

    ax1.set_ylim([-0.01, y_maximum])
    ax1.grid()
    ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(5))
    ax1.set_xlabel('Position', fontsize=14)
    ax1.set_ylabel('Fraction', fontsize=14)
    ax1.set_title('5\' end', fontsize=14)

    ax2.set_ylim([-0.01, y_maximum])
    ax2.grid()
    ax2.xaxis.set_major_locator(mpl.ticker.MultipleLocator(5))
    ax2.set_xticklabels(['', '-20', '-15', '-10', '-5', '0'])
    ax2.set_xlabel('Position', fontsize=14)
    ax2.set_ylabel('Fraction', fontsize=14)
    ax2.set_title('3\' end', fontsize=14)

    legend = plt.gca().get_legend_handles_labels()
    labels, ids = np.unique(legend[1], return_index=True)
    handles = [legend[0][i] for i in ids]

    ax2.legend(handles, labels, bbox_to_anchor=(1.01, 0.5), loc="center left", frameon=False)

def get_mismatches(bam_file, bam_file_out=None):
    bam_data = read_bam(bam_file)

    header_str = bam_data['header']
    sequence_strings = bam_data['sequence_strings']

    res = []
    out_seqs = []
    for seq_string in sequence_strings:
        sequence = seq_string.split('\t')
        pos = sequence[3]
        num_mismatches = int(sequence[[idx for idx, s in enumerate(sequence) if 'NM:' in s][0]].split(':')[-1])
        if num_mismatches > 10:
            res.append(sequence)
        elif bam_file_out:
            out_seqs.append(seq_string)

    if bam_file_out:
        out_seq_str = '\n'.join(out_seqs) + '\n'
        bam_outstr = header_str + out_seq_str
        samtools_str = "samtools view -bh -o {} -".format(bam_file_out)
        subprocess.run(samtools_str, input=bam_outstr, shell=True, capture_output=True, text=True)

    return res

def get_indels_at_positions(positions, bam_file):
    res = []
    for index, row in positions.iterrows():
        pos = row['POS']

        samtools_str = "samtools view {} Y:{}-{}".format(bam_file, pos, pos)

        samtools_output = subprocess.run(samtools_str, shell=True, capture_output=True, text=True)

        output = samtools_output.stdout.split("\n")[:-1]
        sequences = [i.split('\t') for i in output]

        indel_count = 0
        indels = []

        seq_count = len(sequences)
        for seq in sequences:
            sequence = {'CIGAR' : seq[5],
                        'SEQ'   : seq[9],
                        'FLAG'  : int(seq[1]),
                        'MD'    : seq[[idx for idx, s in enumerate(seq) if 'MD:' in s][0]].split(':')[2],
                        'RG'    : seq[[idx for idx, s in enumerate(seq) if 'RG:' in s][0]].split(':')[2]}

            if ('I' in sequence['CIGAR']) or ('D' in sequence['CIGAR']):
                indel_count += 1
                indels.append([*sequence.values()])

        if indel_count:
            res.append([pos, seq_count, indel_count, round(indel_count/seq_count, 3), indels])

    return res
