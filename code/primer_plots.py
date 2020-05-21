import pandas as pd


def file_process(file, primer, primer_seq):
    df = pd.read_csv(file)

    # select the highest scoring alignment for all seqs queried
    df1 = df.sort_values('bit_score').drop_duplicates('subject acc.ver', keep='last')

    output_dict = {'Primer': primer,
                   'Minimum % identity': df1['%_identity'].min(), 'Max mismatches': df1['mismatches'].max(),
                   'Unique sequences': len(df1['subject acc.ver'].unique()), 'Total alignments': len(df1),
                   'Number of alignments with mismatches': len(df1.loc[df1['%_identity'] < 100, :]),
                   '% alignments with 100% identity': 100 * (1 - len(df1.loc[df1['%_identity'] < 100, :]) / len(df1)),
                   'Mean Alignment Length': df1['alignment_length'].mean(),
                   '% Sequences with Alignment Length Equal to full length of primer':
                       (100 * (1 - len(df1.loc[df1['alignment_length'] < len(primer_seq), :]) / len(df1))),
                   'Minimum Alignment Length': df1['alignment_length'].min()}

    # identify accession numbers for mismatches
    mismatch_seqs = df1.loc[df1['%_identity'] < 100, :]
    mismatch_list = mismatch_seqs['subject acc.ver'].to_list()

    # add the mismatch list to the df after transposing
    df_out = pd.DataFrame.from_dict(output_dict, orient='index').T
    df_out['mismatch_accessions'] = [mismatch_list]

    # identify accession numbers for short seqs
    short_seqs = df1.loc[df1['alignment_length'] < len(primer_seq), :]
    short_list = short_seqs['subject acc.ver'].to_list()

    df_out['short_accessions'] = [short_list]
    df_out['all_problem_acessions'] = [list(set(mismatch_list + short_list))]

    return df_out


def problem_seq_process(df_subset, string_list):
    problem_dict = {}

    for string in string_list:
        problem_df = df_subset[df_subset['Primer'].str.contains(string)]
        problem_out_list = list(set([a for b in problem_df.all_problem_acessions.tolist() for a in b]))
        problem_dict[string] = problem_out_list

    problem_out_list = list(set(problem_dict[string_list[0]]) & set(problem_dict[string_list[1]]))

    return problem_out_list


def process_results(in_data_list):
    results_list = []

    for i_list in in_data_list:
        out = file_process(i_list[0], i_list[1], i_list[2])
        results_list.append(out)

    df_append = pd.concat(results_list).reset_index(drop=True)
    df_append.to_csv('../data/primer_BLAST_summary.csv', index=False)
    problem_list = problem_seq_process(df_append, ['N1', 'N2'])

    return problem_list


# process the BLAST alignment results
data_list = [['../data/BZSJCU50114-Alignment-HitTable_N1_Forward.csv', 'N1 Forward', 'GACCCCAAAATCAGCGAAAT'],
             ['../data/BZT5F9NN114-Alignment-HitTable_N1_Reverse.csv', 'N1 Reverse', 'TCTGGTTACTGCCAGTTGAATCTG'],
             ['../data/BZTBUZH7114-Alignment-HitTable_N1_Probe.csv', 'N1 Probe', 'ACCCCGCATTACGTTTGGTGGACC'],
             ['../data/BZTHE240114-Alignment-HitTable_N2_Forward.csv', 'N2 Forward', 'TTACAAACATTGGCCGCAAA'],
             ['../data/BZU005YT114-Alignment-HitTable_N2_Reverse.csv', 'N2 Reverse', 'GCGCGACATTCCGAAGAA'],
             ['../data/BZU4JNT2114-Alignment-HitTable_N2_Probe.csv', 'N2 Probe', 'ACAATTTGCCCCCAGCGCTTCAG']]

problem_accessions = process_results(data_list)

print(problem_accessions)
