import pandas as pd


def file_process(file, primer, target_length):
    df = pd.read_csv(file)

    # select the highest scoring alignment for all seqs queried
    df1 = df.sort_values('bit_score').drop_duplicates('subject acc.ver', keep='last')

    output_dict = {'Primer': primer,
                   'Minimum % identity': df1['%_identity'].min(), 'Max mismatches': df1['mismatches'].max(),
                   'Unique sequences': len(df1['subject acc.ver'].unique()), 'Total alignments': len(df1),
                   'Number of alignments with mismatches': len(df1.loc[df1['%_identity'] < 100, :]),
                   '% alignments with 100% identity': 100 * (1 -len(df1.loc[df1['%_identity'] < 100, :]) / len(df1)),
                   'Mean Alignment Length': df1['alignment_length'].mean(),
                   '% Sequences with Alignment Length Equal to full length of primer': (100 * (1 - len(df1.loc[df1['alignment_length'] < target_length, :]) / len(df1))),
                   'Minimum Alignment Length': df1['alignment_length'].min()}

    # identify accession numbers for mismatches
    mismatch_seqs = df1.loc[df1['%_identity'] < 100, :]
    mismatch_list = mismatch_seqs['subject acc.ver'].to_list()

    df_out = pd.DataFrame.from_dict(output_dict, orient='index').T
    df_out['mismatch_accessions'] = [mismatch_list]

    # identify accession numbers for short seqs
    mismatch_seqs = df1.loc[df1['alignment_length'] < target_length, :]
    print(mismatch_seqs)

    return df_out


output_N1_forward = file_process('../data/BZSJCU50114-Alignment-HitTable_N1_Forward.csv', 'N1 Forward', 20)
output_N1_reverse = file_process('../data/BZT5F9NN114-Alignment-HitTable_N1_Reverse.csv', 'N1 Reverse', 24)
output_N1_probe = file_process('../data/BZTBUZH7114-Alignment-HitTable_N1_Probe.csv', 'N1 Probe', 24)
output_N2_forward = file_process('../data/BZTHE240114-Alignment-HitTable_N2_Forward.csv', 'N2 Forward', 20)
output_N2_reverse = file_process('../data/BZU005YT114-Alignment-HitTable_N2_Reverse.csv', 'N2 Reverse', 18)
output_N2_probe = file_process('../data/BZU4JNT2114-Alignment-HitTable_N2_Probe.csv', 'N2 Probe', 23)


df_append = pd.concat([output_N1_forward, output_N1_reverse, output_N1_probe, output_N2_forward, output_N2_reverse,
                         output_N2_probe])

print(df_append.reset_index(drop=True))

df_append.to_csv('../data/primer_BLAST_summary.csv', index=False)

