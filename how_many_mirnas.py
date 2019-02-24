import pandas as pd


def check_covered_mirnas(output_folder, covered_mirna_file):

    df = pd.read_csv(covered_mirna_file, sep='\t')
    df = df.loc[df['HighCoverage'] == 1, ['TargetID', 'Interval']]
    df['chrom'] = df['Interval'].str.split(':').str[0]
    df['start'] = df['Interval'].str.split(':').str[1].str.split('-').str[0]
    df['stop'] = df['Interval'].str.split(':').str[1].str.split('-').str[1]
    df['TargetID'] = df['TargetID'].str.lower()
    df_results = pd.read_csv(output_folder + '/distinct.csv', sep=',')
    df_results = df_results[df_results['seq_type'] == 'mirna']
    df_results = df_results.groupby(['chrom',
                                     'gene'])[['start',
                                               'stop']].agg('first').reset_index()

    joined_df = df.join(df_results.set_index('gene'), how='outer', on='TargetID',
                        rsuffix='_results')
    outside_1600 = joined_df[joined_df['Interval'].isnull()]
    outside_1600.to_csv(output_folder + '/mirnas_ouside_1600.csv', sep=',', index=False)
