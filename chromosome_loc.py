import pandas as pd


def chrom_localisation(folder_output, chromosome_file):

    mir_coverage = pd.read_csv(chromosome_file, sep='\t')
    mir_coverage['TargetID'] = mir_coverage['TargetID'].apply(lambda x: x.lower())
    mir_coverage = mir_coverage[mir_coverage['HighCoverage'] == 1]
    mir_coverage['chrom'] = mir_coverage['Interval'].str.split(':').str[0]
    chrom_counter = mir_coverage.groupby('chrom')['TargetID'].count()
    chrom_counter.columns = ['high_confidence_mirnas']
    chrom_counter.to_csv(folder_output + '/high_confidence_mirna_per_chrom.csv',
                         sep=','
                         )
    mir_mutations = pd.read_csv(folder_output + '/all_mutations_filtered_merge_programs.csv',
                                sep=',')
    mir_mutations = mir_mutations.join(mir_coverage.set_index('TargetID'), on='gene', how='left', rsuffix='_cover')
    mir_mutations = mir_mutations[mir_mutations['HighCoverage'] == 1]
    mir_mutations_mirna = mir_mutations.groupby('chrom_cover').agg({'gene': 'nunique'})
    mir_mutations_mirna.to_csv(folder_output + '/mirnas_perchrom.csv',
                               sep=','
                               )
    mir_mutations_pat = mir_mutations.groupby(['chrom_cover', 'indiv_name']).agg({'start': 'sum'})
    mir_mutations_pat.columns = ['count']
    mir_mutations_pat['count'] = 1
    mir_mutations_pat = mir_mutations_pat.groupby(['chrom_cover']).agg({'count': 'sum'})
    mir_mutations_pat.to_csv(folder_output + '/patients_perchrom.csv',
                             sep=','
                             )
