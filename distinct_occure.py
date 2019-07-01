import pandas as pd
from distinct_occure_helpers import concat_ints, type_of_mutation, concat_alg, seq_type, if_complex, \
    subst_type, from_end, from_start
import warnings


warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None


def dist_occur(output_folder, localization_file, coordinates_file):

    all_mutations_raw = pd.read_csv(output_folder + '/all_mutations_filtered.csv')
    all_mutations_raw.columns = all_mutations_raw.columns.str.lower()
    all_mutations_raw.drop(['control:mut/norm', 'tumor:mut/norm', 'ratio', 'eval',
                            'qual',
                            'filter', 'info', 'format', 'normal', 'tumor',
                            'indiv_id', 'sample_id_tumor_name',
                            'sample_id_tumor_aliq', 'sample_id_normal_name',
                            'sample_id_normal_aliq'
                            ], axis=1, inplace=True)
    if all_mutations_raw.shape[0] == 0:
        print('no mutations found')
        return 0
    all_mutations_raw['mutation_type'] = all_mutations_raw.apply(lambda x: type_of_mutation(x), axis=1)

    all_mutations_raw.fillna(-1, inplace=True)

    all_mutations = all_mutations_raw.groupby(['chrom', 'pos',
                                               'indiv_name', 'ref', 'alt',
                                               'mutation_type']).agg({'alg': concat_alg,
                                                                      'norm_ref_count': sum,
                                                                      'norm_alt_count': sum,
                                                                      'tumor_ref_count': sum,
                                                                      'tumor_alt_count': sum
                                                                      }).reset_index()
    all_mutations['type_of_subst'] = all_mutations.apply(lambda x: subst_type(x), axis=1)

    all_mutations.to_csv(output_folder + '/all_mutations_algorithms_merged.csv',
                         sep=',',
                         index=False)

    all_mutations_not_mirna = all_mutations.copy()
    coordinates = pd.read_csv(coordinates_file, names=['chrom', 'start', 'stop', 'name'],
                              sep='\t')

    coordinates_not_mirna = coordinates[~coordinates['name'].str.contains('hsa-')]

    all_mutations_not_mirna = all_mutations_not_mirna.join(coordinates_not_mirna.set_index('chrom'), on='chrom', how='left')
    all_mutations_not_mirna = all_mutations_not_mirna[(all_mutations_not_mirna['pos'] >= all_mutations_not_mirna['start'])
                                  & (all_mutations_not_mirna['pos'] <= all_mutations_not_mirna['stop'])]
    all_mutations_not_mirna.to_csv(output_folder + '/all_mutations_not_mirna.csv',
                                   sep=',',
                                   index=False)

    localizations = pd.read_csv(localization_file, sep=',')
    all_mutations = all_mutations.join(localizations.set_index('chrom'), on='chrom', how='left')
    all_mutations = all_mutations[(all_mutations['pos'] >= all_mutations['start'])
                                  & (all_mutations['pos'] <= all_mutations['stop'])]
    double_check = all_mutations.groupby(['chrom', 'pos', 'indiv_name'])[['ref']].count()
    double_check.columns = ['no_of_loc']
    all_mutations = all_mutations.join(double_check, on=['chrom', 'pos', 'indiv_name'], how='left')

    all_mutations['seq_type'] = all_mutations['name'].apply(lambda x: seq_type(x))

    all_mutations['from_start'] = all_mutations.apply(lambda x: from_start(x, 'start', 'stop'), axis=1)
    all_mutations['from end'] = all_mutations.apply(lambda x: from_end(x, 'stop', 'start'), axis=1)

    all_mutations.to_csv(output_folder + '/all_mutations.csv',
                         sep=',',
                         index=False)

    df_complex = all_mutations.groupby(['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'indiv_name']).agg({
        'pos': ['nunique', 'count'],
        'no_of_loc': concat_ints
    }).reset_index()
    df_complex.columns = ['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'indiv_name',
                          'pos_nunique', 'pos_count', 'no_of_loc']
    df_complex['complex'] = df_complex['pos_count'].apply(lambda x: 0 if x < 2 else 1)

    df_complex.to_csv(output_folder + '/complex.csv',
                      sep=',',
                      index=False)

    df_by_gene = all_mutations.groupby(['chrom', 'pre_name', 'id', 'start_pre', 'seq_type']).agg(
        {'indiv_name': ['nunique',
                        'count'],
         'pos': 'nunique',
         'no_of_loc': concat_ints
         }).reset_index()
    df_by_gene.columns = ['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'indiv_name_nunique',
                          'indiv_name_count', 'pos_nunique', 'no_of_loc']
    df_by_gene['if_complex'] = df_by_gene.apply(lambda x: if_complex(x, df_complex), axis=1)
    df_by_gene.to_csv(output_folder + '/occur.csv',
                      sep=',',
                      index=False)

    df_by_mutation = all_mutations.groupby(['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'pos', 'ref', 'alt',
                                            'mutation_type', 'type_of_subst']).agg(
        {'indiv_name': 'nunique',
         'norm_ref_count': [sum, concat_ints],
         'norm_alt_count': [sum, concat_ints],
         'tumor_ref_count': [sum, concat_ints],
         'tumor_alt_count': [sum, concat_ints],
         'no_of_loc': concat_ints
         }).reset_index()
    df_by_mutation['chrom'] = df_by_mutation[('chrom', '')]
    df_by_mutation['name'] = df_by_mutation[('pre_name', '')]
    df_by_mutation['id'] = df_by_mutation[('id', '')]
    df_by_mutation['start_pre'] = df_by_mutation[('start_pre', '')]
    df_by_mutation['seq_type'] = df_by_mutation[('seq_type', '')]
    df_by_mutation['pos'] = df_by_mutation[('pos', '')]
    df_by_mutation['ref'] = df_by_mutation[('ref', '')]
    df_by_mutation['alt'] = df_by_mutation[('alt', '')]
    df_by_mutation['mutation_type'] = df_by_mutation[('mutation_type', '')]
    df_by_mutation['type_of_subst'] = df_by_mutation[('type_of_subst', '')]
    df_by_mutation['indiv_name'] = df_by_mutation[('indiv_name', 'nunique')]
    df_by_mutation['norm_ref_count_sum'] = df_by_mutation[('norm_ref_count', 'sum')]
    df_by_mutation['norm_ref_count_concat'] = df_by_mutation[('norm_ref_count', 'concat_ints')]
    df_by_mutation['norm_alt_count_sum'] = df_by_mutation[('norm_alt_count', 'sum')]
    df_by_mutation['norm_alt_count_concat'] = df_by_mutation[('norm_alt_count', 'concat_ints')]
    df_by_mutation['tumor_ref_count_sum'] = df_by_mutation[('tumor_ref_count', 'sum')]
    df_by_mutation['tumor_ref_count_concat'] = df_by_mutation[('tumor_ref_count', 'concat_ints')]
    df_by_mutation['tumor_alt_count_sum'] = df_by_mutation[('tumor_alt_count', 'sum')]
    df_by_mutation['tumor_alt_count_concat'] = df_by_mutation[('tumor_alt_count', 'concat_ints')]
    df_by_mutation['no_of_loc'] = df_by_mutation[('no_of_loc', 'concat_ints')]
    df_by_mutation.drop([x for x in df_by_mutation.columns if type(x) != str], axis=1, inplace=True)

    df_by_mutation.to_csv(output_folder + '/by_distinct_mutation.csv',
                          sep=',',
                          index=False)

    return 0
