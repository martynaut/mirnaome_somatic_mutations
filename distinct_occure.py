import pandas as pd
import numpy as np
from distinct_occure_helpers import concat_ints, type_of_mutation, take_from_coord, concat_alg, seq_type, if_complex, \
    subst_type, find_localization, from_end, from_start, find_in_mirna, find_arm, set_balance
import warnings


warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None


def dist_occur(output_folder, coordinates, confidence_file, mirgenedb_file, cancer_exon_file,
               localization_file, mirna_reads_file):

    confidence = pd.read_excel(confidence_file)
    mirgenedb = pd.read_csv(mirgenedb_file, sep='\t', skiprows=3,
                            names=['chr', '.', 'type', 'start', 'stop',
                                   '.1', 'orientation', '.2', 'ID'])
    mirgenedb['mirbase_id'] = mirgenedb['ID'].str.split(';').apply(lambda x: x[-1]).str.extract(
        r'Alias=([A-Z0-9]*)', expand=False
    )
    mirgenedb = mirgenedb[~mirgenedb['ID'].str.contains('v2')]
    mirgenedb = mirgenedb[~mirgenedb['ID'].str.contains('v3')]
    confidence = confidence.join(mirgenedb.set_index('mirbase_id')[['ID']], on='id', rsuffix='_mirdb')
    confidence['ID'].fillna('no', inplace=True)
    confidence['ID'] = confidence['ID'].str.replace(';', ':')
    confidence.drop(['score', 'id',
                     'start', 'stop'], inplace=True, axis=1)
    confidence.set_index('name', inplace=True)
    coordinates = pd.read_table(coordinates,
                                names=['chr', 'start', 'stop', 'gene'])
    coordinates['start_ref'] = coordinates['start']
    coordinates['stop_ref'] = coordinates['stop']
    cancer_exons = []
    with open(cancer_exon_file, 'r') as exon_file:
        for line in exon_file.readlines():
            cancer_exons.append(line[:-1])

    df_rows = pd.read_csv(output_folder + '/all_mutations_filtered.csv')
    df_rows.columns = df_rows.columns.str.lower()
    df_rows.drop(['control:mut/norm', 'tumor:mut/norm', 'ratio', 'eval'], axis=1, inplace=True)
    df_rows['mutation_type'] = df_rows.apply(lambda x: type_of_mutation(x), axis=1)
    df_rows['start'] = df_rows.apply(lambda x: take_from_coord(coordinates, 'start', x), axis=1)
    df_rows['stop'] = df_rows.apply(lambda x: take_from_coord(coordinates, 'stop', x), axis=1)
    df_rows['gene'] = df_rows.apply(lambda x: take_from_coord(coordinates, 'gene', x), axis=1)
    df_rows.fillna(-1, inplace=True)
    df_rows.drop(['qual',
                  'filter', 'info', 'format', 'normal', 'tumor',
                  'indiv_id', 'sample_id_tumor_name',
                  'sample_id_tumor_aliq', 'sample_id_normal_name',
                  'sample_id_normal_aliq'], axis=1, inplace=True)

    df_rows['seq_type'] = df_rows['gene'].apply(lambda x: seq_type(x, cancer_exons))

    df_rows3 = df_rows.groupby(['chrom', 'start', 'stop', 'gene', 'seq_type', 'pos',
                                'indiv_name', 'ref', 'alt',
                                'mutation_type']).agg({'alg': concat_alg,
                                                       'norm_ref_count': sum,
                                                       'norm_alt_count': sum,
                                                       'tumor_ref_count': sum,
                                                       'tumor_alt_count': sum
                                                       }).reset_index()
    df_rows3['type_of_subst'] = df_rows3.apply(lambda x: subst_type(x), axis=1)
    df_rows3['ones'] = 1
    df_rows3_conf = df_rows3.join(confidence, on='gene', how='left')
    df_rows3_conf.to_csv(output_folder + '/all_mutations_filtered_merge_programs.csv',
                         sep=',',
                         index=False)

    df_complex = df_rows3.groupby(['chrom', 'start', 'stop', 'gene', 'seq_type', 'indiv_name']).agg({
                                                                          'pos': ['nunique', 'count'],
                                                                          }).reset_index()
    df_complex.columns = ['chrom', 'start', 'stop', 'gene', 'seq_type', 'indiv_name',
                          'pos_nunique', 'pos_count']
    df_complex['complex'] = df_complex['pos_count'].apply(lambda x: 0 if x < 2 else 1)

    df_complex['ones'] = 1
    df_complex = df_complex.join(confidence, on='gene', how='left')
    df_complex.to_csv(output_folder + '/complex.csv',
                      sep=',',
                      index=False)

    df_chrom = df_complex[(df_complex['seq_type'] == 'mirna')]
    df_chrom = df_chrom.groupby(['chrom']).nunique()
    df_chrom[['chrom', 'gene', 'indiv_name']].to_csv(output_folder + '/miRNA_per_chrom.csv',
                                                     sep=',')

    df_by_gene = df_rows3.groupby(['chrom', 'start', 'stop', 'gene', 'seq_type']).agg({'indiv_name': ['nunique',
                                                                                                      'count'],
                                                                                       'pos': 'nunique',
                                                                                       }).reset_index()
    df_by_gene['ones'] = 1
    df_by_gene.columns = ['chrom', 'start', 'stop', 'gene', 'seq_type', 'indiv_name_nunique',
                          'indiv_name_count', 'pos_nunique', 'ones']
    df_by_gene['if_complex'] = df_by_gene.apply(lambda x: if_complex(x, df_complex), axis=1)
    df_by_gene = df_by_gene.join(confidence, on='gene', how='left')
    df_by_gene.to_csv(output_folder + '/occur.csv',
                      sep=',',
                      index=False)

    df_by_mutation = df_rows3.groupby(['chrom', 'start', 'stop', 'gene', 'seq_type', 'pos', 'ref', 'alt',
                                       'mutation_type', 'type_of_subst']).agg({'indiv_name': 'nunique',
                                                                               'norm_ref_count': [sum, concat_ints],
                                                                               'norm_alt_count': [sum, concat_ints],
                                                                               'tumor_ref_count': [sum, concat_ints],
                                                                               'tumor_alt_count': [sum, concat_ints]
                                                                               }).reset_index()

    df_by_mutation = df_by_mutation.join(confidence, on='gene', how='left')

    del confidence
    del mirgenedb
    del coordinates
    del cancer_exons
    del df_rows3
    del df_complex
    del df_rows
    del df_chrom
    del df_by_gene

    df_by_mutation['chrom'] = df_by_mutation[('chrom', '')]
    df_by_mutation['start'] = df_by_mutation[('start', '')]
    df_by_mutation['stop'] = df_by_mutation[('stop', '')]
    df_by_mutation['gene'] = df_by_mutation[('gene', '')]
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
    df_by_mutation.drop([x for x in df_by_mutation.columns if type(x) != str], axis=1, inplace=True)

    distinct = df_by_mutation.copy().reset_index()
    all_mutations = df_rows3_conf.copy().reset_index()

    del df_by_mutation

    localizations = pd.read_csv(localization_file, sep=',')

    mirnas = localizations[localizations['type'].isin(['pre-seed',
                                                       'seed',
                                                       'post-seed'])].copy()
    mirnas['name'] = mirnas['name'].str.extract(r'([a-z0-9\-]*)_[post\-sedr]*', expand=False)
    mirnas = mirnas.groupby(['chrom', 'name']).agg({'start': min, 'stop': max, 'orientation': 'first'}).reset_index()

    # take only miRNAs
    distinct = distinct[distinct['seq_type'] == 'mirna']
    all_mutations = all_mutations[all_mutations['seq_type'] == 'mirna']

    # check which miRNAs we don't have strand information
    without_strand = distinct[distinct['Strand'].isnull()]['gene'].unique()
    for mirna in without_strand:
        print('\n' + mirna + ' - localization is not included as strand is not known\n')

    del without_strand

    distinct['name'] = distinct['gene'].str.lower().str.extract(r'([hsa]*\-[mirlet]*\-[0-9a-z]*)', expand=False)
    all_mutations['name'] = all_mutations['gene'].str.lower().str.extract(r'([hsa]*\-[mirlet]*\-[0-9a-z]*)',
                                                                          expand=False)
    print(distinct.shape)
    if distinct.shape[0] > 0:
        newcols = distinct.apply(lambda x: pd.Series(find_localization(x, localizations)), axis=1)
    else:
        newcols = pd.DataFrame([[np.nan,
                   np.nan,
                   np.nan,
                   np.nan,
                   np.nan,
                   np.nan,
                   np.nan,
                   np.nan]], columns=['chrom_loc',
                       'name_loc',
                       'start_loc',
                       'stop_loc',
                       'orient_loc',
                       'based_on_coordinates',
                       'arm',
                       'type'])
    print(distinct.columns)
    print(newcols.columns)
    newcols.columns = ['chrom_loc',
                       'name_loc',
                       'start_loc',
                       'stop_loc',
                       'orient_loc',
                       'based_on_coordinates',
                       'arm',
                       'type']
    distinct = distinct.join(newcols)

    del newcols

    try:

        distinct['from_start'] = distinct.apply(lambda x: from_start(x, 'start_loc', 'stop_loc'), axis=1)
        distinct['from end'] = distinct.apply(lambda x: from_end(x, 'stop_loc', 'start_loc'), axis=1)

        new_cols2 = distinct.apply(lambda x: pd.Series(find_in_mirna(x, mirnas)), axis=1)
        new_cols2.columns = ['start_mirna', 'stop_mirna']
    except ValueError:
        new_cols2 = pd.DataFrame([[np.nan,
                   np.nan,
                   np.nan,
                   np.nan]], columns=['from_start',
                       'from_end',
                       'start_mirna',
                       'srop_mirna'])

    distinct = distinct.join(new_cols2)

    try:

        newcols_all = all_mutations.apply(lambda x: pd.Series(find_localization(x, localizations)), axis=1)
        newcols_all.columns = ['chrom_loc',
                           'name_loc',
                           'start_loc',
                           'stop_loc',
                           'orient_loc',
                           'based_on_coordinates',
                           'arm',
                           'type']
    except ValueError:
        newcols_all = pd.DataFrame([[np.nan,
                       np.nan,
                       np.nan,
                       np.nan,
                       np.nan,
                       np.nan,
                       np.nan,
                       np.nan]], columns=['chrom_loc',
                                          'name_loc',
                                          'start_loc',
                                          'stop_loc',
                                          'orient_loc',
                                          'based_on_coordinates',
                                          'arm',
                                          'type'])

    all_mutations = all_mutations.join(newcols_all)

    del newcols_all

    try:

        all_mutations['from_start'] = all_mutations.apply(lambda x: from_start(x, 'start_loc', 'stop_loc'), axis=1)
        all_mutations['from end'] = all_mutations.apply(lambda x: from_end(x, 'stop_loc', 'start_loc'), axis=1)
        new_cols2 = all_mutations.apply(lambda x: pd.Series(find_in_mirna(x, mirnas)), axis=1)
        new_cols2.columns = ['start_mirna', 'stop_mirna']

    except ValueError:
        new_cols2 = pd.DataFrame([[np.nan,
                   np.nan,
                   np.nan,
                   np.nan]], columns=['from_start',
                       'from_end',
                       'start_mirna',
                       'srop_mirna'])

    del mirnas

    all_mutations = all_mutations.join(new_cols2)

    del new_cols2

    mirna_reads = pd.read_csv(mirna_reads_file)
    mirna_reads['gene'] = mirna_reads['Name_pre'].str.extract(r'([a-z]*\-[a-z]*\-[a-z0-9]*)', expand=False)
    mirna_reads['arm'] = mirna_reads.apply(lambda x: find_arm(x), axis=1)

    mirna_reads.drop(['start_pre', 'stop_pre', '-/+_pre', 'start', 'stop', 'ID', 'Alias', ],
                     axis=1, inplace=True)
    mirna_reads_3p = mirna_reads[mirna_reads['arm'] == '3p'].copy().groupby(['chr', 'From', 'Name_pre']).first()

    mirna_reads_5p = mirna_reads[mirna_reads['arm'] == '5p'].copy().groupby(['chr', 'From', 'Name_pre']).first()

    mirna_reads_joined = mirna_reads_3p.join(mirna_reads_5p, rsuffix='_5p', lsuffix='_3p', how='outer')
    mirna_reads_joined['ratio 3/(3+5)'] = mirna_reads_joined['reads_3p'] / (mirna_reads_joined['reads_3p'] +
                                                                            mirna_reads_joined['reads_5p'])
    mirna_reads_joined['ratio 5/(3+5)'] = mirna_reads_joined['reads_5p'] / (mirna_reads_joined['reads_3p'] +
                                                                            mirna_reads_joined['reads_5p'])
    # set balance proportion btw 5p i 3p miRNA
    mirna_reads_joined['balance'] = mirna_reads_joined.apply(lambda x: set_balance(x, 0.75), axis=1)
    mirna_reads_joined = mirna_reads_joined.reset_index()[['chr',
                                                           'Name_pre',
                                                           'balance']].copy().set_index(
        [
            'chr',
            'Name_pre'
        ])

    distinct['name'] = distinct['gene'].str.lower().str.extract(r'([hsa]*\-[mirlet]*\-[0-9]*)', expand=False)
    distinct = distinct.join(mirna_reads_joined, on=['chrom', 'gene'])

    all_mutations['name'] = all_mutations['gene'].str.lower().str.extract(r'([hsa]*\-[mirlet]*\-[0-9]*)', expand=False)
    all_mutations = all_mutations.join(mirna_reads_joined, on=['chrom', 'gene'])

    # save files with miRNAs with localization

    all_mutations.to_csv(output_folder + '/all_mutations_filtered_mut_type_gene.csv',
                         sep=',',
                         index=False)
    distinct.to_csv(output_folder + '/distinct_with_loc.csv',
                    sep=',',
                    index=False)

    return 0
