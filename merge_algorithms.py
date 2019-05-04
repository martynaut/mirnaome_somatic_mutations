import pandas as pd
import numpy as np


def filter_and_combine(output_folder):
    df_muse = pd.read_csv(output_folder + '/results_muse_eval.csv')
    df_mutect2 = pd.read_csv(output_folder + '/results_mutect2_eval.csv')
    df_varscan2 = pd.read_csv(output_folder + '/results_varscan2_eval.csv')
    df_ss = pd.read_csv(output_folder + '/results_somaticsniper_eval.csv')

    df_muse = df_muse[df_muse['eval']]
    df_mutect2 = df_mutect2[df_mutect2['eval']]
    df_varscan2 = df_varscan2[df_varscan2['eval']]
    df_ss = df_ss[df_ss['eval']]

    df_muse['alg'] = 'muse'
    df_mutect2['alg'] = 'mutect2'
    df_varscan2['alg'] = 'varscan2'
    df_ss['alg'] = 'somaticsniper'

    df = pd.concat([df_muse, df_mutect2, df_varscan2, df_ss])

    df = df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
             'TUMOR', 'NORMAL', 'indiv_name', 'indiv_id', 'sample_id_tumor_name',
             'sample_id_tumor_aliQ', 'sample_id_normal_name',
             'sample_id_normal_aliQ', 'norm_ref_count', 'norm_alt_count',
             'tumor_ref_count', 'tumor_alt_count', 'BQ_ref_tum', 'BQ_alt_tum',
             'BQ_ref_norm', 'BQ_alt_norm', 'QSS_ref_tum', 'QSS_alt_tum',
             'QSS_ref_nor', 'QSS_alt_nor', 'SSC', 'SPV', 'eval', 'name', 'alg']]

    df['control:mut/norm'] = df['norm_alt_count'] / df['norm_ref_count']
    df['tumor:mut/norm'] = df['tumor_alt_count'] / df['tumor_ref_count']

    df['ratio'] = df['tumor:mut/norm'] / df['control:mut/norm']

    df.replace(np.inf, 0, inplace=True)

    df.to_csv(output_folder + '/all_mutations_filtered.csv',
              sep=',',
              index=False)
    return 0
