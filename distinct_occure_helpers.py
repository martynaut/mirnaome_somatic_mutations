import numpy as np


def set_balance(row, ratio):
    if row['ratio 3/(3+5)'] >= ratio:
        return '3p'
    elif row['ratio 5/(3+5)'] >= ratio:
        return '5p'
    elif np.isnan(row['reads_3p']) and np.isnan(row['reads_5p']):
        return 'unknown'
    elif np.isnan(row['reads_3p']):
        return '5p'
    elif np.isnan(row['reads_5p']):
        return '3p'
    else:
        return 'both'


def find_in_mirna(row, df_loc):

    if df_loc[
              (df_loc['chrom'] == row['chrom']) &
              (df_loc['start'] <= row['pos']) &
              (df_loc['orientation'] == row['orient_loc']) &
              (df_loc['stop'] >= row['pos'])].shape[0] != 0:

        temp = df_loc[
                         (df_loc['chrom'] == row['chrom']) &
                         (df_loc['start'] <= row['pos']) &
                         (df_loc['orientation'] == row['orient_loc']) &
                         (df_loc['stop'] >= row['pos'])].values[0]
        if row['orient_loc'] == '+':
            start = row['pos'] - temp[2] + 1
            stop = row['pos'] - temp[3] - 1
        else:
            start = -(row['pos'] - temp[3] - 1)
            stop = -(row['pos'] - temp[2] + 1)
        localiz = [start, stop]
    else:
        localiz = [np.nan,
                   np.nan]
    return localiz


def find_arm(row):
    if row['-/+'] == '+':
        if row['start'] - row['start_pre'] < row['stop_pre'] - row['stop']:
            return '5p'
        else:
            return '3p'
    if row['-/+'] == '-':
        if row['start'] - row['start_pre'] < row['stop_pre'] - row['stop']:
            return '3p'
        else:
            return '5p'


def from_start(row, column_start, column_stop):
    if row['orientation'] == '+':
        return row['pos'] - row[column_start] + 1
    else:
        return row[column_stop] - row['pos'] + 1


def from_end(row, column_stop, column_start):
    if row['orientation'] == '+':
        return row['pos'] - row[column_stop] - 1
    else:
        return row[column_start] - row['pos'] - 1


def find_localization(row, df_loc):
    print(row)
    localizations = df_loc[(df_loc['chrom'] == row['chrom']) &
                         (df_loc['start'] <= row['pos']) &
                         (df_loc['stop'] >= row['pos'])].values
    print(localizations)
    if len(localizations) > 1:
        raise ValueError
    else:
        return localizations[0]


def if_complex(row, complex_df):
    if complex_df[(complex_df['chrom'] == row['chrom']) &
                  (complex_df['pre_name'] == row['pre_name']) &
                  (complex_df['id'] == row['id']) &
                  (complex_df['start_pre'] == row['start_pre']) &
                  (complex_df['seq_type'] == row['seq_type'])].shape[0] != 0:
        values = complex_df[(complex_df['chrom'] == row['chrom']) &
                            (complex_df['pre_name'] == row['pre_name']) &
                            (complex_df['id'] == row['id']) &
                            (complex_df['start_pre'] == row['start_pre']) &
                            (complex_df['seq_type'] == row['seq_type'])]['complex'].unique()
        if 1 in values:
            return 1
        else:
            return 0
    else:
        return 0


def concat_ints(col):
    row = list(col.values)
    new_row = []
    for x in row:
        new_row.append(str(x))
    return '"' + ':'.join(new_row) + '"'


def concat_alg(col):
    row = list(col.values)
    new_row = []
    for x in row:
        new_row.append(str(x))
    new_row = sorted(set(new_row))
    return '"' + ':'.join(new_row) + '"'


def type_of_mutation(row):
    if len(row['ref']) > len(row['alt']):
        return 'del'
    elif len(row['ref']) == len(row['alt']):
        return 'subst'
    elif ',' in row['alt']:
        return 'subst'
    else:
        return 'ins'


def take_from_coord(coordinates, column_name, row):
    return coordinates[(coordinates['chr'] == row['chrom']) &
                       (coordinates['start'] < int(row['pos'])) &
                       (coordinates['stop'] > int(row['pos']))][column_name].values[0]


def seq_type(value):
    if 'hsa-' in value:
        return 'mirna'
    else:
        return 'not_defined'


def subst_type(row):
    if row['mutation_type'] == 'subst':
        if (((row['ref'] in ['A', 'G']) and (row['alt'] in ['A', 'G'])) or
                ((row['ref'] in ['C', 'T']) and (row['alt'] in ['C', 'T']))):
            return 'transition'
        else:
            return 'transversion'
    else:
        return 'n.a.'
