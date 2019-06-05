import pandas as pd
import RNA
from Bio import SeqIO
import re
from prepare_reads_file_helpers import parse_html
from distinct_occure_helpers import find_arm, set_balance


def add_mirgenedb(mirgenedb_file, output_folder):
    mirgenedb = pd.read_csv(mirgenedb_file, sep='\t', comment='#',
                            names=['chr', '.', 'type', 'start', 'stop',
                                   '.1', 'orientation', '.2', 'ID'])
    mirgenedb['mirbase_id'] = mirgenedb['ID'].str.extract(
        r'Alias=(MI[0-9]+)', expand=False
    )
    mirgenedb['ID'] = mirgenedb['ID'].str.extract(
        r'ID=([a-zA-Z0-9\-_]+)', expand=False
    )
    mirgenedb.dropna(subset=['mirbase_id'], inplace=True)

    mirgenedb = mirgenedb[~mirgenedb['ID'].str.contains('v2')]
    mirgenedb = mirgenedb[~mirgenedb['ID'].str.contains('v3')]
    mirgenedb.rename({'ID': 'mirgenedb_ID'}, axis='columns', inplace=True)
    mirgenedb.drop(['chr', '.', 'type', 'start', 'stop', '.1', 'orientation', '.2'],
                   axis=1, inplace=True)
    mirgenedb.drop_duplicates(inplace=True)
    mirgenedb.to_csv(output_folder + '/.mirgenedb.csv', index=False)


def merge_all(output_folder):
    confidence = pd.read_csv(output_folder + '/.confidence_file.csv')
    localizations = pd.read_csv(output_folder + '/.localizations.csv')
    print(localizations.shape)
    localizations['pre_name'] = localizations['name'].str.split('_').str[0]

    output_file = localizations.join(confidence.set_index(['mirbase_id',
                                                           'Strand', 'mirna_name']), on=['id', 'orientation',
                                                                                       'pre_name'], how='left')
    output_file.loc[output_file['confidence'].isnull(), 'confidence'] = 'Low'
    print(output_file.shape)
    reads = pd.read_csv(output_folder + '/.reads.csv')
    output_file = output_file.join(reads.set_index('From'), on='id', how='left')
    print(output_file.shape)
    mirgenedb = pd.read_csv(output_folder + '/.mirgenedb.csv')
    output_file = output_file.join(mirgenedb.set_index('mirbase_id'), on='id', how='left')
    print(output_file.shape)
    print("miRNA genome regions:")
    print(output_file.shape[0] / 9)
    output_file.to_csv(output_folder + '/localizations.csv',
                       index=False)


def create_reads(output_folder):
    loc_info = pd.read_csv(output_folder + '/.hsa_mirbase_coordinates.csv', sep=',')
    pre_mirbase = loc_info[loc_info['type'] == 'miRNA_primary_transcript'].copy()
    matures = loc_info[loc_info['type'] == 'miRNA'].copy()

    pre_mirbase.drop(['ID', 'type', 'From', 'chr'], inplace=True, axis=1)
    matures.drop(['type'], inplace=True, axis=1)

    joined_df = matures.join(pre_mirbase.set_index('Alias'), on='From', how='left', rsuffix='_pre')
    joined_df = joined_df[(joined_df['start'] >= joined_df['start_pre']) &
                          (joined_df['stop'] <= joined_df['stop_pre'])]
    joined_df['reads'] = joined_df.apply(lambda x: parse_html(x), axis=1)
    joined_df['arm'] = joined_df['Name'].str[-2:]
    joined_df.loc[~joined_df['arm'].isin(['3p', '5p']), 'arm'] = \
        joined_df[~joined_df['arm'].isin(['3p', '5p'])].apply(lambda x: find_arm(x), axis=1)

    joined_df.drop(['start_pre', 'stop_pre', '-/+_pre', 'start', 'stop', 'ID', 'Alias'],
                   axis=1, inplace=True)
    mirna_reads_3p = joined_df[joined_df['arm'] == '3p'].copy().groupby(['From']).first()

    mirna_reads_5p = joined_df[joined_df['arm'] == '5p'].copy().groupby(['From']).first()

    mirna_reads_joined = mirna_reads_3p.join(mirna_reads_5p, rsuffix='_5p', lsuffix='_3p', how='outer')
    mirna_reads_joined['ratio 3/(3+5)'] = mirna_reads_joined['reads_3p'] / (mirna_reads_joined['reads_3p'] +
                                                                            mirna_reads_joined['reads_5p'])
    mirna_reads_joined['ratio 5/(3+5)'] = mirna_reads_joined['reads_5p'] / (mirna_reads_joined['reads_3p'] +
                                                                            mirna_reads_joined['reads_5p'])

    mirna_reads_joined['balance'] = mirna_reads_joined.apply(lambda x: set_balance(x, 0.75), axis=1)
    mirna_reads_joined = mirna_reads_joined.reset_index()[['From',
                                                           'balance']]

    mirna_reads_joined.to_csv(output_folder + '/.reads.csv', sep=',', index=False)


def add_confidence(confidence_file, confidence_score_file,
                   aliases_file, chromosome_build_file, output_folder):
    df_confidence_file = pd.read_csv(confidence_file, sep='\t',
                                     names=['mirna_name', 'id', '1', '2', '3', '4', '5',
                                            '6', '7', '8', '9', '10',
                                            '11', '12', '13', '14', '15'])
    df_confidence_file = df_confidence_file[df_confidence_file['mirna_name'].str.contains('hsa')]
    df_confidence_score_file = pd.read_csv(confidence_score_file, sep='\t',
                                           names=['id', 'score'])

    df_aliases_file = pd.read_csv(aliases_file, sep='\t',
                                  names=['id', 'aliases'])
    df_aliases_file = df_aliases_file[~df_aliases_file['id'].str.contains('MIMAT')]
    df_aliases_file = df_aliases_file[df_aliases_file['aliases'].str.contains('hsa')]
    df_aliases_file['aliases'] = df_aliases_file['aliases'].str.split(';')
    df_aliases_file = df_aliases_file.aliases.apply(pd.Series) \
        .merge(df_aliases_file, left_index=True, right_index=True).\
        drop(["aliases"], axis=1).\
        melt(id_vars=['id'], value_name="alias").\
        drop("variable", axis=1) \
        .dropna()
    df_aliases_file = df_aliases_file[df_aliases_file['alias'] != '']
    df_aliases_file.rename({'id': 'mirbase_id'}, axis='columns', inplace=True)
    df_mirna_chromosome_build = pd.read_csv(chromosome_build_file, sep='\t',
                                            names=['id', '1', 'start', 'stop', 'Strand'])
    df_mirna_chromosome_build = df_mirna_chromosome_build[
        df_mirna_chromosome_build['1'].str.contains('chr')][['id', 'start', 'stop', 'Strand']]

    df_output = df_confidence_file.join(df_confidence_score_file.set_index('id'),
                                        on='id', how='inner')

    df_output = df_output.join(df_mirna_chromosome_build.set_index('id'), on='id', how='inner')
    df_output = df_output.join(df_aliases_file.set_index('alias'), on='mirna_name')

    df_output = df_output[['mirna_name', 'start', 'stop', 'score', 'Strand', 'mirbase_id']]
    df_output.rename({'start': 'start_pre_build',
                      'stop': 'stop_pre_build'}, axis='columns', inplace=True)
    df_output.drop_duplicates(inplace=True)
    df_output['confidence'] = df_output['score'].apply(lambda x: 'High' if x > 0 else 'Low')
    df_output.drop('score', axis=1, inplace=True)
    df_output.to_csv(output_folder + '/.confidence_file.csv',
                     index=False)


def prepare_structures(hairpin_fa_file, output_folder):
    df = pd.DataFrame(columns=['id', 'name', 'seq', 'structure'])
    fasta_sequences = SeqIO.parse(open(hairpin_fa_file), 'fasta')
    for fasta_sequence in fasta_sequences:
        if 'hsa' in fasta_sequence.id:
            name = re.search(' (MI[0-9_]*) ', fasta_sequence.description, re.IGNORECASE).group(1)
            df_temp = pd.DataFrame(data=[[fasta_sequence.id,
                                          name,
                                          str(fasta_sequence.seq),
                                          RNA.fold(str(fasta_sequence.seq))[0]]],
                                   columns=['id', 'name', 'seq', 'structure'])
            df = pd.concat([df, df_temp])
    df.reset_index(inplace=True)
    df.drop('index', axis=1, inplace=True)
    df.to_csv(output_folder + '/.hairpin_structure.csv', sep=',', index=False)


def prepare_hsa_files(hsa_gff_mirbase_file, output_folder):
    names = ['chr', '.', 'type', 'start', 'stop', '.2', '-/+', '.3', 'desc']
    data_chr = pd.read_csv(hsa_gff_mirbase_file, sep="\t", comment='#', names=names)
    data_chr['ID'] = data_chr['desc'].str.extract(r'ID\=([A-Za-z0-9]+)', expand=False)
    data_chr['Alias'] = data_chr['desc'].str.extract(r'Alias\=([A-Za-z0-9]+)', expand=False)
    data_chr['Name'] = data_chr['desc'].str.extract(r'Name\=([A-Za-z0-9_\-]+)', expand=False)
    data_chr['From'] = data_chr['desc'].str.extract(r'Derives_from\=([A-Za-z0-9]+)', expand=False)

    data_chr.drop(['desc', '.', '.2', '.3'], inplace=True, axis=1)
    data_chr.to_csv(output_folder + '/.hsa_mirbase_coordinates.csv', index=False)


def create_loc(output_folder):
    # localization based on mirbase
    loc_info = pd.read_csv(output_folder + '/.hsa_mirbase_coordinates.csv', sep=',')
    structures = pd.read_csv(output_folder + '/.hairpin_structure.csv', sep=',')

    pre_mirbase = loc_info[loc_info['type'] == 'miRNA_primary_transcript'].copy()
    matures = loc_info[loc_info['type'] == 'miRNA'].copy()

    pre_mirbase.drop(['ID', 'type', 'From', 'chr'], inplace=True, axis=1)

    matures.drop(['ID', 'type'], inplace=True, axis=1)
    matures = matures.merge(pre_mirbase, right_on='Alias', left_on='From', how='left', suffixes=('', '_pre'))

    matures = matures[(matures['start'] >= matures['start_pre']) &
                      (matures['stop'] <= matures['stop_pre'])]

    matures['arm'] = matures['Name'].str[-2:]
    matures.loc[~matures['arm'].isin(['3p',
                                      '5p']), 'arm'] = \
        matures[~matures['arm'].isin(['3p',
                                      '5p'])].apply(lambda x: find_arm(x), axis=1)

    matures['Name'] = matures['Name_pre'].map(str) + matures['Name'].str.extract(
        r'(\_[35]{1,1}p)', expand=False).fillna('')
    joined_new = pre_mirbase.join(matures.set_index('From'), on='Alias', how='left', rsuffix='_mirna')
    joined_new = joined_new[(joined_new['start'] <= joined_new['start_mirna']) &
                            (joined_new['stop'] >= joined_new['stop_mirna'])]
    localizations = pd.DataFrame(columns=['chrom', 'name', 'id', 'start', 'stop', 'orientation', 'based_on_coordinates',
                                          'arm'])
    for index, row in matures.iterrows():
        if row['-/+'] == '+':
            localizations = localizations.append(
                {'chrom': row['chr'],
                 'id': row['From'],
                 'name': row['Name'] + '_pre-seed',
                 'start': row['start'],
                 'start_pre': row['start_pre'],
                 'stop': row['start'],
                 'orientation': row['-/+'],
                 'type': 'pre-seed',
                 'based_on_coordinates': 'yes',
                 'arm': row['arm']
                 }, ignore_index=True)
            localizations = localizations.append(
                {'chrom': row['chr'],
                 'name': row['Name'] + '_seed',
                 'id': row['From'],
                 'start': row['start']+1,
                 'start_pre': row['start_pre'],
                 'stop': row['start']+7,
                 'orientation': row['-/+'],
                 'type': 'seed',
                 'based_on_coordinates': 'yes',
                 'arm': row['arm']
                 }, ignore_index=True)
            localizations = localizations.append(
                {'chrom': row['chr'],
                 'name': row['Name'] + '_post-seed',
                 'id': row['From'],
                 'start': row['start']+8,
                 'start_pre': row['start_pre'],
                 'stop': row['stop'],
                 'orientation': row['-/+'],
                 'type': 'post-seed',
                 'based_on_coordinates': 'yes',
                 'arm': row['arm']
                 }, ignore_index=True
            )
        else:
            localizations = localizations.append(
                {'chrom': row['chr'],
                 'name': row['Name'] + '_pre_seed',
                 'id': row['From'],
                 'start': row['stop'],
                 'start_pre': row['start_pre'],
                 'stop': row['stop'],
                 'orientation': row['-/+'],
                 'type': 'pre-seed',
                 'based_on_coordinates': 'yes',
                 'arm': row['arm']
                 }, ignore_index=True)
            localizations = localizations.append(
                {'chrom': row['chr'],
                 'name': row['Name'] + '_seed',
                 'id': row['From'],
                 'start': row['stop']-7,
                 'start_pre': row['start_pre'],
                 'stop': row['stop']-1,
                 'orientation': row['-/+'],
                 'type': 'seed',
                 'based_on_coordinates': 'yes',
                 'arm': row['arm']
                 }, ignore_index=True)
            localizations = localizations.append(
                {'chrom': row['chr'],
                 'name': row['Name'] + '_post-seed',
                 'id': row['From'],
                 'start': row['start'],
                 'start_pre': row['start_pre'],
                 'stop': row['stop']-8,
                 'orientation': row['-/+'],
                 'type': 'post-seed',
                 'based_on_coordinates': 'yes',
                 'arm': row['arm']
                 }, ignore_index=True
            )
    temp = joined_new.groupby(['Alias', 'start', 'stop'], as_index=False)[['stop_mirna']].count()
    temp.columns = ['Alias', 'start', 'stop', 'count']
    both_mirnas = temp[temp['count'] == 2]
    single_mirna = temp[temp['count'] == 1]
    joined_both = joined_new[joined_new['Alias'].isin(
        both_mirnas['Alias'].unique())]
    joined_single = joined_new[joined_new['Alias'].isin(single_mirna['Alias'].unique())]

    for index, row in joined_both.drop_duplicates(['Alias', 'start', 'stop'], keep='first').iterrows():
        data = joined_both[(joined_both['Alias'] == row['Alias']) &
                           (joined_both['start'] == row['start']) &
                           (joined_both['stop'] == row['stop'])].copy().reset_index()
        coordinates = []
        coordinates += list(data['start'].unique())
        coordinates += list(data['stop'].unique())
        coordinates += list(data['start_mirna'].unique())
        coordinates += list(data['stop_mirna'].unique())
        coordinates = sorted(coordinates)

        if data.loc[0, '-/+'] == '+':
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_loop',
                 'id': row['Alias'],
                 'start': coordinates[2]+1,
                 'start_pre': coordinates[0],
                 'stop': coordinates[3]-1,
                 'orientation': '+',
                 'type': 'loop',
                 'based_on_coordinates': 'yes',
                 'arm': 'loop'
                 }, ignore_index=True
            )

        elif data.loc[0, '-/+'] == '-':

            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_loop',
                 'id': row['Alias'],
                 'start': coordinates[2]+1,
                 'start_pre': coordinates[0],
                 'stop': coordinates[3]-1,
                 'orientation': '-',
                 'type': 'loop',
                 'based_on_coordinates': 'yes',
                 'arm': 'loop'
                 }, ignore_index=True
            )

        else:
            print('not all records has the same orientation')
            return 1

    for index, row in joined_single.iterrows():
        data = joined_single[(joined_single['Alias'] == row['Alias']) &
                             (joined_single['start'] == row['start']) &
                             (joined_single['stop'] == row['stop'])].reset_index()
        temp_coordinates = []
        coordinates = []
        coordinates += list(data['start'].unique())
        coordinates += list(data['stop'].unique())
        coordinates += list(data['start_mirna'].unique())
        coordinates += list(data['stop_mirna'].unique())
        coordinates = sorted(coordinates)

        if (abs(coordinates[0] - coordinates[1]) < abs(coordinates[2] - coordinates[3])
                and data.loc[0, '-/+'] == '+'):
            known = 5
        elif (abs(coordinates[0] - coordinates[1]) < abs(coordinates[2] - coordinates[3])
              and data.loc[0, '-/+'] == '-'):
            known = 3
        elif data.loc[0, '-/+'] == '+':
            known = 3
        else:
            known = 5
        structure = structures[structures['name'] == data.loc[0, 'Alias']]
        if data.loc[0, '-/+'] == '+':
            if known == 5:
                predicted_start = coordinates[1] - coordinates[0] + 1
                pre_mirna = str(structure['structure'].values[0][:int(predicted_start)]).count('(')
                pre_mirna = pre_mirna - str(structure['structure'].values[0][:int(predicted_start)]).count(')')
                mirna_len = coordinates[2] - coordinates[1]
                mirna_len_com = str(
                    structure['structure'].values[0][int(predicted_start):int(predicted_start+mirna_len+1)]).count('(')
                mirna_len_com = mirna_len_com - str(
                    structure['structure'].values[0][int(predicted_start):int(predicted_start+mirna_len+1)]).count(')')
                closing = 0
                start_second_strand = coordinates[-1] - coordinates[0] + 2
                while closing < pre_mirna:
                    start_second_strand = start_second_strand - 1
                    closing = str(structure['structure'].values[0][int(start_second_strand):]).count(')')
                    closing = closing - str(structure['structure'].values[0][int(start_second_strand):]).count('(')
                closing2 = 0
                start_second_strand2 = coordinates[-1] - coordinates[0] + 2
                while closing2 < pre_mirna+mirna_len_com:
                    start_second_strand2 = start_second_strand2 - 1
                    closing2 = str(structure['structure'].values[0][int(start_second_strand2):]).count(')')
                    closing2 = closing2 - str(structure['structure'].values[0][int(start_second_strand2):]).count('(')

                stop_second_strand = start_second_strand2
                start_second_strand = start_second_strand + 2 + coordinates[0]
                stop_second_strand = stop_second_strand + 3 + coordinates[0]
                temp_coordinates = [start_second_strand, stop_second_strand]
                temp_coordinates.sort()

            if known == 3:
                predicted_start = coordinates[2] - coordinates[0] + 1
                pre_mirna = str(structure['structure'].values[0][int(predicted_start):]).count(')')
                pre_mirna = pre_mirna - str(structure['structure'].values[0][int(predicted_start):]).count('(')
                mirna_len = coordinates[2] - coordinates[1]
                mirna_len_com = str(
                    structure['structure'].values[0][int(predicted_start-mirna_len):int(predicted_start)]).count(')')
                mirna_len_com = mirna_len_com - str(
                    structure['structure'].values[0][int(predicted_start-mirna_len):int(predicted_start)]).count('(')
                closing = 0
                start_second_strand = 0
                while closing < pre_mirna:
                    start_second_strand = start_second_strand + 1
                    closing = str(structure['structure'].values[0][:int(start_second_strand)]).count('(')
                    closing = closing - str(structure['structure'].values[0][:int(start_second_strand)]).count(')')
                closing2 = 0
                start_second_strand2 = 0
                while closing2 < pre_mirna+mirna_len_com:
                    start_second_strand2 = start_second_strand2 + 1
                    closing2 = str(structure['structure'].values[0][:int(start_second_strand2)]).count('(')
                    closing2 = closing2 - str(structure['structure'].values[0][:int(start_second_strand2)]).count(')')

                stop_second_strand = start_second_strand2
                start_second_strand = start_second_strand + 2 + coordinates[0]
                stop_second_strand = stop_second_strand + 2 + coordinates[0]
                temp_coordinates = [start_second_strand, stop_second_strand]
                temp_coordinates.sort()
        if data.loc[0, '-/+'] == '-':
            if known == 5:
                predicted_start = coordinates[3] - coordinates[2]
                pre_mirna = str(structure['structure'].values[0][:int(predicted_start)]).count('(')
                pre_mirna = pre_mirna - str(structure['structure'].values[0][:int(predicted_start)]).count(')')
                mirna_len = coordinates[2] - coordinates[1]
                mirna_len_com = str(
                    structure['structure'].values[0][int(predicted_start):int(predicted_start+mirna_len+1)]).count('(')
                mirna_len_com = mirna_len_com - str(
                    structure['structure'].values[0][int(predicted_start):int(predicted_start+mirna_len+1)]).count(')')
                closing = 0
                start_second_strand = coordinates[3] - coordinates[0] + 1
                while closing < pre_mirna:
                    start_second_strand = start_second_strand - 1
                    closing = str(structure['structure'].values[0][int(start_second_strand):]).count(')')
                    closing = closing - str(structure['structure'].values[0][int(start_second_strand):]).count('(')
                closing2 = 0
                start_second_strand2 = coordinates[3] - coordinates[0] + 1
                while closing2 < pre_mirna+mirna_len_com:
                    start_second_strand2 = start_second_strand2 - 1
                    closing2 = str(structure['structure'].values[0][int(start_second_strand2):]).count(')')
                    closing2 = closing2 - str(structure['structure'].values[0][int(start_second_strand2):]).count('(')

                stop_second_strand = start_second_strand2
                start_second_strand = coordinates[3] - start_second_strand - 1
                stop_second_strand = coordinates[3] - stop_second_strand - 1
                temp_coordinates = [start_second_strand, stop_second_strand]
                temp_coordinates.sort()
            if known == 3:
                predicted_start = coordinates[3] - coordinates[1]
                pre_mirna = str(structure['structure'].values[0][int(predicted_start):]).count(')')
                pre_mirna = pre_mirna - str(structure['structure'].values[0][int(predicted_start):]).count('(')
                mirna_len = coordinates[2] - coordinates[1]
                mirna_len_com = str(
                    structure['structure'].values[0][int(predicted_start-mirna_len):int(predicted_start)]).count(')')
                mirna_len_com = mirna_len_com - str(
                    structure['structure'].values[0][int(predicted_start - mirna_len):int(predicted_start)]).count('(')
                closing = 0
                start_second_strand = 0
                if pre_mirna == 0:
                    start_second_strand = 0
                else:
                    while closing <= pre_mirna:
                        start_second_strand = start_second_strand + 1
                        closing = str(structure['structure'].values[0][:int(start_second_strand)]).count('(')
                        closing = closing - str(structure['structure'].values[0][:int(start_second_strand)]).count(')')

                closing2 = 0
                start_second_strand2 = 0
                while closing2 < pre_mirna+mirna_len_com:
                    start_second_strand2 = start_second_strand2 + 1
                    closing2 = str(structure['structure'].values[0][:int(start_second_strand2)]).count('(')
                    closing2 = closing2 - str(structure['structure'].values[0][:int(start_second_strand2)]).count(')')
                stop_second_strand = start_second_strand2
                start_second_strand = coordinates[3] - start_second_strand
                stop_second_strand = coordinates[3] - stop_second_strand - 1
                temp_coordinates = [start_second_strand, stop_second_strand]
                temp_coordinates.sort()
        coordinates = []
        coordinates += list(data['start'].unique())
        coordinates += list(data['stop'].unique())
        coordinates += list(data['start_mirna'].unique())
        coordinates += list(data['stop_mirna'].unique())
        coordinates += temp_coordinates
        coordinates = sorted(coordinates)

        if known == 3 and data.loc[0, '-/+'] == '+':
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_pred-5_pre-seed',
                 'id': row['Alias'],
                 'start': coordinates[1],
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[1],
                 'orientation': '+',
                 'type': 'pre-seed',
                 'based_on_coordinates': 'no',
                 'arm': '5p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_pred-5_seed',
                 'id': row['Alias'],
                 'start': coordinates[1]+1,
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[1]+7,
                 'orientation': '+',
                 'type': 'seed',
                 'based_on_coordinates': 'no',
                 'arm': '5p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_pred-5_post-seed',
                 'id': row['Alias'],
                 'start': coordinates[1]+8,
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[2],
                 'orientation': '+',
                 'type': 'post-seed',
                 'based_on_coordinates': 'no',
                 'arm': '5p'
                 }, ignore_index=True
            )
        elif known == 3 and data.loc[0, '-/+'] == '-':
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_pred-5_pre-seed',
                 'id': row['Alias'],
                 'start': coordinates[4],
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[4],
                 'orientation': '-',
                 'type': 'pre-seed',
                 'based_on_coordinates': 'no',
                 'arm': '5p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_pred-5_seed',
                 'id': row['Alias'],
                 'start': coordinates[4]-7,
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[4]-1,
                 'orientation': '-',
                 'type': 'seed',
                 'based_on_coordinates': 'no',
                 'arm': '5p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_pred-5_post-seed',
                 'id': row['Alias'],
                 'start': coordinates[3],
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[4]-8,
                 'orientation': '-',
                 'type': 'post-seed',
                 'based_on_coordinates': 'no',
                 'arm': '5p'
                 }, ignore_index=True
            )
        elif known == 5 and data.loc[0, '-/+'] == '+':
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_pred-3_pre-seed',
                 'id': row['Alias'],
                 'start': coordinates[3],
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[3],
                 'orientation': '+',
                 'type': 'pre-seed',
                 'based_on_coordinates': 'no',
                 'arm': '3p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_pred-3_seed',
                 'id': row['Alias'],
                 'start': coordinates[3]+1,
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[3]+7,
                 'orientation': '+',
                 'type': 'seed',
                 'based_on_coordinates': 'no',
                 'arm': '3p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_pred-3_post-seed',
                 'id': row['Alias'],
                 'start': coordinates[3]+8,
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[4],
                 'orientation': '+',
                 'type': 'post-seed',
                 'based_on_coordinates': 'no',
                 'arm': '3p'
                 }, ignore_index=True
            )
        else:
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_pred-3_pre-seed',
                 'id': row['Alias'],
                 'start': coordinates[2],
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[2],
                 'orientation': '-',
                 'type': 'pre-seed',
                 'based_on_coordinates': 'no',
                 'arm': '3p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_pred-3-seed',
                 'id': row['Alias'],
                 'start': coordinates[2]-7,
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[2]-1,
                 'orientation': '-',
                 'type': 'seed',
                 'based_on_coordinates': 'no',
                 'arm': '3p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_pred-3_post-seed',
                 'id': row['Alias'],
                 'start': coordinates[1],
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[2]-8,
                 'orientation': '-',
                 'type': 'post-seed',
                 'based_on_coordinates': 'no',
                 'arm': '3p'
                 }, ignore_index=True
            )
        if data.loc[0, '-/+'] == '+':

            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_loop',
                 'id': row['Alias'],
                 'start': coordinates[2] + 1,
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[3] - 1,
                 'orientation': '+',
                 'type': 'loop',
                 'based_on_coordinates': 'no',
                 'arm': 'loop'
                 }, ignore_index=True
            )

        elif data.loc[0, '-/+'] == '-':

            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_loop',
                 'id': row['Alias'],
                 'start': coordinates[2] + 1,
                 'start_pre': data['start'].unique()[0],
                 'stop': coordinates[3] - 1,
                 'orientation': '-',
                 'type': 'loop',
                 'based_on_coordinates': 'no',
                 'arm': 'loop'
                 }, ignore_index=True
            )

        else:
            print('not all records has the same orientation')
            return 1

    localizations['name'] = localizations['name'].str.lower()
    localizations.sort_values(['chrom', 'start', 'stop'], inplace=True)

    df_part = localizations.copy()

    df_part_minus = df_part[df_part['orientation'] == '-'].copy()
    df_part_plus = df_part[df_part['orientation'] == '+'].copy()

    df_part_minus3 = df_part_minus[(df_part_minus['arm'] == '3p') &
                                   (df_part_minus['type'] == 'post-seed')].copy()

    df_part_minus3['name'] = df_part_minus3['name'].str.extract(r'([a-z\-0-9]*)_', expand=False)
    df_part_minus3['name'] = df_part_minus3['name'].apply(lambda x: x + '_flanking-3')
    df_part_minus3['arm'] = '3p'
    df_part_minus3['stop'] = df_part_minus3['start'] - 1
    df_part_minus3['start'] = df_part_minus3['start'] - 25
    df_part_minus3['type'] = 'flanking-3'

    df_part_minus5 = df_part_minus[(df_part_minus['arm'] == '5p') &
                                   (df_part_minus['type'] == 'pre-seed')].copy()

    df_part_minus5['name'] = df_part_minus5['name'].str.extract(r'([a-z\-0-9]*)_', expand=False)
    df_part_minus5['name'] = df_part_minus5['name'].apply(lambda x: x + '_flanking-5')
    df_part_minus5['stop'] = df_part_minus5['start'] + 25
    df_part_minus5['start'] = df_part_minus5['start'] + 1
    df_part_minus5['arm'] = '5p'
    df_part_minus5['type'] = 'flanking-5'

    df_part_plus3 = df_part_plus[(df_part_plus['arm'] == '3p') &
                                 (df_part_plus['type'] == 'post-seed')].copy()

    df_part_plus3['name'] = df_part_plus3['name'].str.extract(r'([a-z\-0-9]*)_', expand=False)
    df_part_plus3['name'] = df_part_plus3['name'].apply(lambda x: x + '_flanking-3')
    df_part_plus3['start'] = df_part_plus3['stop'] + 1
    df_part_plus3['stop'] = df_part_plus3['stop'] + 25
    df_part_plus3['arm'] = '3p'
    df_part_plus3['type'] = 'flanking-3'

    df_part_plus5 = df_part_plus[(df_part_plus['arm'] == '5p') &
                                 (df_part_plus['type'] == 'pre-seed')].copy()

    df_part_plus5['name'] = df_part_plus5['name'].str.extract(r'([a-z\-0-9]*)_', expand=False)
    df_part_plus5['name'] = df_part_plus5['name'].apply(lambda x: x + '_flanking-5')
    df_part_plus5['stop'] = df_part_plus5['start'] - 1
    df_part_plus5['start'] = df_part_plus5['start'] - 25
    df_part_plus5['arm'] = '5p'
    df_part_plus5['type'] = 'flanking-5'

    localizations_new = pd.concat([df_part, df_part_minus3, df_part_minus5, df_part_plus3,
                                   df_part_plus5])
    localizations_new.sort_values(['chrom', 'start', 'stop'], inplace=True)
    localizations_new.dropna(inplace=True)
    localizations_new.to_csv(output_folder + '/.localizations.csv', sep=',',
                             index=False)
    coordinates = localizations_new.copy()
    coordinates = coordinates.groupby(['id', 'chrom', 'start_pre'], as_index=False).agg({
        'start': min,
        'stop': max
    })
    coordinates = coordinates.join(pre_mirbase.set_index('Alias'), on='id', how='inner', rsuffix='_temp')
    coordinates = coordinates[['chrom', 'start', 'stop', 'Name']]
    coordinates.sort_values(['chrom', 'start', 'stop'], inplace=True)
    coordinates['start'] = coordinates['start'].astype(int) - 1

    coordinates.drop_duplicates(keep='first').to_csv(output_folder + '/coordinates.bed', sep='\t',
                                                     index=False, header=False)
