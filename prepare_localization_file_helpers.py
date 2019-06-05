import pandas as pd
import RNA
from Bio import SeqIO
from distinct_occure_helpers import find_arm


def prepare_structures(hairpin_fa_file, output_folder):
    df = pd.DataFrame(columns=['id', 'seq', 'structure'])
    fasta_sequences = SeqIO.parse(open(hairpin_fa_file), 'fasta')
    for fasta_sequence in fasta_sequences:
        if 'hsa' in fasta_sequence.id:
            df_temp = pd.DataFrame(data=[[fasta_sequence.id,
                                          str(fasta_sequence.seq),
                                          RNA.fold(str(fasta_sequence.seq))[0]]],
                                   columns=['id', 'seq', 'structure'])
            df = pd.concat([df, df_temp])
    df.reset_index(inplace=True)
    df.drop('index', axis=1, inplace=True)
    df.to_csv(output_folder + '/temp_reference/hairpin_structure.csv', sep=',', index=False)


def prepare_hsa_files(hsa_gff_mirbase_file, output_folder):
    names = ['chr', '.', 'type', 'start', 'stop', '.2', '-/+', '.3', 'desc']
    data_chr = pd.read_csv(hsa_gff_mirbase_file, sep="\t", skiprows=13, names=names)
    data_chr['ID'], \
        data_chr['Alias'], \
        data_chr['Name'] = zip(*data_chr['desc'].apply(lambda x: x.split(';')[:3]))
    data_chr['From'] = data_chr['desc'].apply(lambda x: x.split(';')[3] if len(x.split(';')) >= 4 else None)
    data_chr['From'] = data_chr['From'].str.extract(r'Derives_from\=(.+)', expand=False)
    data_chr['ID'] = data_chr['ID'].str.extract(r'ID\=(.+)', expand=False)
    data_chr['Alias'] = data_chr['Alias'].str.extract(r'Alias\=(.+)', expand=False)
    data_chr['Name'] = data_chr['Name'].str.extract(r'Name\=(.+)', expand=False)
    data_chr.drop('desc', inplace=True, axis=1)
    data_chr.to_csv(output_folder + '/temp_reference/hsa_mirbase_coordinates.csv', index=False)


def create_loc(coordinates_bed_file, output_folder):
    # localization based on mirbase
    loc_info = pd.read_csv(output_folder + '/temp_reference/hsa_mirbase_coordinates.csv', sep=',')
    primirnas = pd.read_csv(coordinates_bed_file,
                            names=['chr', 'start', 'stop', 'gene'], sep='\t')
    structures = pd.read_csv(output_folder + '/temp_reference/hairpin_structure.csv', sep=',')

    primirnas = primirnas[primirnas['gene'].str.contains('hsa')]
    pre_mirbase = loc_info[loc_info['type'] == 'miRNA_primary_transcript'].copy()

    matures = loc_info[loc_info['type'] == 'miRNA'].copy()
    pre_mirbase.drop(['.', '.2', '.3', 'Alias', 'type', 'From', 'chr'], inplace=True, axis=1)
    matures.drop(['.', '.2', '.3', 'Alias', 'type'], inplace=True, axis=1)
    pre_joined = primirnas.join(pre_mirbase.set_index('Name'), on='gene', how='inner', rsuffix='_mirbase')
    matures = matures.join(pre_mirbase.set_index('ID'), on='From', how='left', rsuffix='_pre')
    matures['arm'] = matures['Name'].str[-2:]
    matures.loc[~matures['arm'].isin(['3p',
                                      '5p']), 'arm'] = \
        matures[~matures['arm'].isin(['3p',
                                      '5p'])].apply(lambda x: find_arm(x), axis=1)

    matures['Name'] = matures['Name_pre'].map(str) + matures['Name'].str.extract(
        r'(\-[35]{1,1}p)', expand=False).fillna('')
    joined = pre_joined.join(matures.set_index('From'), on='ID', how='left', rsuffix='_mirna')
    localizations = pd.DataFrame(columns=['chrom', 'name', 'start', 'stop', 'orientation', 'based_on_coordinates',
                                          'arm'])

    for index, row in matures.iterrows():
        if row['-/+'] == '+':
            localizations = localizations.append(
                {'chrom': row['chr'],
                 'name': row['Name'] + '_pre-seed',
                 'start': row['start'],
                 'stop': row['start'],
                 'orientation': row['-/+'],
                 'type': 'pre-seed',
                 'based_on_coordinates': 'yes',
                 'arm': row['arm']
                 }, ignore_index=True)
            localizations = localizations.append(
                {'chrom': row['chr'],
                 'name': row['Name'] + '_seed',
                 'start': row['start']+1,
                 'stop': row['start']+7,
                 'orientation': row['-/+'],
                 'type': 'seed',
                 'based_on_coordinates': 'yes',
                 'arm': row['arm']
                 }, ignore_index=True)
            localizations = localizations.append(
                {'chrom': row['chr'],
                 'name': row['Name'] + '_post-seed',
                 'start': row['start']+8,
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
                 'start': row['stop'],
                 'stop': row['stop'],
                 'orientation': row['-/+'],
                 'type': 'pre-seed',
                 'based_on_coordinates': 'yes',
                 'arm': row['arm']
                 }, ignore_index=True)
            localizations = localizations.append(
                {'chrom': row['chr'],
                 'name': row['Name'] + '_seed',
                 'start': row['stop']-7,
                 'stop': row['stop']-1,
                 'orientation': row['-/+'],
                 'type': 'seed',
                 'based_on_coordinates': 'yes',
                 'arm': row['arm']
                 }, ignore_index=True)
            localizations = localizations.append(
                {'chrom': row['chr'],
                 'name': row['Name'] + '_post-seed',
                 'start': row['start'],
                 'stop': row['stop']-8,
                 'orientation': row['-/+'],
                 'type': 'post-seed',
                 'based_on_coordinates': 'yes',
                 'arm': row['arm']
                 }, ignore_index=True
            )

    temp = pd.DataFrame(joined.groupby('ID')['Name', 'stop_mirna'].count())
    temp.columns = ['how many mature mirnas in pre', 'count']
    both_mirnas = temp[temp['how many mature mirnas in pre'] == 2]
    single_mirna = temp[temp['how many mature mirnas in pre'] == 1]
    joined_both = joined[joined['ID'].isin(both_mirnas.index.values)]
    joined_single = joined[joined['ID'].isin(single_mirna.index.values)]

    for premirna in joined_both['ID'].unique():
        data = joined_both[joined_both['ID'] == premirna].reset_index()
        coordinates = []
        coordinates += list(data['start'].unique())
        coordinates += list(data['stop'].unique())
        coordinates += list(data['start_mirna'].unique())
        coordinates += list(data['stop_mirna'].unique())
        coordinates = sorted(coordinates)
        if data.loc[0, '-/+'] == '+':
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_5-flanking',
                 'start': coordinates[0],
                 'stop': coordinates[1]-1,
                 'orientation': '+',
                 'type': '5-flanking',
                 'based_on_coordinates': 'yes',
                 'arm': '5p'
                 }, ignore_index=True)
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_loop',
                 'start': coordinates[2]+1,
                 'stop': coordinates[3]-1,
                 'orientation': '+',
                 'type': 'loop',
                 'based_on_coordinates': 'yes',
                 'arm': 'loop'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_3-flanking',
                 'start': coordinates[4]+1,
                 'stop': coordinates[5],
                 'orientation': '+',
                 'type': '3-flanking',
                 'based_on_coordinates': 'yes',
                 'arm': '3p'
                 }, ignore_index=True)

        elif data.loc[0, '-/+'] == '-':

            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_5-flanking',
                 'start': coordinates[4]+1,
                 'stop': coordinates[5],
                 'orientation': '-',
                 'type': '5-flanking',
                 'based_on_coordinates': 'yes',
                 'arm': '5p'
                 }, ignore_index=True)
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_loop',
                 'start': coordinates[2] + 1,
                 'stop': coordinates[3] - 1,
                 'orientation': '-',
                 'type': 'loop',
                 'based_on_coordinates': 'yes',
                 'arm': 'loop'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_3-flanking',
                 'start': coordinates[0],
                 'stop': coordinates[1] - 1,
                 'orientation': '-',
                 'type': '3-flanking',
                 'based_on_coordinates': 'yes',
                 'arm': '3p'
                 }, ignore_index=True)

        else:
            print('not all records has the same orientation')
            return 1

    for premirna in joined_single['ID'].unique():
        temp_coordinates = []
        data = joined_single[joined_single['ID'] == premirna].reset_index()
        coordinates = []
        coordinates += list(data['start_mirbase'].unique())
        coordinates += list(data['stop_mirbase'].unique())
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
        structure = structures[structures['id'] == data.loc[0, 'gene']]

        if data.loc[0, '-/+'] == '+':
            if known == 5:
                predicted_start = coordinates[1] - coordinates[0] + 1
                pre_mirna = str(structure['structure'].str[:predicted_start].values).count('(')
                mirna_len = coordinates[2] - coordinates[1]
                mirna_len_com = str(
                    structure['structure'].str[predicted_start:predicted_start+mirna_len+1].values).count('(')
                closing = 0
                start_second_strand = coordinates[-1] - coordinates[0] + 2
                while closing < pre_mirna:
                    start_second_strand = start_second_strand - 1
                    closing = str(structure['structure'].str[start_second_strand:].values).count(')')
                closing2 = 0
                start_second_strand2 = coordinates[-1] - coordinates[0] + 2
                while closing2 < pre_mirna+mirna_len_com:
                    start_second_strand2 = start_second_strand2 - 1
                    closing2 = str(structure['structure'].str[start_second_strand2:].values).count(')')

                stop_second_strand = start_second_strand2
                start_second_strand = start_second_strand + 2 + coordinates[0]
                stop_second_strand = stop_second_strand + 3 + coordinates[0]
                temp_coordinates = [start_second_strand, stop_second_strand]
                temp_coordinates.sort()

            if known == 3:
                predicted_start = coordinates[2] - coordinates[0] + 1
                pre_mirna = str(structure['structure'].str[predicted_start:].values).count(')')

                mirna_len = coordinates[2] - coordinates[1]
                mirna_len_com = str(
                    structure['structure'].str[predicted_start-mirna_len:predicted_start].values).count(')')

                closing = 0
                start_second_strand = 0
                while closing < pre_mirna:
                    start_second_strand = start_second_strand + 1
                    closing = str(structure['structure'].str[:start_second_strand].values).count('(')
                closing2 = 0
                start_second_strand2 = 0
                while closing2 < pre_mirna+mirna_len_com:
                    start_second_strand2 = start_second_strand2 + 1
                    closing2 = str(structure['structure'].str[:start_second_strand2].values).count('(')

                stop_second_strand = start_second_strand2
                start_second_strand = start_second_strand + 2 + coordinates[0]
                stop_second_strand = stop_second_strand + 2 + coordinates[0]
                temp_coordinates = [start_second_strand, stop_second_strand]
                temp_coordinates.sort()
        if data.loc[0, '-/+'] == '-':
            if known == 5:
                predicted_start = coordinates[3] - coordinates[2]
                pre_mirna = str(structure['structure'].str[:predicted_start].values).count('(')
                mirna_len = coordinates[2] - coordinates[1]
                mirna_len_com = str(
                    structure['structure'].str[predicted_start:predicted_start+mirna_len+1].values).count('(')

                closing = 0
                start_second_strand = coordinates[3] - coordinates[0] + 1
                while closing < pre_mirna:
                    start_second_strand = start_second_strand - 1
                    closing = str(structure['structure'].str[start_second_strand:].values).count(')')
                closing2 = 0
                start_second_strand2 = coordinates[3] - coordinates[0] + 1
                while closing2 < pre_mirna+mirna_len_com:
                    start_second_strand2 = start_second_strand2 - 1
                    closing2 = str(structure['structure'].str[start_second_strand2:].values).count(')')

                stop_second_strand = start_second_strand2
                start_second_strand = coordinates[3] - start_second_strand - 1
                stop_second_strand = coordinates[3] - stop_second_strand - 1
                temp_coordinates = [start_second_strand, stop_second_strand]
                temp_coordinates.sort()
            if known == 3:
                predicted_start = coordinates[3] - coordinates[1]
                pre_mirna = str(structure['structure'].str[predicted_start:].values).count(')')
                mirna_len = coordinates[2] - coordinates[1]
                mirna_len_com = str(
                    structure['structure'].str[predicted_start-mirna_len:predicted_start].values).count(')')

                closing = 0
                start_second_strand = 0
                while closing < pre_mirna+1:
                    start_second_strand = start_second_strand + 1
                    closing = str(structure['structure'].str[:start_second_strand].values).count('(')
                if pre_mirna == 0:
                    start_second_strand = 0
                closing2 = 0
                start_second_strand2 = 0
                while closing2 < pre_mirna+mirna_len_com:
                    start_second_strand2 = start_second_strand2 + 1
                    closing2 = str(structure['structure'].str[:start_second_strand2].values).count('(')

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
                 'name': data.loc[0, 'Name'] + '-pred-5_pre-seed',
                 'start': coordinates[1],
                 'stop': coordinates[1],
                 'orientation': '+',
                 'type': 'pre-seed',
                 'based_on_coordinates': 'no',
                 'arm': '5p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '-pred-5_seed',
                 'start': coordinates[1]+1,
                 'stop': coordinates[1]+7,
                 'orientation': '+',
                 'type': 'seed',
                 'based_on_coordinates': 'no',
                 'arm': '5p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '-pred-5_post-seed',
                 'start': coordinates[1]+8,
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
                 'name': data.loc[0, 'Name'] + '-pred-5_pre-seed',
                 'start': coordinates[4],
                 'stop': coordinates[4],
                 'orientation': '-',
                 'type': 'pre-seed',
                 'based_on_coordinates': 'no',
                 'arm': '5p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '-pred-5_seed',
                 'start': coordinates[4]-7,
                 'stop': coordinates[4]-1,
                 'orientation': '-',
                 'type': 'seed',
                 'based_on_coordinates': 'no',
                 'arm': '5p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '-pred-5_post-seed',
                 'start': coordinates[3],
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
                 'name': data.loc[0, 'Name'] + '-pred-3_pre-seed',
                 'start': coordinates[3],
                 'stop': coordinates[3],
                 'orientation': '+',
                 'type': 'pre-seed',
                 'based_on_coordinates': 'no',
                 'arm': '3p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '-pred-3_seed',
                 'start': coordinates[3]+1,
                 'stop': coordinates[3]+7,
                 'orientation': '+',
                 'type': 'seed',
                 'based_on_coordinates': 'no',
                 'arm': '3p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '-pred-3_post-seed',
                 'start': coordinates[3]+8,
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
                 'name': data.loc[0, 'Name'] + '-pred-3_pre-seed',
                 'start': coordinates[2],
                 'stop': coordinates[2],
                 'orientation': '-',
                 'type': 'pre-seed',
                 'based_on_coordinates': 'no',
                 'arm': '3p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '-pred-3-seed',
                 'start': coordinates[2]-7,
                 'stop': coordinates[2]-1,
                 'orientation': '-',
                 'type': 'seed',
                 'based_on_coordinates': 'no',
                 'arm': '3p'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '-pred-3_post-seed',
                 'start': coordinates[1],
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
                 'name': data.loc[0, 'Name'] + '_5-flanking',
                 'start': coordinates[0],
                 'stop': coordinates[1] - 1,
                 'orientation': '+',
                 'type': '5-flanking',
                 'based_on_coordinates': 'no',
                 'arm': '5p'
                 }, ignore_index=True)
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_loop',
                 'start': coordinates[2] + 1,
                 'stop': coordinates[3] - 1,
                 'orientation': '+',
                 'type': 'loop',
                 'based_on_coordinates': 'no',
                 'arm': 'loop'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_3-flanking',
                 'start': coordinates[4] + 1,
                 'stop': coordinates[5],
                 'orientation': '+',
                 'type': '3-flanking',
                 'based_on_coordinates': 'no',
                 'arm': '3p'
                 }, ignore_index=True)

        elif data.loc[0, '-/+'] == '-':
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_5-flanking',
                 'start': coordinates[4] + 1,
                 'stop': coordinates[5],
                 'orientation': '-',
                 'type': '5-flanking',
                 'based_on_coordinates': 'no',
                 'arm': '5p'
                 }, ignore_index=True)
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_loop',
                 'start': coordinates[2] + 1,
                 'stop': coordinates[3] - 1,
                 'orientation': '-',
                 'type': 'loop',
                 'based_on_coordinates': 'no',
                 'arm': 'loop'
                 }, ignore_index=True
            )
            localizations = localizations.append(
                {'chrom': data.loc[0, 'chr'],
                 'name': data.loc[0, 'Name'] + '_3-flanking',
                 'start': coordinates[0],
                 'stop': coordinates[1] - 1,
                 'orientation': '-',
                 'type': '3-flanking',
                 'based_on_coordinates': 'no',
                 'arm': '3p'
                 }, ignore_index=True)
        else:
            print('not all records has the same orientation')
            return 1

    localizations['name'] = localizations['name'].str.lower()
    localizations.sort_values(['chrom', 'start', 'stop'], inplace=True)

    # change localization for +/- 25nt

    df_part = localizations[localizations['type'].isin(['seed', 'pre-seed',
                                                        'post-seed', 'loop'])].copy()

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

    localizations_new.to_csv(output_folder + '/localizations_test.csv', sep=',',
                             index=False)
