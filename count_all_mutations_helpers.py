import os
from prepare_vcf_files_helpers import update_dict_with_file
import pandas as pd
import gzip
import numpy as np
from extract_results_for_mirnaome_helpers import retract_counts, retract_info
import threading
from queue import Queue


def file_merge_algorithm(file):
    temp_df = pd.read_csv(file)
    temp_df = temp_df[['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format',
                       'tumor', 'normal', 'indiv_name', 'indiv_id', 'sample_id_tumor_name',
                       'sample_id_tumor_aliQ', 'sample_id_normal_name',
                       'sample_id_normal_aliQ', 'norm_ref_count', 'norm_alt_count',
                       'tumor_ref_count', 'tumor_alt_count', 'BQ_ref_tum', 'BQ_alt_tum',
                       'BQ_ref_norm', 'BQ_alt_norm', 'QSS_ref_tum', 'QSS_alt_tum',
                       'QSS_ref_nor', 'QSS_alt_nor', 'SSC', 'SPV', 'eval', 'algorithm']]

    temp_df['control:mut/norm'] = temp_df['norm_alt_count'] / temp_df['norm_ref_count']
    temp_df['tumor:mut/norm'] = temp_df['tumor_alt_count'] / temp_df['tumor_ref_count']

    temp_df['ratio'] = temp_df['tumor:mut/norm'] / temp_df['control:mut/norm']

    temp_df.replace(np.inf, 0, inplace=True)

    temp_df.columns = temp_df.columns.str.lower()
    temp_df.drop(['control:mut/norm', 'tumor:mut/norm', 'ratio', 'eval',
                  'qual',
                  'filter', 'info', 'format', 'normal', 'tumor',
                  'indiv_id', 'sample_id_tumor_name',
                  'sample_id_tumor_aliq', 'sample_id_normal_name',
                  'sample_id_normal_aliq'
                  ], axis=1, inplace=True)
    if temp_df.shape[0] == 0:
        print('no mutations found')
        return 0
    temp_df['mutation_type'] = temp_df.apply(lambda x: type_of_mutation(x), axis=1)

    temp_df.fillna(-1, inplace=True)

    all_mutations = temp_df.groupby(['chrom', 'pos',
                                     'indiv_name', 'ref', 'alt',
                                     'mutation_type']).agg({'algorithm': concat_alg,
                                                            'norm_ref_count': sum,
                                                            'norm_alt_count': sum,
                                                            'tumor_ref_count': sum,
                                                            'tumor_alt_count': sum
                                                            }).reset_index()
    all_mutations['type_of_subst'] = all_mutations.apply(lambda x: subst_type(x), axis=1)

    all_mutations.to_csv(file.split('.')[0] + '_algorithms_merged.csv',
                         sep=',',
                         index=False)


def subst_type(row):
    if row['mutation_type'] == 'subst':
        if (((row['ref'] in ['A', 'G']) and (row['alt'] in ['A', 'G'])) or
                ((row['ref'] in ['C', 'T']) and (row['alt'] in ['C', 'T']))):
            return 'transition'
        else:
            return 'transversion'
    else:
        return 'n.a.'


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


def each_file_processing(filename, dict_with_files):

    print(filename)

    dataframe_records = pd.DataFrame()

    with gzip.open(filename) as f:
        columns = []
        for line in f.readlines():
            line = line.decode('ascii')
            if line[:2] == '##':
                pass
            elif line[:1] == '#':
                columns = line.replace('#', '').strip().lower().split('\t')
            else:
                new_record = pd.DataFrame([line.replace('\n',
                                                        '').replace(';',
                                                                    ':').replace('"',
                                                                                 '').split('\t')],
                                          columns=columns)
                new_record['indiv_name'] = dict_with_files[filename]['indiv_name']
                new_record['indiv_id'] = dict_with_files[filename]['indiv_id']
                new_record['sample_id_tumor_name'] = dict_with_files[filename]['sample_id_tumor_name']
                new_record['sample_id_tumor_aliQ'] = dict_with_files[filename]['sample_id_tumor_aliQ']
                new_record['sample_id_normal_name'] = dict_with_files[filename]['sample_id_normal_name']
                new_record['sample_id_normal_aliQ'] = dict_with_files[filename]['sample_id_normal_aliQ']
                (new_record['norm_ref_count'], new_record['norm_alt_count'],
                    new_record['tumor_ref_count'], new_record['tumor_alt_count'],
                    new_record['BQ_ref_tum'],
                    new_record['BQ_alt_tum'],
                    new_record['BQ_ref_norm'],
                    new_record['BQ_alt_norm'],
                    new_record['QSS_ref_tum'],
                    new_record['QSS_alt_tum'],
                    new_record['QSS_ref_nor'],
                    new_record['QSS_alt_nor'],
                    new_record['SSC']) = retract_counts(
                        new_record['normal'], new_record['tumor'], new_record['format'],
                        new_record['ref'], new_record['alt'])

                if np.isnan(float(new_record['SSC'].values[0])):
                    new_record['SSC'], new_record['SPV'] = retract_info(
                        new_record['info']
                    )
                else:
                    new_record['SPV'] = np.nan
                dataframe_records = pd.concat([dataframe_records, new_record])

    return dataframe_records


def handle_patient(input_tuple):
    output_folder, patient, rerun, temp_df, dict_with_files = input_tuple
    if os.path.isfile(output_folder + '/patients/results_count_all_{}.csv'.format(patient)) \
            and rerun:
        print('skipping files for: {}'.format(patient))
    else:
        print('starting analysis for {}'.format(patient))
        results_df = pd.DataFrame()
        count = 1
        lendf = temp_df.shape[0]

        for file in list(temp_df.index):
            print("{}: {}/{}".format(patient, count, lendf))
            dataframe = each_file_processing(file,
                                             dict_with_files
                                             )
            dataframe['algorithm'] = temp_df.loc[file, 'type_of_file']
            results_df = pd.concat([results_df, dataframe])
            count = count + 1
        results_df = results_df[results_df['filter'].str.contains('PASS')]
        results_df.to_csv(output_folder + '/patients/results_count_all_{}.csv'.format(patient),
                          sep=',',
                          index=False)
        print('finishing analysis for {}'.format(patient))


print_lock = threading.Lock()
url_queue = Queue()


def process_queue():
    while True:
        input_tuple = url_queue.get()
        handle_patient(input_tuple)
        url_queue.task_done()


def count_mutations(input_folder,  output_folder, rerun):

    if os.path.isfile(output_folder + '/files_summary_count_mutations.csv') and rerun:
        print('Skipping files check')
        df_dict_with_files = pd.read_csv(output_folder + '/files_summary_count_mutations.csv')
        df_dict_with_files.set_index('filename', inplace=True)
        dict_with_files = df_dict_with_files.to_dict(orient='index')

        if not os.path.exists(output_folder + '/patients'):
            os.makedirs(output_folder + '/patients')

    else:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        if not os.path.exists(output_folder + '/patients'):
            os.makedirs(output_folder + '/patients')

        files = [[x[0] + '/' + y for y in x[2]] for x in os.walk(input_folder)]
        files = [file for sublist in files for file in sublist]

        files_temp = [[x[0] + '/' + y for y in x[2]] for x in os.walk(output_folder + '/temp')]
        files_temp = [file for sublist in files_temp for file in sublist]

        files.extend(files_temp)

        with open(output_folder + '/do_not_use.txt', 'r') as file_dont:
            do_not_use = []
            for line in file_dont.readlines():
                do_not_use.append(line[:-1].split('/')[-1])
            gz_files = []
            for file in files:
                if file[-6:] == 'vcf.gz' and file.split('/')[-1] not in do_not_use:
                    gz_files.append(file)
                else:
                    pass

        # files summary
        dict_with_files = {}
        for gz_file in gz_files:
            dict_with_files = update_dict_with_file(gz_file, dict_with_files)

        df_dict_with_files = pd.DataFrame.from_dict(dict_with_files, orient='index')

        df_dict_with_files.index.name = 'filename'
        df_dict_with_files.to_csv(output_folder + '/files_summary_count_mutations.csv', sep=',')

    for i in range(10):
        t = threading.Thread(target=process_queue)
        t.daemon = True
        t.start()

    for patient in list(df_dict_with_files['indiv_name'].unique()):
        temp_df = df_dict_with_files.loc[df_dict_with_files['indiv_name'] == patient, :].copy()
        url_queue.put((output_folder, patient, rerun, temp_df, dict_with_files))

    url_queue.join()


def post_analyse(input_file, output_file):
    df = pd.read_csv(input_file)
    try:
        df['eval'] = df.apply(validation_function, axis=1)
    except ValueError:
        df['eval'] = True
    df.to_csv(output_file,
              sep=',',
              index=False)


def validation_function(row):
    results_ssc = False
    results_bq = False
    results_qss = False
    if row['tumor_alt_count'] >= 2 and row['norm_alt_count'] == 0:
        results_count = True
    elif row['tumor_alt_count'] >= 2 and \
            row['tumor_ref_count'] == 0 and \
            row['norm_alt_count'] / (row['norm_ref_count'] + row['norm_alt_count']) <= 0.2:
        results_count = True
    elif (row['tumor_alt_count'] >= 2 and
          (row['tumor_alt_count'] / (row['tumor_ref_count']
                                     + row['tumor_alt_count']))
          >= (row['norm_alt_count'] / (row['norm_ref_count']
                                       + row['norm_alt_count']))*5):
        results_count = True
    else:
        results_count = False

    if np.isnan(row['SSC']):
        results_ssc = True
    elif int(row['SSC']) > 30:
        results_ssc = True

    if np.isnan(row['BQ_alt_tum']):
        results_bq = True
    elif int(row['BQ_alt_tum']) > 20:
        results_bq = True
    try:
        if np.isnan(row['QSS_alt_tum']):
            results_qss = True
        elif int(row['QSS_alt_tum']) / row['tumor_alt_count'] > 20:
            results_qss = True
    except ZeroDivisionError:
        results_qss = False

    return results_count and results_ssc and results_bq and results_qss
