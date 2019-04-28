import gzip
import pandas as pd
import os
import numpy as np
from extract_results_for_mirnaome_helpers import retract_counts, retract_info
from prepare_vcf_files_helpers import update_dict_with_file


pd.options.mode.chained_assignment = None


def each_file_processing(filename, reference, dict_with_files):

    print(filename)

    dataframe_records = pd.DataFrame()

    try:
        f = gzip.open(filename)
        for line in f.readlines():
            try:
                line = line.decode('ascii')
            except AttributeError:
                continue
            if line[:1] == '#':
                pass
            else:
                position = line.split('\t')[:5]
                if reference[(reference['chr'] == position[0]) &
                             (reference['start_ref'] < int(position[1])) &
                             (reference['stop_ref'] > int(position[1]))].shape[0] > 0:

                    new_record = pd.DataFrame([line.replace('\n',
                                                            '').replace(';',
                                                                        ':').replace('"',
                                                                                     '').split('\t')],
                                              columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                                                       'INFO', 'FORMAT', 'NORMAL', 'TUMOR'])
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
                        new_record['NORMAL'], new_record['TUMOR'], new_record['FORMAT'],
                        new_record['REF'], new_record['ALT'])

                    if np.isnan(float(new_record['SSC'].values[0])):
                        new_record['SSC'], new_record['SPV'] = retract_info(
                            new_record['INFO']
                        )
                    else:
                        new_record['SPV'] = np.nan
                    dataframe_records = pd.concat([dataframe_records, new_record])
    except OSError:
        f = open(filename)
        for line in f.readlines():
            try:
                line = line.decode('ascii')
            except AttributeError:
                continue
            if line[:1] == '#':
                pass
            else:
                position = line.split('\t')[:5]
                if reference[(reference['chr'] == position[0]) &
                             (reference['start_ref'] < int(position[1])) &
                             (reference['stop_ref'] > int(position[1]))].shape[0] > 0:

                    new_record = pd.DataFrame([line.replace('\n',
                                                            '').replace(';',
                                                                        ':').replace('"',
                                                                                     '').split('\t')],
                                              columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                                                       'INFO', 'FORMAT', 'NORMAL', 'TUMOR'])
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
                        new_record['NORMAL'], new_record['TUMOR'], new_record['FORMAT'],
                        new_record['REF'], new_record['ALT'])

                    if np.isnan(float(new_record['SSC'].values[0])):
                        new_record['SSC'], new_record['SPV'] = retract_info(
                            new_record['INFO']
                        )
                    else:
                        new_record['SPV'] = np.nan
                    dataframe_records = pd.concat([dataframe_records, new_record])


        f.close()
    return dataframe_records


def all_files_processing(input_folder, output_folder, coordinates_file):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    coordinates = pd.read_table(coordinates_file,
                                names=['chr', 'start', 'stop', 'gene'])

    coordinates['start_ref'] = coordinates['start']
    coordinates['stop_ref'] = coordinates['stop']

    files = [[x[0] + '/' + y for y in x[2]] for x in os.walk(input_folder)]
    files = [file for sublist in files for file in sublist]
    print(files)
    files_temp = [[x[0] + '/' + y for y in x[2]] for x in os.walk(output_folder + '/temp')]
    files_temp = [file for sublist in files_temp for file in sublist]
    print(files_temp)

    files.extend(files_temp)

    with open(output_folder + '/do_not_use.txt', 'r') as file_dont:
        do_not_use = []
        for line in file_dont.readlines():
            do_not_use.append(line[:-1])
        gz_files = []
        for file in files:
            if (file.endswith('vcf.gz') or file.endswith('vcf')) and file not in do_not_use:
                gz_files.append(file)
            else:
                pass

    # files summary
    dict_with_files = {}
    for gz_file in gz_files:
        dict_with_files = update_dict_with_file(gz_file, dict_with_files)

    df_dict_with_files = pd.DataFrame.from_dict(dict_with_files, orient='index')

    df_dict_with_files.index.name = 'filename'
    df_dict_with_files.to_csv(output_folder + '/files_summary.csv', sep=',')

    df_dict_with_files_grouped = df_dict_with_files.reset_index().groupby(['indiv_name',
                                                                           'type_of_file']).agg('nunique')
    df_dict_with_files_grouped.to_csv(output_folder + '/files_summary_count_per_patient.csv', sep=',')

    df_file_type = df_dict_with_files.reset_index().groupby(['type_of_file'])[['filename']].agg('nunique')
    df_file_type.columns = ['count_filename']
    df_file_type.to_csv(output_folder + '/files_count_per_type.csv', sep=',')

    counter = 0
    all_files = df_dict_with_files.shape[0]
    for file_type in list(df_dict_with_files['type_of_file'].unique()):

        results_df = pd.DataFrame()

        for file in list(df_dict_with_files.loc[df_dict_with_files['type_of_file'] == file_type, :].index):
            counter += 1
            print(str(counter) + ' / ' + str(all_files) + '\n')
            dataframe = each_file_processing(file,
                                             coordinates,
                                             dict_with_files
                                             )
            results_df = pd.concat([results_df, dataframe])
        results_df = results_df[results_df['FILTER'].str.contains('PASS')]
        results_df.to_csv(output_folder + '/results_{}.csv'.format(file_type),
                          sep=',',
                          index=False)
