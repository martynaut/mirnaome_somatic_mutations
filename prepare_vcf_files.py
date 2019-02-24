import gzip
import pandas as pd
import os
import shutil
from prepare_vcf_files_helpers import update_dict_with_file, change_format, change_info


pd.options.mode.chained_assignment = None


def make_unique_files(input_folder, output_folder):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if not os.path.exists(output_folder + '/temp'):
        os.makedirs(output_folder + '/temp')

    files = [[x[0] + '/' + y for y in x[2]] for x in os.walk(input_folder)]
    flat_files = [file for sublist in files for file in sublist]
    gz_files = [file for file in flat_files if file.endswith('vep.vcf.gz')]

    dict_with_files = {}
    for gz_file in gz_files:
        dict_with_files = update_dict_with_file(gz_file, dict_with_files)

    df_dict_with_files = pd.DataFrame.from_dict(dict_with_files, orient='index')
    df_dict_with_files.index.name = 'filename'

    df_dict_with_files_grouped = df_dict_with_files.reset_index().groupby(['indiv_name',
                                                                           'type_of_file']).agg('nunique')
    df_dict_with_files_grouped.to_csv(output_folder + '/files_summary_count_per_patient.csv', sep=',')

    df_not_unique_patients = df_dict_with_files_grouped.loc[df_dict_with_files_grouped['filename'] != 1, :]
    df_not_unique_patients.to_csv(output_folder + '/not_unique_patients.csv', sep=',')
    with open(output_folder+'/do_not_use.txt', 'w+') as do_not_use_file:
        for patient in list(df_not_unique_patients.unstack().index.unique()):
            this_patient = df_not_unique_patients.xs(patient, level=0)
            for file_type in list(this_patient.index.unique()):
                first = True
                with open(output_folder + '/temp/' + patient+'_'+file_type+'.vcf', 'wb') as combined:
                    temp_df = df_dict_with_files.loc[(df_dict_with_files['indiv_name'] == patient) &
                                                     (df_dict_with_files['type_of_file'] == file_type), :]

                    lines_df = pd.DataFrame()
                    for filename in list(temp_df.index.unique()):
                        print(filename)
                        do_not_use_file.write(filename+'\n')
                        with gzip.open(filename) as f:
                            for line in f.readlines():
                                dline = line.decode('ascii')
                                if dline.startswith('#') and first:
                                    combined.write(line)
                                elif dline.startswith('#'):
                                    pass
                                else:
                                    new_record = \
                                        pd.DataFrame([dline.replace('\n',
                                                                    '').replace(';',
                                                                                ':').replace('"',
                                                                                             '').split('\t')],
                                                     columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                                                              'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR'])
                                    lines_df = pd.concat([lines_df, new_record])
                        first = False
                    lines_df = lines_df[lines_df['FILTER'].str.contains('PASS')]
                    group_dict = {
                            'QUAL': 'first',
                            'FILTER': 'first',
                            'INFO': lambda x: change_info(x, file_type),
                            'FORMAT': 'first',
                            'NORMAL': lambda x: change_format(x, file_type),
                            'TUMOR': lambda x: change_format(x, file_type)
                        }
                    lines_df_grouped = lines_df.groupby(['CHROM', 'POS', 'ID',
                                                         'REF',
                                                         'ALT']).agg(group_dict).reset_index()
                    temp = []
                    for index, row in lines_df_grouped.iterrows():
                        temp.append(row.tolist())
                    for line in temp:
                        combined.write('\t'.join(line).encode('utf-8')+b'\n')
                with open(output_folder + '/temp/' + patient+'_'+file_type+'.vcf', 'rb') as combined, \
                        gzip.open(output_folder + '/temp/' + patient+'_'+file_type+'.vcf.gz', 'wb') as f_out:
                    shutil.copyfileobj(combined, f_out)
