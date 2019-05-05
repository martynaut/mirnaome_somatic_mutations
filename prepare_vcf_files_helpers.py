import re
import gzip


def change_info(col, file_type):
    result = list(col.values)[0]
    try:
        if file_type == 'varscan2':
            new_col = col.str.extract(r':SSC=([0-9]{0,})', expand=True)
            new_col2 = col.str.extract(r':SPV=([0-9/.Ee/-]{0,})', expand=True)
            ssc_max = new_col.astype('float').max().values[0]
            spv_min = new_col2.astype('float').min().values[0]
            result = re.sub(r':SSC=([0-9.]*)', ':SSC='+str(ssc_max), result)
            result = re.sub(r':SPV=([0-9.Ee\-]*)', ':SPV='+str(spv_min), result)
        else:
            result = list(col.values)[0]
    except Exception as e:
        print(e)
        print('sth went wrong, look it up please')
        print(col)
    return result


def change_format(col, file_type, column_name):
    first = list(col)[0].split(':')
    try:
        if file_type == 'muse':
            if column_name in ['AD_normal', 'AD_tumor']:
                ad_values = col.str.split(',', expand=True)
                ad_values = ad_values.astype(float).sum(axis=0)
                ad_values = ','.join(map(lambda x: str(int(x)), ad_values.values))
                return ad_values
            elif column_name in ['BQ_normal', 'BQ_tumor']:
                bq_values = col.str.split(',', expand=True)
                bq_values = bq_values.astype(float).max(axis=0)
                bq_values = ','.join(map(lambda x: str(int(x)), bq_values.values))
                return bq_values

        elif file_type == 'mutect2':

            if column_name in ['AD_normal', 'AD_tumor']:
                ad_values = col.str.split(',', expand=True)
                if ad_values.values is [['.']]:
                    return '.'
                else:
                    ad_values = ad_values.astype(float).sum(axis=0)
                    ad_values = ','.join(map(lambda x: str(int(x)), ad_values.values))
                    return ad_values
            elif column_name in ['QSS_normal', 'QSS_tumor']:
                qss_values = col.str.split(',', expand=True)
                qss_values = qss_values.astype(float).sum(axis=0)
                qss_values = ','.join(map(lambda x: str(int(x)), qss_values.values))
                return qss_values

        elif file_type == 'somaticsniper':

            if column_name in ['BCOUNT_normal', 'BCOUNT_tumor']:
                bc_values = col.str.replace('.', '0').str.split(',', expand=True)
                bc_values = bc_values.astype(float).sum(axis=0)
                bc_values = ','.join(map(lambda x: str(int(x)), bc_values.values))
                return bc_values
            elif column_name in ['SSC_normal', 'SSC_tumor']:
                ssc_values = col.str.replace('.', '0').str.split(',', expand=True)
                ssc_values = ssc_values.astype(float).max(axis=0)
                ssc_values = ','.join(map(lambda x: str(int(x)), ssc_values.values))
                return ssc_values

        else:

            if column_name in ['AD_normal', 'AD_tumor']:
                ad_values = col.str.replace('.', '0').str.split(',', expand=True)
                ad_values = ad_values.astype(float).sum(axis=0)
                ad_values = ','.join(map(lambda x: str(int(x)), ad_values.values))
                return ad_values
            elif column_name in ['RD_normal', 'RD_tumor']:
                rd_values = col.str.split(',', expand=True)
                rd_values = rd_values.astype(float).sum(axis=0)
                rd_values = ','.join(map(lambda x: str(int(x)), rd_values.values))
                return rd_values

    except Exception as e:
        raise e

    return ':'.join(first)


def update_dict_with_file(filename, dict_of_files):
    dict_entry = {
        'indiv_name': '',
        'indiv_id': '',
        'sample_id_tumor_name': '',
        'sample_id_tumor_aliQ': '',
        'sample_id_normal_name': '',
        'sample_id_normal_aliQ': '',
        'type_of_file': ''
        }

    regex_indiv_name = r'##INDIVIDUAL=<NAME=([A-Za-z0-9\-]+),'
    regex_indiv_id = r'##INDIVIDUAL=<NAME=[A-Za-z0-9\-]+,ID=([A-Za-z0-9\-]+)>'
    regex_sample_tumor_name = r'##SAMPLE=<ID=TUMOR,NAME=([A-Za-z0-9\-]+),'
    regex_sample_tumor_aliq = r'##SAMPLE=<ID=TUMOR,NAME=[A-Za-z0-9\-]+,ALIQUOT_ID=([A-Za-z0-9\-]+),'
    regex_sample_normal_name = r'##SAMPLE=<ID=NORMAL,NAME=([A-Za-z0-9\-]+),'
    regex_sample_normal_aliq = r'##SAMPLE=<ID=NORMAL,NAME=[A-Za-z0-9\-]+,ALIQUOT_ID=([A-Za-z0-9\-]+),'
    regex_type_of_file = r'##gdcWorkflow=<ID=somatic_mutation_calling_workflow,Name=([A-Za-z0-9\-]+),'

    with gzip.open(filename) as f:
        print('checking:')
        print(filename)
        for line in f.readlines():
            line = line.decode('ascii')
            if line.startswith('##INDIVIDUAL='):
                name_search = re.search(regex_indiv_name, line)
                if name_search:
                    dict_entry['indiv_name'] = name_search.group(1)
                id_search = re.search(regex_indiv_id, line)
                if id_search:
                    dict_entry['indiv_id'] = id_search.group(1)

            elif line.startswith('##SAMPLE=<ID=NORMAL'):
                sample_normal_name_search = re.search(regex_sample_normal_name, line)
                if sample_normal_name_search:
                    dict_entry['sample_id_normal_name'] = sample_normal_name_search.group(1)
                sample_normal_aliq_search = re.search(regex_sample_normal_aliq, line)
                if sample_normal_aliq_search:
                    dict_entry['sample_id_normal_aliQ'] = sample_normal_aliq_search.group(1)

            elif line.startswith('##SAMPLE=<ID=TUMOR'):
                sample_tumor_name_search = re.search(regex_sample_tumor_name, line)
                if sample_tumor_name_search:
                    dict_entry['sample_id_tumor_name'] = sample_tumor_name_search.group(1)
                sample_tumor_aliq_search = re.search(regex_sample_tumor_aliq, line)
                if sample_tumor_aliq_search:
                    dict_entry['sample_id_tumor_aliQ'] = sample_tumor_aliq_search.group(1)

            elif line.startswith('##gdcWorkflow=<ID=somatic_mutation_calling_workflow'):
                file_type_search = re.search(regex_type_of_file, line)
                if file_type_search:
                    dict_entry['type_of_file'] = file_type_search.group(1)

    dict_of_files[filename] = dict_entry
    return dict_of_files
