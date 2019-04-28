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
            result = re.sub(r':SSC=([0-9.]*)', ':SSC='+str(ssc_max), list(col.values)[0])
            result = re.sub(r':SPV=([0-9.Ee\-]*)', ':SPV='+str(spv_min), result)
        else:
            result = list(col.values)[0]
    except Exception as e:
        print(e)
        print('sth went wrong, look it up please')
        print(col)
    return result


def change_format(col, file_type):
    first = list(col)[0].split(':')
    try:
        if file_type == 'muse':

            if len(list(col)[0].split(':')[1].split(',')) > 2:

                t_ad_value_1 = col.apply(lambda x: x.split(':')[1].split(',')[0])
                t_ad_value_2 = col.apply(lambda x: x.split(':')[1].split(',')[1])
                t_ad_value_3 = col.apply(lambda x: x.split(':')[1].split(',')[2])

                t_bq_1 = col.apply(lambda x: x.split(':')[2].split(',')[0])
                t_bq_2 = col.apply(lambda x: x.split(':')[2].split(',')[1])
                t_bq_3 = col.apply(lambda x: x.split(':')[2].split(',')[2])

                ad_value_1 = t_ad_value_1.astype('int').sum()
                ad_value_2 = t_ad_value_2.astype('int').sum()
                ad_value_3 = t_ad_value_3.astype('int').sum()

                bq_1 = t_bq_1.astype('float').max()
                bq_2 = t_bq_2.astype('float').max()
                bq_3 = t_bq_3.astype('float').max()

                first[1] = str(ad_value_1) + ',' + str(ad_value_2) + ',' + str(ad_value_3)
                first[2] = str(bq_1) + ',' + str(bq_2) + ',' + str(bq_3)

            else:

                t_ad_value_1 = col.apply(lambda x: x.split(':')[1].split(',')[0])
                t_ad_value_2 = col.apply(lambda x: x.split(':')[1].split(',')[1])

                t_bq_1 = col.apply(lambda x: x.split(':')[2].split(',')[0])
                t_bq_2 = col.apply(lambda x: x.split(':')[2].split(',')[1])

                ad_value_1 = t_ad_value_1.astype('int').sum()
                ad_value_2 = t_ad_value_2.astype('int').sum()

                bq_1 = t_bq_1.astype('float').max()
                bq_2 = t_bq_2.astype('float').max()

                first[1] = str(ad_value_1) + ',' + str(ad_value_2)
                first[2] = str(bq_1) + ',' + str(bq_2)

        elif file_type == 'mutect2':

            if len(list(col)[0].split(':')[1].split(',')) > 2:

                t_ad_value_1 = col.apply(lambda x: x.split(':')[1].split(',')[0])
                t_ad_value_2 = col.apply(lambda x: x.split(':')[1].split(',')[1])
                t_ad_value_3 = col.apply(lambda x: x.split(':')[1].split(',')[2])

                t_qss_1 = col.apply(lambda x: x.split(':')[-3].split(',')[0])
                t_qss_2 = col.apply(lambda x: x.split(':')[-3].split(',')[1])

                ad_value_1 = t_ad_value_1.astype('int').sum()
                ad_value_2 = t_ad_value_2.astype('int').sum()
                ad_value_3 = t_ad_value_3.astype('int').sum()

                qss_1 = t_qss_1.astype('int').sum()
                qss_2 = t_qss_2.astype('int').sum()

                first[1] = str(ad_value_1) + ',' + str(ad_value_2) + ',' + str(ad_value_3)
                first[-3] = str(qss_1) + ',' + str(qss_2)

            else:

                if list(col)[0].split(':')[1] == '.':
                    pass

                else:
                    t_ad_value_1 = col.apply(lambda x: x.split(':')[1].split(',')[0])
                    t_ad_value_2 = col.apply(lambda x: x.split(':')[1].split(',')[1])

                    ad_value_1 = t_ad_value_1.astype('int').sum()
                    ad_value_2 = t_ad_value_2.astype('int').sum()

                    first[1] = str(ad_value_1) + ',' + str(ad_value_2)

                t_qss_1 = col.apply(lambda x: x.split(':')[-3].split(',')[0])
                t_qss_2 = col.apply(lambda x: x.split(':')[-3].split(',')[1])

                qss_1 = t_qss_1.astype('int').sum()
                qss_2 = t_qss_2.astype('int').sum()

                first[-3] = str(qss_1) + ',' + str(qss_2)

        elif file_type == 'somaticsniper':

            t_ad_value_1 = col.apply(lambda x: x.split(':')[4].split(',')[0])
            t_ad_value_2 = col.apply(lambda x: x.split(':')[4].split(',')[1])
            t_ad_value_3 = col.apply(lambda x: x.split(':')[4].split(',')[2])
            t_ad_value_4 = col.apply(lambda x: x.split(':')[4].split(',')[3])

            t_ssc = col.apply(lambda x: x.split(':')[-1].replace('.', '0'))

            ad_value_1 = t_ad_value_1.astype('int').sum()
            ad_value_2 = t_ad_value_2.astype('int').sum()
            ad_value_3 = t_ad_value_3.astype('int').sum()
            ad_value_4 = t_ad_value_4.astype('int').sum()

            ssc = t_ssc.astype('int').max()

            first[4] = str(ad_value_1) + ',' + str(ad_value_2) + ',' + str(ad_value_3) + ',' + str(ad_value_4)
            first[-1] = str(ssc)

        else:

            t_rd = col.apply(lambda x: x.split(':')[3])
            t_ad = col.apply(lambda x: x.split(':')[4])

            rd = t_rd.astype('int').sum()
            ad = t_ad.astype('int').sum()

            first[3] = str(rd)
            first[4] = str(ad)
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
    regex_sample_tumor_name2 = r'##TUMOR="Sample=([A-Za-z0-9\-]+),'
    regex_sample_tumor_name = r'##SAMPLE=<ID=TUMOR,[A-Za-z]+=([A-Za-z0-9\-]+),'
    regex_sample_tumor_aliq = r'##SAMPLE=<ID=TUMOR,[A-Za-z]+=[A-Za-z0-9\-]+,ALIQUOT_ID=([A-Za-z0-9\-]+),'
    regex_sample_normal_name2 = r'##NORMAL="Sample=([A-Za-z0-9\-]+),'
    regex_sample_normal_name = r'##SAMPLE=<ID=NORMAL,[A-Za-z]+=([A-Za-z0-9\-]+),'
    regex_sample_normal_aliq = r'##SAMPLE=<ID=NORMAL,[A-Za-z]+=[A-Za-z0-9\-]+,ALIQUOT_ID=([A-Za-z0-9\-]+),'
    regex_type_of_file = r'##gdcWorkflow=<ID=somatic_mutation_calling_workflow,Name=([A-Za-z0-9\-]+),'

    print('checking:')
    print(filename)

    try:
        f = gzip.open(filename)

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

            elif line.startswith('##NORMAL'):
                sample_normal_name_search = re.search(regex_sample_normal_name2, line)
                if sample_normal_name_search:
                    dict_entry['sample_id_normal_name'] = sample_normal_name_search.group(1)
                dict_entry['indiv_name'] = '-'.join(sample_normal_name_search.group(1).split('-')[:3])

            elif line.startswith('##TUMOR'):
                sample_tumor_name_search = re.search(regex_sample_tumor_name2, line)
                if sample_tumor_name_search:
                    dict_entry['sample_id_tumor_name'] = sample_tumor_name_search.group(1)

            elif line.startswith('##MuSE_version'):
                dict_entry['type_of_file'] = 'muse'

            elif line.startswith('##GATKCommandLine.MuTect2'):
                dict_entry['type_of_file'] = 'mutect2'

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
        if dict_entry['indiv_name'] == '' and dict_entry['sample_id_normal_name'] != '':
            dict_entry['indiv_name'] = '-'.join(dict_entry['sample_id_normal_name'].split('-')[:3])
        if dict_entry['indiv_name'] == '' and dict_entry['sample_id_tumor_name'] != '':
            dict_entry['indiv_name'] = '-'.join(dict_entry['sample_id_tumor_name'].split('-')[:3])
    except OSError:
        f = open(filename)

        for line in f.readlines():
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

            elif line.startswith('##NORMAL'):
                sample_normal_name_search = re.search(regex_sample_normal_name2, line)
                if sample_normal_name_search:
                    dict_entry['sample_id_normal_name'] = sample_normal_name_search.group(1)
                    dict_entry['indiv_name'] = '-'.join(sample_normal_name_search.group(1).split('-')[:3])

            elif line.startswith('##TUMOR'):
                sample_tumor_name_search = re.search(regex_sample_tumor_name2, line)
                if sample_tumor_name_search:
                    dict_entry['sample_id_tumor_name'] = sample_tumor_name_search.group(1)

            elif line.startswith('##MuSE_version'):
                dict_entry['type_of_file'] = 'muse'

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
    f.close()

    dict_of_files[filename] = dict_entry
    return dict_of_files
