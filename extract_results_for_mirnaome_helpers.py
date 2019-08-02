import numpy as np
import re


def retract_counts(normal, tumor, data_format, ref, alt):

    ref = ref.values[0]
    alt = alt.values[0]

    # somaticsnipper
    if 'BCOUNT' in data_format.values[0]:
        try:
            bcount_pos = data_format.values[0].strip().split(':').index('BCOUNT')
            seq_dict = {
                'A': 0,
                'C': 1,
                'G': 2,
                'T': 3
            }
            norm_bcount = normal.values[0].split(':')[bcount_pos].split(',')
            tumor_bcount = tumor.values[0].split(':')[bcount_pos].split(',')
            norm_ref = int(norm_bcount[seq_dict.get(ref, seq_dict.get(ref[0], 0))])
            norm_alt = int(norm_bcount[seq_dict.get(alt, seq_dict.get(alt[0], 0))])
            tumor_ref = int(tumor_bcount[seq_dict.get(ref, seq_dict.get(ref[0], 0))])
            tumor_alt = int(tumor_bcount[seq_dict.get(alt, seq_dict.get(alt[0], 0))])
        except:
            norm_ref = 0
            norm_alt = 0
            tumor_ref = 0
            tumor_alt = 0

    # varscan2
    elif 'RD' in data_format.values[0]:
        try:
            rd_pos = data_format.values[0].strip().split(':').index('RD')
            ad_pos = data_format.values[0].strip().split(':').index('AD')
            norm_rd = normal.values[0].split(':')[rd_pos]
            norm_ad = normal.values[0].split(':')[ad_pos]
            norm_ref = int(norm_rd)
            norm_alt = int(norm_ad)
            tumor_rd = tumor.values[0].split(':')[rd_pos]
            tumor_ad = tumor.values[0].split(':')[ad_pos]
            tumor_ref = int(tumor_rd)
            tumor_alt = int(tumor_ad)
        except:
            norm_ref = 0
            norm_alt = 0
            tumor_ref = 0
            tumor_alt = 0

    # muse & mutect2
    elif 'AD' in data_format.values[0]:
        try:
            ad_pos = data_format.values[0].strip().split(':').index('AD')
            norm_ad = normal.values[0].split(':')[ad_pos].split(',')
            tumor_ad = tumor.values[0].split(':')[ad_pos].split(',')
            norm_ref = int(norm_ad[0])
            norm_alt = int(norm_ad[1])
            tumor_ref = int(tumor_ad[0])
            tumor_alt = int(tumor_ad[1])
        except:
            norm_ref = 0
            norm_alt = 0
            tumor_ref = 0
            tumor_alt = 0

    else:
        norm_ref = 0
        norm_alt = 0
        tumor_ref = 0
        tumor_alt = 0

    bq_ref_tum, bq_alt_tum, \
        bq_ref_norm, bq_alt_norm, qss_ref_tum, qss_alt_tum, qss_ref_norm, \
        qss_alt_not, ssc = [np.nan] * 9

    if 'SSC' in data_format.values[0]:
        ssc_pos = data_format.values[0].strip().split(':').index('SSC')
        ssc = float(tumor.values[0].split(':')[ssc_pos])

    if ':BQ:' in data_format.values[0]:
        bq_pos = data_format.values[0].strip().split(':').index('BQ')
        try:
            bq_ref_tum, bq_alt_tum = tumor.values[0].split(':')[bq_pos].split(',')[:2]
            bq_ref_norm, bq_alt_norm = normal.values[0].split(':')[bq_pos].split(',')[:2]
            bq_ref_tum, bq_alt_tum, bq_ref_norm, bq_alt_norm = \
                float(bq_ref_tum), float(bq_alt_tum), float(bq_ref_norm), float(bq_alt_norm)
        except ValueError:
            pass

    if 'QSS' in data_format.values[0]:
        qss_pos = data_format.values[0].strip().split(':').index('QSS')
        qss_ref_tum, qss_alt_tum = tumor.values[0].split(':')[qss_pos].split(',')[:2]
        qss_ref_norm, qss_alt_not = normal.values[0].split(':')[qss_pos].split(',')[:2]
        qss_ref_tum, qss_alt_tum, qss_ref_norm, qss_alt_not = \
            float(qss_ref_tum), float(qss_alt_tum), float(qss_ref_norm), float(qss_alt_not)

    return norm_ref, norm_alt, tumor_ref, tumor_alt, bq_ref_tum, bq_alt_tum, \
        bq_ref_norm, bq_alt_norm, qss_ref_tum, qss_alt_tum, qss_ref_norm, \
        qss_alt_not, ssc


def retract_info(info):
    ssc, spv = np.nan, np.nan

    ssc_search = re.search(r':SSC=([0-9.]*)', info.values[0], re.IGNORECASE)

    if ssc_search:
        ssc = float(ssc_search.group(1))

    spv_search = re.search(r':SPV=([0-9.Ee\-]*)', info.values[0], re.IGNORECASE)

    if spv_search:
        spv = float(spv_search.group(1))

    return ssc, spv
