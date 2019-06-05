import numpy as np


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

    if np.isnan(row['QSS_alt_tum']):
        results_qss = True
    elif int(row['QSS_alt_tum']) / row['tumor_alt_count'] > 20:
        results_qss = True

    return results_count and results_ssc and results_bq and results_qss


def add_coordinates(row, coordinates):
    result = coordinates[(coordinates['chr'] == row['chrom']) &
                         (coordinates['start_ref'] < row['pos']) &
                         (coordinates['stop_ref'] > row['pos'])]['gene'].values[0]
    return result
