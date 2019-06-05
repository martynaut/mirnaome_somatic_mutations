from lxml import html
import requests
import numpy as np
import pandas as pd


def parse_html(row):
    page = requests.get('http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc={}'.format(row['From']))
    tree = html.fromstring(page.content)

    text1 = tree.xpath('//a[@href="/cgi-bin/get_read.pl?acc={}"]/text()'.format(row['Alias']))

    try:
        return int(text1[0])
    except IndexError:
        return np.nan


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
