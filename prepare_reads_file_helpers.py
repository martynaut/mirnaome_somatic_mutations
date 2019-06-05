from lxml import html
import requests
import numpy as np


def parse_html(row):
    page = requests.get('http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc={}'.format(row['From']))
    tree = html.fromstring(page.content)

    text1 = tree.xpath('//a[@href="/cgi-bin/get_read.pl?acc={}"]/text()'.format(row['Alias']))

    try:
        return int(text1[0])
    except IndexError:
        return np.nan

