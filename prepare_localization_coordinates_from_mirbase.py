import pandas as pd
import click
import os
from Bio import SeqIO
import RNA


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
    df.to_csv(output_folder + '/hairpin_structure.csv', sep=',', index=False)


def prepare_hsa_files(hairpin_fa_file, mature_fa_file, chrom_coord, output_folder):
    with open(hairpin_fa_file) as hairpin, \
         open(mature_fa_file) as mature:
        names = ['chr', '.', 'type', 'start', 'stop', '.2', '-/+', '.3', 'desc']
        data_chr = pd.read_csv(chrom_coord, sep="\t", skiprows=13, names=names)
        data_chr['ID'], \
            data_chr['Alias'], \
            data_chr['Name'] = zip(*data_chr['desc'].apply(lambda x: x.split(';')[:3]))
        data_chr['ID'] = data_chr['ID'].str.extract(r'ID\=(.+)', expand=False)
        data_chr['Alias'] = data_chr['Alias'].str.extract(r'Alias\=(.+)', expand=False)
        data_chr['Name'] = data_chr['Name'].str.extract(r'Name\=(.+)', expand=False)
        data_chr.drop('desc', inplace=True, axis=1)
        data_chr.to_csv(output_folder + '/hsa_mirbase_coordinates.csv', index=False)
        data_hairpin = pd.DataFrame(columns=['id', 'seq'])
        for record in SeqIO.parse(hairpin, "fasta"):
            if 'hsa' in record.id:
                data_hairpin = data_hairpin.append(pd.DataFrame([[record.id, record.seq]],
                                                                columns=['id', 'seq']))
        data_hairpin.to_csv(output_folder + '/hsa_mirbase_hairpins.csv', index=False)
        data_mature = pd.DataFrame(columns=['id', 'seq'])
        for record in SeqIO.parse(mature, "fasta"):
            if 'hsa' in record.id:
                data_mature = data_mature.append(pd.DataFrame([[record.id, record.seq]],
                                                              columns=['id', 'seq']))
        data_mature.to_csv(output_folder + '/hsa_mirbase_matures.csv', index=False)


@click.command()
@click.argument('path_to_ViennaRNA')
@click.argument('hairpin_fa_file')
@click.argument('mature_fa_file')
@click.argument('chrom_coord')
@click.argument('output_folder')
def main(path_to_viennarna, hairpin_fa_file, mature_fa_file, chrom_coord, output_folder):
    click.echo('ViennaRNA used from {}'.format(path_to_viennarna))

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    prepare_hsa_files(hairpin_fa_file, mature_fa_file, chrom_coord, output_folder)
    prepare_structures(hairpin_fa_file, output_folder)

    return 0


if __name__ == "__main__":
    main()
