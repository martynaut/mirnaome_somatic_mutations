import pandas as pd
import click
from prepare_reads_file_helpers import parse_html, prepare_hsa_files
import os


@click.command()
@click.argument('hsa_gff_mirbase_file')
@click.argument('output_file')
@click.argument('output_folder')
def main(hsa_gff_mirbase_file, output_file,
         output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder + '/temp_reference'):
        os.makedirs(output_folder + '/temp_reference')
    prepare_hsa_files(hsa_gff_mirbase_file, output_folder)
    loc_info = pd.read_csv(output_folder + '/temp_reference/hsa_mirbase_coordinates.csv', sep=',')
    pre_mirbase = loc_info[loc_info['type'] == 'miRNA_primary_transcript'].copy()
    matures = loc_info[loc_info['type'] == 'miRNA'].copy()

    pre_mirbase.drop(['.', '.2', '.3', 'Alias', 'type', 'From', 'chr'], inplace=True, axis=1)
    matures.drop(['.', '.2', '.3', 'type'], inplace=True, axis=1)

    joined_df = matures.join(pre_mirbase.set_index('ID'), on='From', how='left', rsuffix='_pre')

    joined_df['reads'] = joined_df.apply(lambda x: parse_html(x), axis=1)

    joined_df.to_csv(output_file, sep=',', index=False)

    click.echo("Reads file created")


if __name__ == "__main__":
    main()
