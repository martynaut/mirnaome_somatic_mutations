import sys
import click
import os
sys.path.append(sys.argv[1] + "/interfaces/Python3")
from prepare_localization_file_helpers import prepare_structures, prepare_hsa_files, create_loc


@click.command()
@click.argument('path_to_ViennaRNA')
@click.argument('hairpin_fa_file')
@click.argument('hsa_gff_mirbase_file')
@click.argument('coordinates_bed_file')
@click.argument('output_folder')
def main(path_to_viennarna, hairpin_fa_file, hsa_gff_mirbase_file,
         coordinates_bed_file, output_folder):
    click.echo('ViennaRNA used from {}'.format(path_to_viennarna))

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_folder + '/temp_reference'):
        os.makedirs(output_folder + '/temp_reference')
    prepare_structures(hairpin_fa_file, output_folder)
    prepare_hsa_files(hsa_gff_mirbase_file, output_folder)
    create_loc(coordinates_bed_file, output_folder)
    click.echo('Localization_file_created')
    return 0


if __name__ == "__main__":
    main()
