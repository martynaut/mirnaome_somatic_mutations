import sys
import click
import os
sys.path.append(sys.argv[1] + "/interfaces/Python3")
from prepare_localization_file_helpers import prepare_structures, prepare_hsa_files, create_loc, \
    add_confidence, merge_all, create_reads, add_mirgenedb


@click.command()
@click.argument('path_to_ViennaRNA')
@click.argument('hairpin_fa_file')
@click.argument('hsa_gff_mirbase_file')
@click.argument('confidence_file')
@click.argument('confidence_score_file')
@click.argument('aliases_file')
@click.argument('chromosome_build_file')
@click.argument('mirgenedb_file')
@click.argument('output_folder')
def main(path_to_viennarna, hairpin_fa_file, hsa_gff_mirbase_file,
         confidence_file, confidence_score_file,
         aliases_file, chromosome_build_file,
         mirgenedb_file,
         output_folder):
    click.echo('ViennaRNA used from {}'.format(path_to_viennarna))

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    add_mirgenedb(mirgenedb_file, output_folder)
    add_confidence(confidence_file, confidence_score_file,
                   aliases_file, chromosome_build_file, output_folder)

    prepare_structures(hairpin_fa_file, output_folder)
    prepare_hsa_files(hsa_gff_mirbase_file, output_folder)
    create_reads(output_folder)
    create_loc(output_folder)

    merge_all(output_folder)
    click.echo('Localization_file_created')
    return 0


if __name__ == "__main__":
    main()
