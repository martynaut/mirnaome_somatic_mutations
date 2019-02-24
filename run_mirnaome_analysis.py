import click
from prepare_vcf_files import make_unique_files
from extract_results_for_mirnaome import all_files_processing
from post_extraction_analysis import post_analyse, add_name
from merge_algorithms import filter_and_combine
from distinct_occure import dist_occur
from chromosome_loc import chrom_localisation
from how_many_mirnas import check_covered_mirnas


@click.command()
@click.argument('input_folder')
@click.argument('output_folder')
@click.argument('coordinates_file')
@click.argument('confidence_file')
@click.argument('mirgenedb_file')
@click.argument('cancer_exon_file')
@click.argument('localization_file')
@click.argument('mirna_reads_file')
@click.option('--from_step', '-s')
@click.option('--chromosome_coverage', '-c')
def main(input_folder,  output_folder, coordinates_file,
         confidence_file, mirgenedb_file, cancer_exon_file,
         localization_file, mirna_reads_file,
         from_step='', chromosome_coverage=''):
    if not from_step:
        from_step = 0
    from_step = int(from_step)
    if from_step <= 1:
        click.echo("Step 1: Analysis started")
        make_unique_files(input_folder=input_folder, output_folder=output_folder)
        click.echo("Not unique files combined")
    else:
        click.echo("Skipping step 1")
    if from_step <= 2:
        click.echo("Step 2: Extract results for mirnaome")
        all_files_processing(input_folder, output_folder, coordinates_file)
    else:
        click.echo("Skipping step 2")
    if from_step <= 3:
        click.echo("Step 3: Post extraction data wrangling")
        post_analyse(output_folder + '/results_muse.csv',
                     output_folder + '/results_muse_eval.csv')
        post_analyse(output_folder + '/results_mutect2.csv',
                     output_folder + '/results_mutect2_eval.csv')
        post_analyse(output_folder + '/results_varscan2.csv',
                     output_folder + '/results_varscan2_eval.csv')
        post_analyse(output_folder + '/results_somaticsniper.csv',
                     output_folder + '/results_somaticsniper_eval.csv')
        add_name(output_folder + '/results_muse_eval.csv',
                 output_folder + '/results_muse_eval.csv', coordinates_file)
        add_name(output_folder + '/results_mutect2_eval.csv',
                 output_folder + '/results_mutect2_eval.csv', coordinates_file)
        add_name(output_folder + '/results_varscan2_eval.csv',
                 output_folder + '/results_varscan2_eval.csv', coordinates_file)
        add_name(output_folder + '/results_somaticsniper_eval.csv',
                 output_folder + '/results_somaticsniper_eval.csv', coordinates_file)
    else:
        click.echo("Skipping step 3")
    if from_step <= 4:
        click.echo("Step 4: Merge all algorithms")
        filter_and_combine(output_folder)
    else:
        click.echo("Skipping step 4")
    if from_step <= 5:
        click.echo("Step 5: Make distinct and occure files")
        dist_occur(output_folder, coordinates_file, confidence_file, mirgenedb_file, cancer_exon_file,
                   localization_file, mirna_reads_file)
    else:
        click.echo("Skipping step 5")
    if chromosome_coverage:
        chrom_localisation(output_folder, chromosome_coverage)
        check_covered_mirnas(output_folder, chromosome_coverage)
    click.echo("Analysis finished")


if __name__ == "__main__":
    main()
