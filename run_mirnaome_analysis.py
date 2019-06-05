import click
from prepare_vcf_files import make_unique_files
from extract_results_for_mirnaome import all_files_processing
from merge_algorithms import filter_and_combine
from distinct_occure import dist_occur


@click.command()
@click.argument('input_folder')
@click.argument('output_folder')
@click.argument('coordinates_file')
@click.argument('localization_file')
@click.option('--from_step', '-s')
def main(input_folder,  output_folder, coordinates_file,
         localization_file,
         from_step=''):
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
        click.echo("Step 3: Merge all algorithms")
        filter_and_combine(output_folder)
    else:
        click.echo("Skipping step 3")
    if from_step <= 4:
        click.echo("Step 4: Make distinct and occure files")
        dist_occur(output_folder, localization_file)
    else:
        click.echo("Skipping step 4")

    click.echo("Analysis finished")


if __name__ == "__main__":
    main()
