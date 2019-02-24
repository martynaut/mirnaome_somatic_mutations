import pandas as pd
import click


@click.command()
@click.argument('confidence_file')
@click.argument('confidence_score_file')
@click.argument('aliases_file')
@click.argument('chromosome_build_file')
@click.argument('output_file')
def main(confidence_file, confidence_score_file,
         aliases_file, chromosome_build_file, output_file):
    df_confidence_file = pd.read_csv(confidence_file, sep='\t',
                                     names=['name', 'id', '1', '2', '3', '4', '5',
                                            '6', '7', '8', '9', '10',
                                            '11', '12', '13', '14', '15'])
    df_confidence_score_file = pd.read_csv(confidence_score_file, sep='\t',
                                           names=['id', 'score'])

    df_aliases_file = pd.read_csv(aliases_file, sep='\t',
                                  names=['id', 'aliases'])
    df_mirna_chromosome_build = pd.read_csv(chromosome_build_file, sep='\t',
                                            names=['id', '1', 'start', 'stop', 'Strand'])
    df_confidence_file = df_confidence_file[df_confidence_file['name'].str.contains('hsa')]
    df_aliases_file = df_aliases_file[~df_aliases_file['id'].str.contains('MIMAT')]
    df_aliases_file = df_aliases_file[df_aliases_file['aliases'].str.contains('hsa')]
    df_aliases_file['alias'] = df_aliases_file['aliases'].str.split(';').str[-2]
    df_aliases_file.drop('aliases', axis=1, inplace=True)

    df_aliases_file.drop_duplicates(inplace=True, keep='first')
    df_mirna_chromosome_build = df_mirna_chromosome_build[df_mirna_chromosome_build['1'].str.contains('chr')][['id',
                                                                                                               'start',
                                                                                                               'stop',
                                                                                                               'Strand']
                                                                                                              ]
    df_mirna_chromosome_build.drop_duplicates(keep='first', inplace=True)
    df_confidence_file = df_confidence_file.join(df_mirna_chromosome_build.set_index('id'), on='id')
    df_confidence_file = df_confidence_file.join(df_confidence_score_file.set_index('id'), on='id')[['name', 'score',
                                                                                                     'start',
                                                                                                     'stop',
                                                                                                     'Strand'
                                                                                                     ]]
    df_confidence_file = df_confidence_file.join(df_aliases_file.set_index('alias'), 'name')
    df_confidence_file['confidence'] = df_confidence_file['score'].apply(lambda x: 'High' if x > 0 else 'Low')
    df_confidence_file.to_excel(output_file)
    click.echo("Confidence file created")


if __name__ == "__main__":
    main()
