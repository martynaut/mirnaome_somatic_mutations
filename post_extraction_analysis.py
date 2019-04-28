import pandas as pd
from post_extraction_analysis_helpers import validation_function, add_coordinates


def post_analyse(input_file, output_file):
    df = pd.read_csv(input_file)
    try:
        df['eval'] = df.apply(validation_function, axis=1)
    except ValueError:
        df['eval'] = True
    df.to_csv(output_file,
              sep=',',
              index=False)


def add_name(input_file, output_file, coordinates):
    df = pd.read_csv(input_file)
    coordinates = pd.read_table(coordinates,
                                names=['chr', 'start', 'stop', 'gene'])
    coordinates['start_ref'] = coordinates['start'].apply(lambda x: x - 0)
    coordinates['stop_ref'] = coordinates['stop'].apply(lambda x: x + 0)
    try:
        df['name'] = df.apply(lambda x: add_coordinates(x, coordinates), axis=1)
    except ValueError:
        df['name'] = ''
    df.to_csv(output_file,
              sep=',',
              index=False)
