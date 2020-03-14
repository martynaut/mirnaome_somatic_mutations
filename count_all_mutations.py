import click
from count_all_mutations_helpers import count_mutations, file_merge_algorithm
from count_all_mutations_helpers import post_analyse
import os
import pandas as pd
import threading
from queue import Queue


print_lock = threading.Lock()
url_queue = Queue()


def process_queue():
    while True:
        input_tuple = url_queue.get()
        print(input_tuple)
        file_merge_algorithm(input_tuple)
        url_queue.task_done()


@click.command()
@click.argument('input_folder')
@click.argument('output_folder')
@click.option('--from_step', '-s')
@click.option('--rerun', '-r')
def main(input_folder,  output_folder,
         from_step='', rerun=''):
    if not from_step:
        from_step = 0
    if not rerun:
        rerun = 0
    rerun = bool(rerun)
    from_step = int(from_step)
    print(from_step, rerun)
    click.echo("Analysis starting for:")
    click.echo(input_folder)
    if from_step <= 1:
        print("Initializing step 1")
        count_mutations(input_folder=input_folder, output_folder=output_folder,
                        rerun=rerun)
        print("Finished step 1")
    else:
        click.echo("Skipping step 1")

    if from_step <= 2:
        print("Initializing step 2")
        files_temp = [[x[0] + '/' + y for y in x[2] if ('_eval' not in y and 'results_count_all' in y)] for x
                      in os.walk(output_folder + '/patients')
                      ]
        files_temp = [file for sublist in files_temp for file in sublist]
        # print(files_temp)
        for file in files_temp:
            if not os.path.isfile(file.split('.')[0] + '_eval.csv'):

                post_analyse(file,
                             file.split('.')[0] + '_eval.csv')
            if not os.path.isfile(file.split('.')[0] + '_evaluated.csv'):
                file_df = pd.read_csv(file.split('.')[0] + '_eval.csv')

                file_df = file_df[file_df['eval']]
                file_df.to_csv(file.split('.')[0] + '_evaluated.csv')
        print("Finished step 2")
    else:
        click.echo("Skipping step 2")

    if from_step <= 3:
        print("Initializing step 3")
        files_temp = [[x[0] + '/' + y for y in x[2] if ('_evaluated.csv' in y and 'results_count_all' in y
                                                        )] for x
                      in os.walk(output_folder + '/patients')
                      ]
        files_temp = [file for sublist in files_temp for file in sublist]

        results = pd.DataFrame()
        for i in range(20):
            t = threading.Thread(target=process_queue)
            t.daemon = True
            t.start()

        for file in files_temp:
            # print(file)
            url_queue.put(file)

        url_queue.join()

        files_temp = [[x[0] + '/' + y for y in x[2] if ('_algorithms_merged.csv' in y and 'results_count_all' in y
                                                            )] for x
                          in os.walk(output_folder + '/patients')
                          ]
        files_temp = [file for sublist in files_temp for file in sublist]

        for file in files_temp:
            temp_df = pd.read_csv(file)
            user = temp_df['indiv_name'].unique()[0]
            count = temp_df.shape[0]
            new_record = pd.DataFrame([[user, count]], columns=['patient_id', 'mutation_count'])
            results = pd.concat([results, new_record])
        results.to_csv(output_folder + '/patient_mutation_count.csv', index=False)
        print("Finished step 3")
    else:
        click.echo("Skipping step 3")

    click.echo("Analysis finished")


if __name__ == "__main__":
    main()
