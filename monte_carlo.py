import subprocess
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description='Monte carlo simulation of error correction of hybrid surface-color codes')
parser.add_argument(
    'file', help='path to csv input file containing run parameters')
parser.add_argument(
    'start_row', type=int, help='first row of parameters to run')
parser.add_argument(
    'end_row', type=int, help='last row of parameters to run')

args = parser.parse_args()
path = args.file
start_row = args.start_row
end_row = args.end_row

df = pd.read_csv(path)
for i in range(start_row, end_row+1):
    row = df.iloc[[i]].squeeze()
    disorder_trials = row['disorder_avg']
    for j in range(disorder_trials):
        args = ['./build/hybridDecoder', row['L'], row['uP'], row['eP'], str(row['randomize_unencoding']).lower(), row['trials'], row['job_number']]
        args_string = [str(a) for a in args]
        print(args_string, 'disorder trial', j + 1, 'of', disorder_trials)
        subprocess.run(args_string)

