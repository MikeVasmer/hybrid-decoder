import csv
import subprocess

Ls = [9, 12, 18, 27]
uPs = [0.2, 0.4, 0.6, 0.8, 1]
ePs = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13]
trials = 1e3
job = 0
disorder = 100
git_hash = subprocess.check_output(['git', 'rev-parse',  'HEAD']).strip()
randomize_unencoding = True

with open('03_11_20.csv', 'w') as csv_file:
    writer = csv.writer(csv_file, delimiter=',')
    writer.writerow(['L', 'uP', 'eP', 'randomize_unencoding', 'trials', 'disorder_avg', 'job_number', 'git hash (metadata)'])
    for L in Ls:
        for uP in uPs:
            for eP in ePs:
                writer.writerow([L, uP, eP, randomize_unencoding, trials, disorder, job, git_hash])
                job += 1