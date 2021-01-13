import csv
import subprocess

Ls = [9, 12, 18, 27]
# uPs = [0.2, 0.4, 0.6, 0.8]
uPs = [0, 1]
# ePs = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13]
ePs = { \
        0 : [0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09], \
        0.2 : [0.076, 0.077, 0.078, 0.079, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086], \
        0.4 : [0.079, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089], \
        0.6 : [0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.091, 0.092], \
        0.8 : [0.086, 0.087, 0.088, 0.089, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096], \
        1 : [0.096, 0.097, 0.098, 0.099, 0.1, 0.101, 0.102, 0.103, 0.104, 0.105, 0.106] \
}
trials = 1e4
job = 0
disorder = {0 : 1, 0.2 : 100, 0.4 : 100, 0.6 : 100, 0.8 : 100, 1 : 1}
git_hash = subprocess.check_output(['git', 'rev-parse',  'HEAD']).strip()

with open('19_10_20b.csv', 'w') as csv_file:
    writer = csv.writer(csv_file, delimiter=',')
    writer.writerow(['L', 'uP', 'eP', 'trials', 'disorder_avg', 'job_number', 'git hash (metadata)'])
    for uP in uPs:
        for L in Ls:
            for eP in ePs[uP]:
                writer.writerow([L, uP, eP, trials, disorder[uP], job, git_hash])
                job += 1
