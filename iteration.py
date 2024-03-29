import os
import pandas as pd
import numpy as np
import subprocess
import shutil

import sys

pd.set_option('display.max_rows', 500)


def read_times(pth):
    times = []

    for run in os.listdir(pth):

        if not run.startswith('pwave'):
            continue
        run_num = int(run.split('_')[2])

        # Run legion prof.
        prof_path = os.path.join(pth, run, 'run', 'legion_prof')
        subprocess.run(['legion_prof', os.path.join(pth, run, 'run', 'prof_0.gz'), '-o', prof_path],
                       stdout = subprocess.DEVNULL,
                       stderr = subprocess.DEVNULL)

        # Scrape the data we need.
        df = pd.read_csv(os.path.join(prof_path, "legion_prof_processor.tsv"), delimiter='\t')

        proc_file = df[df['text'] == 'CPU Proc 3'].iloc[0]['tsv']
        df = pd.read_csv(os.path.join(prof_path, proc_file), delimiter='\t')

        weija = df[df.apply(lambda r : 'AwaitMPI' in r['title'] or 'HandoffToMPI' in r['title'], axis=1)]
        weija = weija.sort_values('start')
        weija = weija.drop_duplicates('op_id', keep="last")
        weija = weija.reset_index()
        weija = weija[['title', 'start', 'end']]
        weija = weija.iloc[20:]
        weija = weija.reset_index()

        iters_to_skip = 3
        weija = weija.iloc[iters_to_skip * 11:]
        weija.reset_index()

        iters = len(weija) // 11
        iteration_times = []
        for it in range(iters):
            start = weija.iloc[0]['end']
            end = weija.iloc[10]['start']
            iteration_times.append((end - start) * 1e-6)
            weija = weija.iloc[11:]
            weija.reset_index()

        times.append((run_num, iteration_times))

        # Clean up after ourselves.
        shutil.rmtree(prof_path)

    times.sort(key=lambda x: x[0])
    for t, l in times:
        print(t, l, np.average(l))


read_times(sys.argv[1])
