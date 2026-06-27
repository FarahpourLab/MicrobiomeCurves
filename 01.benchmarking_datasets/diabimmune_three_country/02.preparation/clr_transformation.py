"""CLR-transform a relative-abundance table (CLI helper).

Usage:
    python clr_transformation.py <relative_abundance.csv> <out_dir>

Zeros are replaced by a pseudo-count (half the smallest non-zero value in the
table) so the centred-log-ratio (CLR) is defined, then the transform is applied
per sample. Writes ``<out_dir>/clr_transformation_<input_name>_<pseudo_count>.csv``
and prints the pseudo-count used.
"""
import pandas as pd
import numpy as np
from skbio.stats.composition import clr
import os
import sys

file_path = sys.argv[1]
out_dir = sys.argv[2]

file_name = os.path.basename(file_path).split('.csv')[0]

data = pd.read_csv(file_path, index_col=0)
data.reset_index(inplace=True, drop=False)

sample_list = data.columns.tolist()
sample_list.pop(0)

tmp_min = 100
for sample in sample_list:
    for i in range(len(data)):
        if (tmp_min > data[sample][i]) and (data[sample][i] != 0.0):
            tmp_min = data[sample][i]
    pseudo_count = tmp_min / 2

data[sample_list] = data[sample_list].replace(0.0, 0.0 + pseudo_count)

data.set_index(data.columns[0], inplace=True)
data = data.T
clr_data = pd.DataFrame(clr(data))
clr_data.columns = data.columns
clr_data.index = data.index
clr_data = clr_data.T
print(pseudo_count)

os.makedirs(f"{out_dir}/", exist_ok=True)

clr_data.to_csv(f"{out_dir}/clr_transformation_{file_name}_{pseudo_count}.csv")
