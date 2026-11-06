import pandas as pd
import os

expr_dir = "/scratch/users/gsantamaria/test_large_files/metabognn_dev/results/labels/A549/expr/"
mc_dir = "/home/users/gsantamaria/projects/metabognn_dev/results/labels/A549/metchange/rank/"
e = pd.read_csv(os.path.join(expr_dir, "fit_info_nsig.csv"), index_col=0).n_sig
m = pd.read_csv(os.path.join(mc_dir, "fit_info_nsig.csv"), index_col=0).n_sig.reindex(e.index)
print(pd.crosstab(pd.cut(e, [-1, 0, 10, 100, 1000, 5000, 1e9]),
                  pd.cut(m, [-1, 0, 100, 1000, 4000, 1e9]), margins=True))
print("spearman:", round(e.corr(m, method="spearman"), 2))