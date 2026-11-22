import numpy as np
from config import *
from common import load_split, pair_matrices, gate_stats

tr, te = load_split("train"), load_split("test")
cov, worst = pair_matrices(tr)
np.fill_diagonal(worst, 1.0)   # diagonal = single genes; skip them here
cov_ok = np.where(worst <= OFFTARGET_CAP, cov, -1)

flat = np.argsort(cov_ok, axis=None)[::-1]
seen, shown = set(), 0
for f in flat:
    i, j = np.unravel_index(f, cov_ok.shape)
    if i >= j or cov_ok[i, j] < 0: continue
    genes = [tr.genes[i], tr.genes[j]]
    c_te, w_te = gate_stats(te, genes)
    print(f"{genes[0]} AND {genes[1]}: train cov {cov[i,j]:.3f}, train worst {worst[i,j]:.3f} "
          f"| TEST cov {c_te:.3f}, TEST worst {w_te:.3f}")
    shown += 1
    if shown == 10: break