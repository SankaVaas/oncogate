import numpy as np, matplotlib.pyplot as plt
from config import *
from common import load_split, pair_matrices, gate_stats, tissue_rates

tr, te = load_split("train"), load_split("test")
s_cov, s_worst = tr.T.mean(0), tissue_rates(tr.H, tr.tissues).max().values
p_cov, p_worst = pair_matrices(tr)
iu = np.triu_indices(len(tr.genes), k=1)

caps = [0.001, 0.003, 0.01, 0.03, 0.05, 0.1]
rows = {"single": [], "pair": []}
for cap in caps:
    ok = np.where(s_worst <= cap, s_cov, -1)
    g = [tr.genes[ok.argmax()]] if ok.max() >= 0 else None
    rows["single"].append(gate_stats(te, g)[0] if g else 0)

    ok = np.where(p_worst[iu] <= cap, p_cov[iu], -1)
    k = ok.argmax()
    g = [tr.genes[iu[0][k]], tr.genes[iu[1][k]]] if ok.max() >= 0 else None
    rows["pair"].append(gate_stats(te, g)[0] if g else 0)
    print(f"cap {cap}: single {rows['single'][-1]:.3f} | pair {rows['pair'][-1]:.3f}  {g}")

for name, v in rows.items():
    plt.plot(caps, v, marker="o", label=name)
plt.xscale("log"); plt.xlabel("worst-tissue off-target cap (train)")
plt.ylabel("tumor coverage on held-out patients"); plt.legend()
plt.savefig(DATA / "pareto.png", dpi=150, bbox_inches="tight")
print("saved data/pareto.png")