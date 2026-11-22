# Step 3: BASELINE. Best single marker = highest tumor coverage whose WORST healthy tissue
# stays under the cap. Everything later must beat this.
import pandas as pd
from config import *
from common import load_split, tissue_rates

tr, te = load_split("train"), load_split("test")

res = pd.DataFrame({
    "train_cov":   tr.T.mean(0),
    "train_worst": tissue_rates(tr.H, tr.tissues).max().values,
    "test_cov":    te.T.mean(0),
    "test_worst":  tissue_rates(te.H, te.tissues).max().reindex(tr.genes).values,
}, index=tr.genes)

ok = res[res.train_worst <= OFFTARGET_CAP].sort_values("train_cov", ascending=False)
print(f"Single markers with worst-tissue off-target <= {OFFTARGET_CAP}:")
print(ok.head(10).round(3))
print("\nMost common tumor markers overall (note how unsafe they are):")
print(res.sort_values("train_cov", ascending=False).head(5).round(3))