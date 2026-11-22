# Shared helpers.
from collections import namedtuple
import numpy as np
import pandas as pd
from config import *

Data = namedtuple("Data", "T H tissues genes")  # T: tumor cells x genes, H: healthy cells x genes (bool)


def load_surface_genes():
    """Gene symbols for surface proteins. Uses data/surfaceome_genes.txt if present."""
    p = DATA / "surfaceome_genes.txt"
    if not p.exists():
        sheets = pd.read_excel(SURFY_URL, sheet_name=None, header=1)
        for name, df in sheets.items():
            gcols = [c for c in df.columns if "gene" in str(c).lower()]
            print("sheet:", name, "| gene columns:", gcols)
            if "surfaceome" in name.lower() and gcols:
                genes = df[gcols[0]].dropna().astype(str).unique()
                p.write_text("\n".join(genes))
                break
        if not p.exists():
            raise RuntimeError("Auto-download failed. Put one gene symbol per line in "
                               "data/surfaceome_genes.txt (from the SURFY table) and rerun.")
    return set(p.read_text().split())


def load_split(split):
    df = pd.read_parquet(DATA / "cells.parquet")
    df = df[df["split"] == split]
    genes = [c for c in df.columns if c not in META_COLS]
    t, h = df[df.label == "tumor"], df[df.label == "healthy"]
    return Data(t[genes].values, h[genes].values, h["tissue_general"].values, genes)


META_COLS = ["soma_joinid", "dataset_id", "donor_id", "tissue_general",
             "cell_type", "disease", "label", "split"]


def tissue_rates(H, tissues):
    """Fraction of cells positive, per tissue. Rows = tissues, cols = genes (or 1 col for a gate)."""
    return pd.DataFrame(H).groupby(tissues).mean()


def gate_stats(d, gene_list):
    """AND-gate over gene_list. Returns (tumor coverage, worst-tissue off-target rate)."""
    idx = [d.genes.index(g) for g in gene_list]
    cov = d.T[:, idx].all(axis=1).mean()
    hit = d.H[:, idx].all(axis=1)
    worst = tissue_rates(hit.reshape(-1, 1), d.tissues).max().iloc[0]
    return cov, worst


def pair_matrices(d):
    """Coverage and worst-tissue off-target for ALL gene pairs at once (matrix multiply trick)."""
    T = d.T.astype(np.float32)
    cov = (T.T @ T) / len(T)
    worst = np.zeros_like(cov)
    for t in np.unique(d.tissues):
        Ht = d.H[d.tissues == t].astype(np.float32)
        worst = np.maximum(worst, (Ht.T @ Ht) / len(Ht))
    return cov, worst