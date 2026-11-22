# oncogate

Finding logic-gated (AND / NOT / OR) surface-marker panels that cover tumor cells
while keeping off-target risk low across a whole-body healthy atlas.

**Status:** computational candidates only. Nothing here says a marker is safe in a patient.

## Plan
1. `src/01_explore_census.py` - see what tumor/healthy data exists (no download of expression yet)
2. Build a surfaceome-filtered single-cell table (breast first)
3. Baseline: best single marker
4. Best AND-pair, then triples
5. Pareto curve: tumor coverage vs worst-tissue off-target rate, tested on held-out patients