# Official IMPPAT datasets — download to this folder for local analysis

The public GitHub package does **not** republish bulk IMPPAT structure and activity tables or 18,030-row structure-level derived screening results. Their public redistribution rights were not verified.

For complete local reproducibility, obtain the **September 30, 2026 release of IMPPAT 3.0** from https://cb.imsc.res.in/imppat/download and place these two files here:

- `imppat_bioactivity.tsv`: the original bioactivity table (or extract its tabular source and retain a copy of the original).
- `imppat_chemical_normalized.csv`: normalized chemical export as used in this study (see `docs/DATA_SOURCES_AND_LICENSES.md`). Alternatively use the official chemical export in `code/prepare_imppat.py` after adjusting the notebook path and parser as appropriate.

Expected SHA256 checksums for *our original study inputs* (may differ for future database releases):

- Bioactivity TSV: `6a8aeaa8d7350c2513decf52733da50fbe1917e88cd6f187299a8b88b007e33e`
- Normalized chemical CSV: `95b7cb4855238c86051b44b9f65e0969962e1d6766f2d7f1773d7eef85f42a66`

The `research_v2/results/` public package includes only summary and ChEMBL-based results. To restore full cached natural-product results locally, see the **separate private local-data archive** supplied alongside the GitHub-ready ZIP.
