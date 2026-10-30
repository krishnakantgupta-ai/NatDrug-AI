# Data sources, reuse and provenance

- **ChEMBL:** Measured IC50 rows targeting `CHEMBL203` (EGFR, UniProt `P00533`) used for the source training file `data/CHEMBL203_IC50.csv`. ChEMBL distributes its data under **CC BY-SA 3.0** with attribution and share-alike obligations: https://chembl.gitbook.io/chembl-interface-documentation/about#data-licensing . Cite the relevant ChEMBL database paper and exact release where known. The training export's specific ChEMBL release is **not frozen/verified**, so do not invent one.
- **IMPPAT:** `IMPPAT 3.0`, release `2026-09-30`, official download https://cb.imsc.res.in/imppat/download . Bulk data used for research-only local processing here. Public redistribution permission has **not been confirmed**, so the GitHub upload package omits its bulk exports and full individual-screening records.
- **BindingDB/NPASS:** These are identified as secondary sources in some IMPPAT rows. Some measurements match after unit conversion. Database labels are not proof of independent experimental publication.
- **External EGFR validation:** 13 candidate natural products were identified, but their original DOI/PMID and assay IDs were unavailable in the source export. **No verified independent natural-product validation metric is reported.**

This code package does not assert that a separate model R² on 13 provisional records represents external validation.

For a complete scientific release, obtain original per-measurement DOIs/PMIDs, assay IDs, comparable IC50 assay types, and a training-article nonoverlap audit.
