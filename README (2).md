# NatDrug-AI: Research Studio

**An interactive, target-specific molecular screening platform for natural-product drug discovery**

**Version:** 1.0 (research prototype)  
**Developer:** Krishna Kant Gupta  
**Affiliation:** Rajiv Gandhi Institute of Information Technology and Biotechnology, Bharati Vidyapeeth (Deemed to be University), Pune, Maharashtra, India  
**Status:** Research software under development; not clinically validated

## Overview

NatDrug-AI combines molecular cheminformatics, experimental bioactivity data, machine learning, and an interactive web interface to explore natural compounds as potential leads against selected molecular targets. The Streamlit application uses RDKit to calculate physicochemical descriptors and molecular fingerprints, retrieves target and IC50 records from the ChEMBL REST API, and trains target-specific regression models using either a Random Forest or a PyTorch multilayer perceptron (MLP). Users can upload their own molecular libraries and export screening results.

**Important:** Predicted bioactivity is not experimental confirmation, clinical efficacy, safety, or drug approval. The software does not establish that any particular natural compound is an approved medicine.

## Features

- **Target discovery:** Search ChEMBL targets and inspect identifiers, names, organisms, and target types.
- **Approved-drug mechanism catalog:** Retrieve ChEMBL mechanism-of-action associations for approved drugs as evidence-linked target candidates. This is **not** an exhaustive list of all clinically validated targets, and associations must be checked individually.
- **Experimental bioactivity:** Retrieve exact IC50 measurements (nM) for a chosen ChEMBL target or upload your own measured IC50 CSV.
- **Molecular properties:** Calculate molecular weight, calculated LogP, H-bond donors and acceptors, topological polar surface area, rotatable bonds, QED, Lipinski rule violations, and Veber screening status.
- **Model training:** Train a Random Forest regressor or fingerprint-based PyTorch MLP to predict pIC50.
- **Cross-validation:** Use Murcko scaffold-grouped folds when enough distinct scaffolds are available; otherwise use random K-fold and clearly label the fallback.
- **Compound explorer:** Visualize a 2D molecular structure, inspect descriptors, and generate a target-specific prediction if a trained model is selected.
- **Batch screening:** Upload CSV or SDF molecular libraries, inspect predictions, and export results to CSV.
- **Model reports:** View saved model metadata, cross-validation metrics, and out-of-fold prediction plots.

## How it works

```text
ChEMBL IC50 records OR user-supplied measured IC50 CSV
                         |
               Structure standardization
                         |
            RDKit Morgan fingerprints
                  (radius 2, 2048 bits)
                         |
          Random Forest OR PyTorch MLP
                         |
           Cross-validation and refit
                         |
              Saved target model
                         |
         New natural compounds (CSV/SDF)
                         |
    Predicted pIC50 + descriptors + similarity
```

For experimental IC50 values in nanomolar units:

`pIC50 = 9 - log10(IC50_nM)`

The corresponding model-derived IC50 estimate is `IC50_nM = 10^(9 - predicted_pIC50)`.

**Drug-likeness is descriptor-based, not learned from approved-drug labels in this version.** QED, Lipinski and Veber metrics are screening aids; natural-product-derived drugs can fall outside these heuristic ranges.

## Requirements

- macOS or Linux (the project may also run on Windows with an appropriate Python environment)
- **Python 3.12 recommended**
- Internet access for live ChEMBL target and bioactivity retrieval
- Python packages from `requirements.txt`: Streamlit, RDKit, pandas, NumPy, scikit-learn, requests, joblib, Plotly and PyTorch

## Installation on macOS

### 1. Install Python 3.12 (if needed)

```bash
brew install python@3.12
```

Check the version:

```bash
$(brew --prefix python@3.12)/bin/python3.12 --version
```

### 2. Open the project folder

```bash
cd /Users/krishnagupta/NatDrug_AI
```

Replace this path with your actual project location if different.

### 3. Create a project-local virtual environment

```bash
$(brew --prefix python@3.12)/bin/python3.12 -m venv .venv
source .venv/bin/activate
```

### 4. Install dependencies

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### 5. Start the application

```bash
python -m streamlit run app.py
```

Open **http://localhost:8501** in your browser. Stop the server with **Control+C** in Terminal.

**If Anaconda interferes**, bypass shell activation and run:

```bash
./.venv/bin/python -m streamlit run app.py
```

## Quick start

### A. Find a target

1. Open **Clinical targets**.
2. Choose a therapeutic-area preset or type a target name (for example, EGFR).
3. Click **Search targets**.
4. Verify the organism, target type and **ChEMBL target ID** before proceeding.
5. Optionally build and download the approved-drug mechanism catalog. Mechanism association is not by itself proof of the target's suitability for small-molecule IC50 regression.

### B. Train a target-specific model

1. Open **Train AI**.
2. Choose **Live ChEMBL IC50** and enter the verified ChEMBL target ID; alternatively select **Upload measured IC50 CSV**.
3. Retrieve or upload measured bioactivity data.
4. Review the number of usable compounds and unique Murcko scaffolds.
5. Select **Random Forest** or **PyTorch MLP (deep learning)**.
6. Click **Train, cross-validate and save**.
7. Inspect MAE, RMSE, R² and the actual-versus-predicted pIC50 plot.
8. Download the curated dataset and out-of-fold predictions when available.

The application requires **at least 20 unique compounds** for exploratory training. Datasets with fewer than 50 unique compounds are flagged as exploratory. This threshold is a software safeguard, **not** a criterion for publication readiness.

### C. Analyze one natural compound

1. Select a saved model from the sidebar if you want a bioactivity prediction.
2. Open **Compound explorer**.
3. Paste a valid SMILES string and click **Analyze molecule**.
4. Review the molecular drawing, descriptors and optional target-specific pIC50 prediction.

### D. Upload your molecular library

1. Open **Batch screening**.
2. Upload a `.csv` file containing a `SMILES` column, or an `.sdf` file containing molecular structures.
3. Click **Screen uploaded compounds**.
4. Review valid and invalid structures and download `NatDrug_AI_predictions.csv`.

Example screening CSV:

```csv
Name,SMILES
Caffeine,Cn1c(=O)c2c(ncn2C)n(C)c1=O
Resveratrol,Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1
Aspirin,CC(=O)Oc1ccccc1C(=O)O
```

**Note:** A CSV containing only compound names and SMILES can be screened, but **cannot train a target-specific activity model** without experimentally measured IC50 values.

## Training dataset format

To train using your own measured bioactivity data, provide at least:

```csv
SMILES,IC50_nM
CCO,1500
CCN,420
```

These are **format examples only**, not genuine experimental measurements and not suitable as scientific training data. Replace them with verified experimental values and retain source, assay, organism and target annotations in your research records.

The current app canonicalizes valid SMILES, discards nonpositive or nonnumeric IC50 values, converts measurements to pIC50, and aggregates repeated structures using the median pIC50. Because assay conditions may differ, review and harmonize the source data before drawing scientific conclusions.

## Model architectures and validation

| Component | Implementation |
|---|---|
| Molecular representation | RDKit Morgan fingerprints, radius 2, 2,048 bits |
| Classical ML | RandomForestRegressor |
| Deep learning | PyTorch multilayer perceptron on fingerprints |
| Regression target | Experimental pIC50 derived from IC50 (nM) |
| Validation | Scaffold GroupKFold if at least 3 distinct scaffolds; random KFold otherwise |
| Metrics | MAE, RMSE and R² on cross-validated predictions |
| Structural similarity | Maximum Morgan-fingerprint Tanimoto similarity to training compounds |

The MLP is **not** a graph neural network, protein-language model or structure-based docking model. The maximum Tanimoto score is a similarity indicator, **not** a calibrated confidence interval or validated applicability-domain boundary.

## Understanding the results

- **QED:** An empirical quantitative estimate of drug-likeness (0–1), not a probability of approval.
- **Lipinski violations:** Count of common oral drug-likeness heuristic violations.
- **Veber pass:** Rule-based flexibility and polar surface area screen.
- **Predicted pIC50:** Estimated target-specific inhibitory potency learned from measured training compounds.
- **Estimated IC50 (nM):** Back-transformation of predicted pIC50; not a laboratory measurement.
- **Maximum Tanimoto:** Highest fingerprint similarity to the training set; low similarity may signal extrapolation.
- **MAE/RMSE/R²:** Cross-validation diagnostics, not external validation or proof of generalizability.

## Project files

```text
NatDrug_AI/
├── app.py                  # Streamlit interface and modeling workflow
├── requirements.txt        # Python dependencies
├── example_compounds.csv   # Example molecules for upload testing
├── README.md               # Documentation
├── data/                   # Downloaded IC50 datasets (created by app)
└── models/                 # Saved model artifacts and metadata (created by app)
```

## Troubleshooting

**`No module named streamlit`** — The wrong Python interpreter may be active. Check `which python` and install requirements with `./.venv/bin/python -m pip install -r requirements.txt`.

**`python3.11: command not found`** — Use the installed Homebrew Python 3.12 executable shown in the installation instructions.

**`source: no such file or directory: .venv/bin/activate`** — Create `.venv` in the current project folder first.

**ChEMBL search or download fails** — Check internet access and try again. Live API availability, rate limits, query filters and schema changes can affect retrieval.

**Training fails due to too few compounds** — Confirm at least 20 distinct valid structures with positive measured IC50 values. For meaningful benchmarking, substantially larger, well-curated datasets are preferable.

**Model file load error** — Load only models created by you or another trusted source. Serialized model files may be unsafe if obtained from untrusted parties.

## Research limitations and publication checklist

NatDrug-AI v1.0 is a **prototype**. Before making performance or novelty claims in a journal manuscript, the following should be completed and reported:

1. Curate target-specific assays, including species, target identity, assay type, confidence and measurement provenance.
2. Remove duplicates and prevent leakage between training and test structures; assess scaffold-disjoint external holdouts.
3. Compare Random Forest and MLP with appropriately tuned baseline models and report uncertainty across repeated splits.
4. Evaluate a separately curated natural-product test set and analyze chemical-space coverage.
5. Report dataset and ChEMBL release identifiers, software/package versions, random seeds and training settings.
6. Include failure cases, applicability-domain analyses, calibration where applicable, and reproducible tables/figures.
7. Confirm prioritized candidates using appropriate experimental assays before making biological or clinical claims.
8. Provide a public, versioned code repository, documented license and reproducibility instructions.

The software has not been demonstrated to outperform established tools, and no benchmark scores or experimental outcomes should be inferred from this README.

## Data sources and dependencies

- **ChEMBL:** https://www.ebi.ac.uk/chembl/ — target metadata, drug mechanism records and measured bioactivity.
- **RDKit:** https://www.rdkit.org/ — molecular parsing, descriptors, fingerprints and visualization.
- **Streamlit:** https://streamlit.io/ — interactive web application.
- **scikit-learn:** https://scikit-learn.org/ — Random Forest and validation metrics.
- **PyTorch:** https://pytorch.org/ — neural network implementation.

Consult the original database and software documentation for licensing, citation requirements and appropriate attribution.

## Citation and availability

**Suggested project reference (placeholder; not a published paper):**

Gupta, K. K. *NatDrug-AI: An Interactive Machine-Learning and Deep-Learning Research Platform for Target-Specific Natural-Compound Screening.* Software project, version 1.0 (unpublished).

**Source code:** Public repository URL to be added when released.  
**Manuscript:** Under preparation; no DOI assigned.  
**License:** Not yet specified. Do not assume an open-source license until one is explicitly added.

## Contact

**Krishna Kant Gupta**  
Rajiv Gandhi Institute of Information Technology and Biotechnology  
Bharati Vidyapeeth (Deemed to be University), Pune, Maharashtra, India

*For research use only. Not intended for diagnosis, treatment decisions, or clinical deployment.*
