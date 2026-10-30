"""Correct empty Murcko scaffold flags without changing molecule predictions.

Usage:
    python code/correct_scaffolds.py --screen results/v2/reliability_screening/natural_products_with_intervals.csv \
        --train results/legacy_benchmark/curated_activities.csv \
        --out results/v2/reliability_screening/natural_products_scaffolds_corrected.csv

Acyclic structures do not have a Bemis-Murcko ring framework; their status
is NO_MURCKO_FRAMEWORK, not SHARED_FRAMEWORK.
"""
import argparse
from pathlib import Path
import pandas as pd


def correct_scaffolds(screen: pd.DataFrame, train: pd.DataFrame) -> pd.DataFrame:
    if "Scaffold" not in screen or "Scaffold" not in train:
        raise ValueError("Both screening and training tables need a Scaffold column")
    data = screen.copy()
    valid = data["Scaffold"].notna() & data["Scaffold"].astype(str).str.strip().ne("")
    known = set(train["Scaffold"].dropna().astype(str).str.strip()) - {""}
    data["Scaffold_status_full_curated_CHEMBL203"] = "NO_MURCKO_FRAMEWORK"
    data.loc[valid, "Scaffold_status_full_curated_CHEMBL203"] = "NOVEL_FRAMEWORK"
    data.loc[valid & data["Scaffold"].isin(known), "Scaffold_status_full_curated_CHEMBL203"] = "SHARED_FRAMEWORK"
    if "Scaffold_seen_in_train" in data:
        reported = data["Scaffold_seen_in_train"].fillna(False).astype(bool)
        data["Scaffold_status_model_training_subset"] = "NO_MURCKO_FRAMEWORK"
        data.loc[valid & reported, "Scaffold_status_model_training_subset"] = "SHARED_FRAMEWORK"
        data.loc[valid & ~reported, "Scaffold_status_model_training_subset"] = "NOVEL_FRAMEWORK"
        corrected = pd.Series(pd.NA, index=data.index, dtype="boolean")
        corrected.loc[valid] = reported.loc[valid].to_numpy()
        data["Scaffold_seen_in_train_CORRECTED"] = corrected
    return data


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--screen",type=Path,required=True)
    p.add_argument("--train",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=correct_scaffolds(pd.read_csv(a.screen),pd.read_csv(a.train))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(a.out,index=False)
    print("Corrected scaffold categories:", result["Scaffold_status_full_curated_CHEMBL203"].value_counts().to_dict())
    print("Written:", a.out)

if __name__ == "__main__":main()
