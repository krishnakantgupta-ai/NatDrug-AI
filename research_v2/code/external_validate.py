#!/usr/bin/env python3
"""Blind, independently sourced measured-IC50 external evaluation.

Researcher must freeze the external dataset *before* model selection and must
provide its provenance. Uses fixed, prespecified RF architecture; no tuning on
external results. External compounds overlapping training are excluded/reported.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error,mean_squared_error,r2_score
from benchmark_natdrug_ai import read_curate,features,GEN
from rdkit import DataStructs

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--train-csv',required=True)
 p.add_argument('--test-csv',required=True)
 p.add_argument('--test-provenance',required=True,help='Published/independent source and acquisition date')
 p.add_argument('--target',default='CHEMBL203');p.add_argument('--out',default='results/external')
 p.add_argument('--model',choices=['RF_fingerprint','RF_fingerprint_descriptors'],default='RF_fingerprint')
 p.add_argument('--standardization',choices=['none','parent'],default='parent');p.add_argument('--seed',type=int,default=42)
 args=p.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
 if Path(args.train_csv).resolve()==Path(args.test_csv).resolve():raise ValueError('Cannot use the training file as independent test data')
 train,traincounts=read_curate(args.train_csv,args.standardization,'pooled')
 test,testcounts=read_curate(args.test_csv,args.standardization,'pooled')
 trainset=set(train.Canonical_SMILES);test['overlap_training']=test.Canonical_SMILES.isin(trainset)
 n_overlap=int(test.overlap_training.sum());test=test[~test.overlap_training].copy()
 if len(test)<10:raise ValueError('At least 10 nonoverlapping measured test structures required for this script')
 rep='fp' if args.model=='RF_fingerprint' else 'fp_desc'
 x,fps=features(train.Canonical_SMILES);xt,fpst=features(test.Canonical_SMILES)
 mdl=RandomForestRegressor(n_estimators=250,min_samples_leaf=2,max_features='sqrt' if rep=='fp' else 0.5,random_state=args.seed,n_jobs=-1)
 mdl.fit(x[rep],train.pIC50);pred=mdl.predict(xt[rep]);sim=[max(DataStructs.BulkTanimotoSimilarity(fp,fps)) for fp in fpst]
 test['Prediction_pIC50']=pred;test['Max_train_Tanimoto']=sim;test['Absolute_error_pIC50']=abs(test.pIC50-test.Prediction_pIC50)
 train_scafs=set(train.Scaffold)
 test['Scaffold_in_training']=test.Scaffold.isin(train_scafs)
 test.to_csv(out/'external_predictions.csv',index=False)
 result={'target':args.target,'test_provenance':args.test_provenance,'model':args.model,
  'n_train':len(train),'n_test':len(test),'excluded_exact_overlap':n_overlap,
  'test_scaffolds':int(test.Scaffold.nunique()),'scaffolds_shared_with_training':int(test.Scaffold_in_training.sum()),
  'MAE_pIC50':float(mean_absolute_error(test.pIC50,pred)),
  'RMSE_pIC50':float(np.sqrt(mean_squared_error(test.pIC50,pred))),
  'R2':float(r2_score(test.pIC50,pred)),
  'train_sha256':hashlib.sha256(Path(args.train_csv).read_bytes()).hexdigest(),
  'test_sha256':hashlib.sha256(Path(args.test_csv).read_bytes()).hexdigest(),
  'train_curation_counts':traincounts,'test_curation_counts':testcounts,
  'note':'Test dataset provenance and temporal independence must be justified by the investigator. The test data are never used for hyperparameter tuning.'}
 (out/'external_metrics.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
