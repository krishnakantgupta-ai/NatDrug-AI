"""Predeclared scaffold-separated robustness studies for EGFR pIC50.
Three seeds for full-batch MLP with early stopping and fold-level curves.
Assay-policy/standardization RF sensitivity uses the same 250-tree configuration.
"""
import argparse,json,time,hashlib,platform,sys
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error,mean_squared_error,r2_score
from benchmark_natdrug_ai import read_curate,features,EarlyStoppingMLP

def metrics(y,yp):
 return dict(MAE_pIC50=float(mean_absolute_error(y,yp)),RMSE_pIC50=float(np.sqrt(mean_squared_error(y,yp))),R2=float(r2_score(y,yp)))

def run_rf(csv,out,policies,folds,seed,trees=250):
 rows=[];foldrows=[];oof=[];counts={}
 for mode,assay in policies:
  key=f'{mode}__{assay}'
  print('RF sensitivity:',key,flush=True)
  d,cnt=read_curate(csv,mode,assay);counts[key]=cnt
  if len(d)<50 or d.Scaffold.nunique()<folds:
   rows.append({'Condition':key,'Status':'INSUFFICIENT_DATA','n':len(d),'scaffolds':d.Scaffold.nunique()});continue
  x,_=features(d.Canonical_SMILES);y=d.pIC50.to_numpy();groups=d.Scaffold.to_numpy();pred=np.full(len(d),np.nan)
  for k,(tr,te) in enumerate(GroupKFold(n_splits=folds).split(x['fp_desc'],y,groups),start=1):
   mdl=RandomForestRegressor(n_estimators=trees,min_samples_leaf=2,max_features=0.5,random_state=seed,n_jobs=-1)
   mdl.fit(x['fp_desc'][tr],y[tr]);p=mdl.predict(x['fp_desc'][te]);pred[te]=p
   foldrows.append({'Condition':key,'Fold':k,'N_test':len(te),**metrics(y[te],p)})
   oof.extend({'Condition':key,'SMILES':d.iloc[i].Canonical_SMILES,'Observed_pIC50':y[i],'Predicted_pIC50':p[j],'Fold':k} for j,i in enumerate(te))
   print('  fold',k,'done',flush=True)
  rows.append({'Condition':key,'Status':'COMPLETED','n':len(d),'scaffolds':d.Scaffold.nunique(),**metrics(y,pred)})
 pd.DataFrame(rows).to_csv(out/'curation_sensitivity_summary.csv',index=False)
 pd.DataFrame(foldrows).to_csv(out/'curation_sensitivity_folds.csv',index=False)
 pd.DataFrame(oof).to_csv(out/'curation_sensitivity_oof.csv',index=False)
 (out/'curation_counts.json').write_text(json.dumps(counts,indent=2))
 return rows

def run_mlp(csv,out,seeds,folds):
 import torch
 torch.set_num_threads(min(2,torch.get_num_threads()))
 d,_=read_curate(csv,'parent','pooled')
 x,_=features(d.Canonical_SMILES);x=x['fp'];y=d.pIC50.to_numpy();groups=d.Scaffold.to_numpy()
 rows=[];oof=[];curves=[]
 for seed in seeds:
  pred=np.full(len(d),np.nan)
  for k,(tr,te) in enumerate(GroupKFold(n_splits=folds).split(x,y,groups),start=1):
   mdl=EarlyStoppingMLP(seed=seed+k,patience=18,epochs=220)
   mdl.fit(x[tr],y[tr],groups[tr]);p=mdl.predict(x[te]);pred[te]=p
   print('MLP seed',seed,'fold',k,'epochs',mdl.epochs_trained,flush=True)
   curves.extend({'Seed':seed,'Fold':k,**r} for r in mdl.curve)
   oof.extend({'Seed':seed,'Fold':k,'SMILES':d.iloc[i].Canonical_SMILES,'Observed_pIC50':y[i],'Predicted_pIC50':p[j]} for j,i in enumerate(te))
  rows.append({'Seed':seed,'Scaffolds':d.Scaffold.nunique(),'Compounds':len(d),**metrics(y,pred)})
  # Save checkpoints between seeds to survive interruptions
  pd.DataFrame(rows).to_csv(out/'mlp_3seed_summary.csv',index=False)
  pd.DataFrame(oof).to_csv(out/'mlp_3seed_oof.csv',index=False)
  pd.DataFrame(curves).to_csv(out/'mlp_learning_curves.csv',index=False)
 return rows

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--csv',required=True);p.add_argument('--out',required=True);p.add_argument('--mode',choices=['all','mlp','curation'],default='all');p.add_argument('--seeds',default='42,52,62');p.add_argument('--folds',type=int,default=5)
 a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 provenance={'sha256':hashlib.sha256(Path(a.csv).read_bytes()).hexdigest(),'models':'MLP 2048-256-64-1, dropout 0.2, AdamW, early stopping 18 epochs with inner scaffold split; RF 250 trees, 0.5 max_features, leaf 2, five scaffold GroupKFold outer folds','seeds':a.seeds,'python':sys.version,'platform':platform.platform()}
 (out/'robustness_provenance.json').write_text(json.dumps(provenance,indent=2))
 if a.mode in ('all','curation'):run_rf(a.csv,out,[('parent','pooled'),('none','pooled'),('parent','dominant')],a.folds,42)
 if a.mode in ('all','mlp'):run_mlp(a.csv,out,list(map(int,a.seeds.split(','))),a.folds)
