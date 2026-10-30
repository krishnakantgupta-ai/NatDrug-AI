#!/usr/bin/env python3
"""Prospective *in silico* screening of an externally sourced natural-product library.

Outputs are computational hypotheses ONLY. The script does not assert measured
activity, efficacy, NP authenticity or externally verified model accuracy.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor
from benchmark_natdrug_ai import read_curate,standardize,features,GEN

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--train-csv',required=True,help='Same frozen measured IC50 data used for benchmark')
 p.add_argument('--library-csv',required=True,help='Sourced natural product CSV: Name,SMILES, Source_ID recommended')
 p.add_argument('--library-provenance',required=True,help='Database, version/date, license and download URL to record in manifest')
 p.add_argument('--out',default='results/NP_screening');p.add_argument('--seed',type=int,default=42)
 p.add_argument('--model',choices=['RF_fingerprint','RF_fingerprint_descriptors'],default='RF_fingerprint')
 p.add_argument('--standardization',choices=['none','parent'],default='parent')
 args=p.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
 d,c=read_curate(args.train_csv,args.standardization,'pooled');raw=pd.read_csv(args.library_csv)
 if not 'SMILES' in raw:raise ValueError('Library must contain SMILES')
 if 'Name' not in raw:raw['Name']=['NP_'+str(i+1) for i in range(len(raw))]
 raw['Canonical_SMILES']=raw.SMILES.map(lambda s:standardize(s,args.standardization))
 n_invalid=int(raw.Canonical_SMILES.isna().sum());valid=raw.dropna(subset=['Canonical_SMILES']).copy()
 n_duplicates=int(valid.Canonical_SMILES.duplicated().sum());valid=valid.drop_duplicates('Canonical_SMILES')
 training=set(d.Canonical_SMILES);n_overlap=int(valid.Canonical_SMILES.isin(training).sum())
 valid=valid[~valid.Canonical_SMILES.isin(training)].copy()
 if valid.empty:raise ValueError('No valid independent, non-duplicate library molecules after exclusions')
 x,trainfps=features(d.Canonical_SMILES);xt,testfps=features(valid.Canonical_SMILES)
 rep='fp' if args.model=='RF_fingerprint' else 'fp_desc'
 mdl=RandomForestRegressor(n_estimators=250,min_samples_leaf=2,max_features='sqrt' if rep=='fp' else .5,random_state=args.seed,n_jobs=-1)
 mdl.fit(x[rep],d.pIC50)
 pred=mdl.predict(xt[rep]);sim=np.array([max(DataStructs.BulkTanimotoSimilarity(fp,trainfps)) for fp in testfps])
 nps=set(d.Scaffold)
 from rdkit.Chem.Scaffolds import MurckoScaffold
 valid['Scaffold']=valid.Canonical_SMILES.map(lambda s:MurckoScaffold.MurckoScaffoldSmiles(smiles=s))
 valid['Predicted_pIC50']=np.round(pred,4)
 valid['Predicted_IC50_nM']=np.power(10., np.clip(9-pred,-15,15))
 valid['Max_Train_Tanimoto']=np.round(sim,4)
 valid['Scaffold_in_training']=valid.Scaffold.isin(nps)
 valid['MW']=valid.Canonical_SMILES.map(lambda s:Descriptors.MolWt(Chem.MolFromSmiles(s)))
 valid['LogP']=valid.Canonical_SMILES.map(lambda s:Crippen.MolLogP(Chem.MolFromSmiles(s)))
 valid['TPSA']=valid.Canonical_SMILES.map(lambda s:rdMolDescriptors.CalcTPSA(Chem.MolFromSmiles(s)))
 # No arbitrary reliability cutoff is treated as validated. Output full evidence.
 valid=valid.sort_values(['Predicted_pIC50'],ascending=False)
 valid.to_csv(out/'natural_products_ranked_all.csv',index=False)
 metadata={'library_provenance':args.library_provenance,'library_sha256':hashlib.sha256(Path(args.library_csv).read_bytes()).hexdigest(),
           'training_sha256':hashlib.sha256(Path(args.train_csv).read_bytes()).hexdigest(),
           'model':args.model,'seed':args.seed,'raw_rows':len(raw),'invalid_rows':n_invalid,
           'duplicate_structures':n_duplicates,'training_structure_overlap_excluded':n_overlap,
           'independent_molecules_screened':len(valid),'unique_scaffolds':int(valid.Scaffold.nunique()),
           'scaffolds_shared_with_training':int(valid.Scaffold_in_training.sum()),
           'mean_max_similarity':float(sim.mean()),'fraction_similarity_below_0_4_unvalidated_cutoff':float((sim<.4).mean()),
           'rdkit':rdBase.rdkitVersion,'warning':'Predictions are not biological evidence or external validation. Library NP origin is supplied and must be verified. Activity enrichment requires measured labels.'}
 (out/'screening_manifest.json').write_text(json.dumps(metadata,indent=2))
 # Save publication-oriented chemical-space and predicted activity distributions.
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 for field,title,xlabel in [('Predicted_pIC50','Prediction distribution (in silico only)','Predicted pIC50'),('Max_Train_Tanimoto','Chemical-space coverage of NP library','Nearest-training Tanimoto similarity')]:
  fig,ax=plt.subplots(figsize=(7,4.5))
  ax.hist(valid[field].dropna(),bins=40)
  ax.set(xlabel=xlabel,ylabel='Number of compounds',title=title)
  if field=='Max_Train_Tanimoto':ax.set_xlim(0,1)
  fig.tight_layout();fig.savefig(out/(field+'_distribution.png'),dpi=200);plt.close(fig)
 metadata['predicted_pIC50_median']=float(valid.Predicted_pIC50.median())
 metadata['predicted_pIC50_q25']=float(valid.Predicted_pIC50.quantile(.25))
 metadata['predicted_pIC50_q75']=float(valid.Predicted_pIC50.quantile(.75))
 metadata['median_max_similarity']=float(valid.Max_Train_Tanimoto.median())
 (out/'screening_manifest.json').write_text(json.dumps(metadata,indent=2))
 print(json.dumps(metadata,indent=2))
 print('Output',out.resolve())
if __name__=='__main__':main()
