#!/usr/bin/env python3
"""Nested scaffold CV for NatDrug-AI; no synthetic labels or test-fold tuning."""
import argparse, json, hashlib, platform, sys
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator, Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
import sklearn

def curate(df):
    if not {'SMILES','IC50_nM'}.issubset(df): raise ValueError('CSV must contain SMILES and IC50_nM')
    d=df.copy();d['IC50_nM']=pd.to_numeric(d.IC50_nM,errors='coerce')
    d=d[np.isfinite(d.IC50_nM)&(d.IC50_nM>0)].copy()
    def canonical(s):
        try:
            m=Chem.MolFromSmiles(str(s));return Chem.MolToSmiles(m) if m and m.GetNumHeavyAtoms()>0 else None
        except Exception:return None
    d['SMILES']=d.SMILES.map(canonical);d=d.dropna(subset=['SMILES'])
    d['pIC50']=9-np.log10(d.IC50_nM)
    d=d.groupby('SMILES',as_index=False).agg(pIC50=('pIC50','median'),Replicates=('pIC50','size'),Assay_SD=('pIC50','std'))
    d['Scaffold']=d.SMILES.map(lambda s:MurckoScaffold.MurckoScaffoldSmiles(smiles=s))
    return d

def make_features(smiles):
    fps={size:rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=size) for size in (2048,4096)}
    arrays={f'fp{size}':np.zeros((len(smiles),size),dtype=np.float32) for size in fps}
    desc=np.zeros((len(smiles),8),dtype=np.float32)
    for i,s in enumerate(smiles):
        m=Chem.MolFromSmiles(s)
        for size,gen in fps.items():DataStructs.ConvertToNumpyArray(gen.GetFingerprint(m),arrays[f'fp{size}'][i])
        desc[i]=[Descriptors.MolWt(m),Crippen.MolLogP(m),rdMolDescriptors.CalcTPSA(m),Lipinski.NumHDonors(m),Lipinski.NumHAcceptors(m),Lipinski.NumRotatableBonds(m),rdMolDescriptors.CalcFractionCSP3(m),rdMolDescriptors.CalcNumRings(m)]
    arrays['desc']=desc
    for size in fps:arrays[f'fp{size}_desc']=np.hstack((arrays[f'fp{size}'],desc))
    return arrays

def candidates():
    return [
        dict(name='RF_baseline',algo='rf',representation='fp2048',n_estimators=250,min_samples_leaf=2,max_features='sqrt'),
        dict(name='RF_optimized',algo='rf',representation='fp2048_desc',n_estimators=400,min_samples_leaf=2,max_features=0.5),
        dict(name='ExtraTrees',algo='et',representation='fp2048_desc',n_estimators=400,min_samples_leaf=2,max_features=0.5),
        dict(name='ExtraTrees_4096',algo='et',representation='fp4096_desc',n_estimators=400,min_samples_leaf=2,max_features=0.5),
        dict(name='MLR_OLS',algo='ols',representation='desc'),
        dict(name='MLR_Ridge_a1',algo='ridge',representation='fp2048_desc',alpha=1.0),
        dict(name='MLR_Ridge_a10',algo='ridge',representation='fp2048_desc',alpha=10.0),
        dict(name='MLR_Ridge_a100',algo='ridge',representation='fp2048_desc',alpha=100.0),
        dict(name='MLR_Lasso_a0p01',algo='lasso',representation='fp2048_desc',alpha=0.01),
        dict(name='MLR_ElasticNet_a0p01',algo='elastic',representation='fp2048_desc',alpha=0.01),
    ]

def estimator(c):
    algo=c['algo']
    if algo in ('rf','et'):
        cls=RandomForestRegressor if algo=='rf' else ExtraTreesRegressor
        return cls(n_estimators=c['n_estimators'],min_samples_leaf=c['min_samples_leaf'],max_features=c['max_features'],random_state=42,n_jobs=-1)
    if algo=='ols': return make_pipeline(StandardScaler(),LinearRegression())
    if algo=='ridge': return make_pipeline(StandardScaler(),Ridge(alpha=c['alpha'],solver='lsqr'))
    if algo=='lasso': return make_pipeline(StandardScaler(),Lasso(alpha=c['alpha'],max_iter=2500,selection='cyclic',random_state=42))
    if algo=='elastic': return make_pipeline(StandardScaler(),ElasticNet(alpha=c['alpha'],l1_ratio=0.5,max_iter=2500,random_state=42))
    raise ValueError(f'Unknown model: {algo}')

def score(y,p):return dict(MAE_pIC50=float(mean_absolute_error(y,p)),RMSE_pIC50=float(np.sqrt(mean_squared_error(y,p))),R2=float(r2_score(y,p)))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--csv',required=True);ap.add_argument('--target',default='CHEMBL203');ap.add_argument('--out',default='results/enhanced')
    ap.add_argument('--outer-folds',type=int,default=5);ap.add_argument('--inner-folds',type=int,default=3)
    ap.add_argument('--family',choices=['all','linear','trees'],default='all',help='all: compare linear and tree models; linear: only MLR variants; trees: original tree candidates')
    args=ap.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    raw=Path(args.csv).read_bytes();d=curate(pd.read_csv(args.csv));groups=d.Scaffold.to_numpy();y=d.pIC50.to_numpy()
    n_groups=len(np.unique(groups))
    if n_groups<max(3,args.outer_folds):raise ValueError('Insufficient scaffolds for requested scaffold-only outer CV; no random fallback')
    X=make_features(d.SMILES);cs=candidates()
    if args.family=='linear':cs=[c for c in cs if c['algo'] in ('ols','ridge','lasso','elastic')]
    if args.family=='trees':cs=[c for c in cs if c['algo'] in ('rf','et')]
    baseline=dict(name='RF_baseline',algo='rf',representation='fp2048',n_estimators=250,min_samples_leaf=2,max_features='sqrt')
    outer=GroupKFold(n_splits=args.outer_folds)
    selected=np.empty(len(d),dtype=object);fold_ids=np.zeros(len(d),dtype=int)
    oof=np.full(len(d),np.nan);base_oof=np.full(len(d),np.nan);selection=[];fold_metrics=[]
    for fold,(tr,te) in enumerate(outer.split(X['fp2048'],y,groups),1):
        inner_groups=groups[tr];inner_n=min(args.inner_folds,len(np.unique(inner_groups)))
        if inner_n<2:raise ValueError('Not enough training scaffolds for inner CV')
        inner=list(GroupKFold(n_splits=inner_n).split(X['fp2048'][tr],y[tr],inner_groups))
        trials=[]
        for c in cs:
            pred=np.full(len(tr),np.nan)
            for itr,ite in inner:
                m=estimator(c);m.fit(X[c['representation']][tr[itr]],y[tr[itr]])
                pred[ite]=m.predict(X[c['representation']][tr[ite]])
            trials.append((mean_squared_error(y[tr],pred),c))
        _,best=min(trials,key=lambda t:t[0]);m=estimator(best);m.fit(X[best['representation']][tr],y[tr]);oof[te]=m.predict(X[best['representation']][te])
        base=estimator(baseline);base.fit(X['fp2048'][tr],y[tr]);base_oof[te]=base.predict(X['fp2048'][te])
        fold_ids[te]=fold;selected[te]=best['name'];selection.append({'Fold':fold,'Selected':best['name'],'Inner_RMSE':float(np.sqrt(min(t[0] for t in trials))),'Test_N':len(te),'Test_scaffolds':len(np.unique(groups[te]))})
        for name,p in [('Nested_selected',oof[te]),('RF_baseline',base_oof[te])]:fold_metrics.append({'Fold':fold,'Model':name,**score(y[te],p)})
        print(f'Outer fold {fold}/{args.outer_folds}: selected {best["name"]}; completed',flush=True)
    preds=d[['SMILES','Scaffold','pIC50','Replicates','Assay_SD']].copy();preds['Fold']=fold_ids;preds['Selected_model']=selected;preds['Nested_predicted_pIC50']=oof;preds['RF_baseline_predicted_pIC50']=base_oof
    preds.to_csv(out/'out_of_fold_predictions.csv',index=False);pd.DataFrame(selection).to_csv(out/'model_selection.csv',index=False);pd.DataFrame(fold_metrics).to_csv(out/'fold_metrics.csv',index=False)
    metrics=pd.DataFrame([{'Target':args.target,'Model':name,'N':len(d),'Scaffolds':n_groups,'Validation':'Nested scaffold GroupKFold',**score(y,p)} for name,p in [('Nested_selected',oof),('RF_baseline',base_oof)]])
    metrics.to_csv(out/'validation_metrics.csv',index=False)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for name,p in [('Nested_selected',oof),('RF_baseline',base_oof)]:
        fig,ax=plt.subplots(figsize=(5.5,5.5));ax.scatter(y,p,s=9,alpha=.35);low=min(y.min(),p.min());high=max(y.max(),p.max());ax.plot([low,high],[low,high],'k--',lw=1)
        ax.set(xlabel='Experimental pIC50',ylabel='Out-of-fold predicted pIC50',title=f'{args.target}: {name}');fig.tight_layout();fig.savefig(out/f'{name}_scatter.png',dpi=250);plt.close(fig)
    (out/'provenance.json').write_text(json.dumps({'target':args.target,'sha256':hashlib.sha256(raw).hexdigest(),'python':sys.version,'rdkit':rdBase.rdkitVersion,'sklearn':sklearn.__version__,'outer_folds':args.outer_folds,'inner_folds':args.inner_folds,'candidate_models':cs,'baseline_model':baseline,'family':args.family,'notes':['Selection performed exclusively on outer training folds','No external validation or assay harmonization','Use same dataset and code release for manuscript comparisons']},indent=2))
    print('\n',metrics.to_string(index=False,float_format=lambda z:f'{z:.4f}'));print('Saved to',out.resolve())
if __name__=='__main__':main()
