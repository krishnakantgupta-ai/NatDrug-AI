#!/usr/bin/env python3
"""NatDrug-AI publication benchmark: measured IC50 only, leakage-safe CV.

Data provenance is supplied by the investigator. This script does not download
activities or generate synthetic labels. It will refuse scaffold CV when there
are insufficient scaffolds. No evaluation output is fabricated.
"""
from __future__ import annotations
import argparse, hashlib, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator, Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold, KFold
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import sklearn

GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
CLEAN = rdMolStandardize.CleanupParameters()

def standardize(s, mode):
    try:
        m = Chem.MolFromSmiles(str(s))
        if m is None or not m.GetNumHeavyAtoms(): return None
        if mode == 'parent':
            m = rdMolStandardize.Cleanup(m)
            m = rdMolStandardize.FragmentParent(m)
            m = rdMolStandardize.Uncharger().uncharge(m)
            m = rdMolStandardize.TautomerEnumerator().Canonicalize(m)
        return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)
    except Exception:
        return None

def read_curate(path, mode, assay_policy):
    raw = pd.read_csv(path, low_memory=False)
    if not {'SMILES','IC50_nM'}.issubset(raw.columns):
        raise ValueError('Need CSV columns SMILES,IC50_nM; optional Assay_ID and standard_relation')
    counts={'input_rows':len(raw)}
    d=raw.copy()
    if 'standard_relation' in d:
        d=d[d.standard_relation.astype(str).str.strip().eq('=')].copy()
        counts['exact_relation_rows']=len(d)
    d['IC50_nM']=pd.to_numeric(d.IC50_nM,errors='coerce')
    d=d[np.isfinite(d.IC50_nM) & (d.IC50_nM>0)].copy()
    counts['positive_numeric_exact_rows']=len(d)
    d['Canonical_SMILES']=d.SMILES.map(lambda s:standardize(s,mode))
    d=d.dropna(subset=['Canonical_SMILES']).copy()
    counts['valid_structures']=len(d)
    d['Assay_ID']=d.Assay_ID.astype(str) if 'Assay_ID' in d else 'unknown'
    counts['assays_before_filter']=int(d.Assay_ID.nunique())
    if assay_policy=='dominant':
        if d.empty: raise ValueError('No valid records')
        assay=d.Assay_ID.value_counts().index[0]
        d=d[d.Assay_ID.eq(assay)].copy()
        counts['dominant_assay']=assay
    d['pIC50']=9-np.log10(d.IC50_nM)
    # Retain source assay provenance with each aggregate, allowing audit.
    agg=d.groupby('Canonical_SMILES',as_index=False).agg(
        pIC50=('pIC50','median'),n_measurements=('pIC50','size'),
        assay_sd=('pIC50','std'),n_assays=('Assay_ID','nunique'),
        Assay_IDs=('Assay_ID',lambda z:';'.join(sorted(set(z)))))
    agg['Scaffold']=agg.Canonical_SMILES.map(lambda s:MurckoScaffold.MurckoScaffoldSmiles(smiles=s,includeChirality=False))
    counts['curated_compounds']=len(agg)
    counts['distinct_scaffolds']=int(agg.Scaffold.nunique())
    counts['multi_assay_structures']=int((agg.n_assays>1).sum())
    counts['assays_after_filter']=int(d.Assay_ID.nunique())
    return agg, counts

def features(smiles):
    x=np.zeros((len(smiles),2048),dtype=np.float32)
    descriptors=np.zeros((len(smiles),8),dtype=np.float32)
    fps=[]
    for i,s in enumerate(smiles):
        m=Chem.MolFromSmiles(s);fp=GEN.GetFingerprint(m);fps.append(fp)
        DataStructs.ConvertToNumpyArray(fp,x[i])
        descriptors[i]=[Descriptors.MolWt(m),Crippen.MolLogP(m),rdMolDescriptors.CalcTPSA(m),
          Lipinski.NumHDonors(m),Lipinski.NumHAcceptors(m),Lipinski.NumRotatableBonds(m),
          rdMolDescriptors.CalcFractionCSP3(m),rdMolDescriptors.CalcNumRings(m)]
    return {'fp':x,'fp_desc':np.hstack([x,descriptors]),'desc':descriptors},fps

def make_model(name, seed, n_trees):
    if name=='median':return DummyRegressor(strategy='median'),'fp'
    if name=='ridge':return make_pipeline(StandardScaler(),Ridge(alpha=10.0,solver='lsqr')),'fp_desc'
    if name=='RF_fingerprint':return RandomForestRegressor(n_estimators=n_trees,min_samples_leaf=2,max_features='sqrt',random_state=seed,n_jobs=-1),'fp'
    if name=='RF_fingerprint_descriptors':return RandomForestRegressor(n_estimators=n_trees,min_samples_leaf=2,max_features=0.5,random_state=seed,n_jobs=-1),'fp_desc'
    if name=='MLP_early_stop':return None,'fp'
    raise ValueError(name)

class EarlyStoppingMLP:
    """Train on outer training fold; internal validation split by scaffold only."""
    def __init__(self, seed, patience=18, epochs=220):
        self.seed=seed;self.patience=patience;self.epochs=epochs;self.curve=[]
    def fit(self,x,y,groups):
        import torch
        torch.manual_seed(self.seed)
        unique=np.unique(groups)
        if len(unique)<3:raise ValueError('MLP needs >=3 distinct training scaffolds for internal validation')
        rng=np.random.default_rng(self.seed)
        pick=rng.permutation(unique)
        valgroups=set(pick[:max(1,int(np.ceil(.15*len(unique))) )])
        val=np.asarray([i for i,g in enumerate(groups) if g in valgroups]); tr=np.asarray([i for i,g in enumerate(groups) if g not in valgroups]);
        if len(tr)<5 or len(val)<2:raise ValueError('Not enough samples for stable MLP internal scaffold split')
        self.mean=float(y[tr].mean());self.sd=float(max(y[tr].std(),1e-6))
        t=lambda a:torch.from_numpy(np.asarray(a,dtype=np.float32))
        xx=t(x[tr]); xv=t(x[val]);yy=t(((y[tr]-self.mean)/self.sd).reshape(-1,1)); yv=t(((y[val]-self.mean)/self.sd).reshape(-1,1))
        import copy
        net=torch.nn.Sequential(torch.nn.Linear(x.shape[1],256),torch.nn.ReLU(),torch.nn.Dropout(.2),torch.nn.Linear(256,64),torch.nn.ReLU(),torch.nn.Linear(64,1))
        opt=torch.optim.AdamW(net.parameters(),lr=0.001,weight_decay=.01)
        best=np.inf;state=None;counter=0
        for ep in range(self.epochs):
            net.train();opt.zero_grad();loss=torch.nn.functional.mse_loss(net(xx),yy);loss.backward();opt.step()
            net.eval()
            with torch.no_grad():vl=float(torch.nn.functional.mse_loss(net(xv),yv).item())
            self.curve.append({'epoch':ep+1,'train_loss':float(loss.item()),'val_loss':vl})
            if vl<best-1e-5:best=vl;state=copy.deepcopy(net.state_dict());counter=0
            else:counter+=1
            if counter>=self.patience:break
        net.load_state_dict(state);net.eval();self.model=net
        self.epochs_trained=len(self.curve)
        return self
    def predict(self,x):
        import torch
        with torch.no_grad():
            return self.model(torch.tensor(np.asarray(x,dtype=np.float32))).numpy().ravel()*self.sd+self.mean

def score(y,p):
    return {'MAE_pIC50':float(mean_absolute_error(y,p)), 'RMSE_pIC50':float(np.sqrt(mean_squared_error(y,p))),
            'R2':float(r2_score(y,p)) if len(y)>1 else None}

def plot_outputs(oof,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for (split,model),d in oof.groupby(['Split','Model']):
        fig,ax=plt.subplots(figsize=(5.5,5))
        ax.scatter(d.Observed_pIC50,d.Predicted_pIC50,s=8,alpha=.35)
        lo=min(d.Observed_pIC50.min(),d.Predicted_pIC50.min())
        hi=max(d.Observed_pIC50.max(),d.Predicted_pIC50.max())
        ax.plot([lo,hi],[lo,hi],'k--',lw=1)
        ax.set(xlabel='Observed pIC50',ylabel='OOF predicted pIC50',title=f'{split} | {model}')
        fig.tight_layout();fig.savefig(out/f'scatter_{split}_{model}.png',dpi=200);plt.close(fig)
    d=oof[oof.Split.eq('scaffold')].copy()
    d['abs_error']=(d.Observed_pIC50-d.Predicted_pIC50).abs()
    fig,ax=plt.subplots(figsize=(6,4))
    for m,z in d.groupby('Model'):
        ax.scatter(z.Max_train_Tanimoto,z.abs_error,label=m,s=7,alpha=.2)
    ax.set(xlabel='Maximum similarity to training fold',ylabel='Absolute pIC50 error',title='Out-of-fold applicability analysis')
    ax.legend(fontsize=7,loc='upper right');fig.tight_layout()
    fig.savefig(out/'error_vs_similarity.png',dpi=220);plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--csv',required=True,help='Measured activity CSV (SMILES,IC50_nM, optional Assay_ID)')
    ap.add_argument('--target',default='CHEMBL203');ap.add_argument('--out',default='results/MolInf_benchmark')
    ap.add_argument('--standardization',choices=['none','parent'],default='parent')
    ap.add_argument('--assay-policy',choices=['pooled','dominant'],default='pooled')
    ap.add_argument('--folds',type=int,default=5);ap.add_argument('--seed',type=int,default=42)
    ap.add_argument('--trees',type=int,default=250);ap.add_argument('--include-mlp',action='store_true')
    ap.add_argument('--mlp-seeds',default='42,52,62',help='Use >=3 seeds for stability analysis; adds separate replicate result sets')
    ap.add_argument('--only-scaffold',action='store_true')
    args=ap.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    d,cnt=read_curate(args.csv,args.standardization,args.assay_policy)
    d.to_csv(out/'curated_activities.csv',index=False)
    if len(d)<max(50,args.folds):raise ValueError('Insufficient compounds for publication benchmark: require >=50')
    if d.Scaffold.nunique()<args.folds:raise ValueError('Insufficient distinct scaffolds; no optimistic random fallback')
    x,fps=features(d.Canonical_SMILES);y=d.pIC50.to_numpy();groups=d.Scaffold.to_numpy()
    splitters={'scaffold':GroupKFold(n_splits=args.folds)}
    if not args.only_scaffold:splitters['random']=KFold(n_splits=args.folds,shuffle=True,random_state=args.seed)
    models=['median','ridge','RF_fingerprint','RF_fingerprint_descriptors']
    if args.include_mlp:models+=['MLP_early_stop']
    rec=[];fold_metrics=[];curves=[]
    for splitname,splitter in splitters.items():
        splits=splitter.split(x['fp'],y,groups) if splitname=='scaffold' else splitter.split(x['fp'],y)
        for fold,(train,test) in enumerate(splits,1):
            sims=np.asarray([max(DataStructs.BulkTanimotoSimilarity(fps[i], [fps[j] for j in train])) for i in test])
            for model_name in models:
                seeds=[int(v) for v in args.mlp_seeds.split(',')] if model_name=='MLP_early_stop' else [args.seed]
                for seed in seeds:
                    mdl,rep=make_model(model_name,seed,args.trees)
                    if model_name=='MLP_early_stop':
                        mdl=EarlyStoppingMLP(seed).fit(x[rep][train],y[train],groups[train])
                        for row in mdl.curve:curves.append({'Split':splitname,'Fold':fold,'Seed':seed,**row})
                    else: mdl.fit(x[rep][train],y[train])
                    pred=mdl.predict(x[rep][test]);m=f'{model_name}_seed{seed}' if model_name=='MLP_early_stop' else model_name
                    fold_metrics.append({'Split':splitname,'Fold':fold,'Model':m,'Test_N':len(test),'Train_N':len(train),**score(y[test],pred)})
                    for i,p,sim in zip(test,pred,sims):
                        rec.append({'Split':splitname,'Fold':fold,'Model':m,'SMILES':d.Canonical_SMILES[i],
                          'Scaffold':groups[i],'Observed_pIC50':float(y[i]),'Predicted_pIC50':float(p),
                          'Max_train_Tanimoto':float(sim),'Train_n':len(train),'Test_n':len(test)})
            print(f'Completed {splitname} fold {fold}/{args.folds}',flush=True)
    oof=pd.DataFrame(rec);oof.to_csv(out/'oof_predictions.csv',index=False)
    pd.DataFrame(fold_metrics).to_csv(out/'per_fold_metrics.csv',index=False)
    if curves:pd.DataFrame(curves).to_csv(out/'mlp_learning_curves.csv',index=False)
    summ=[]
    for (split,model),z in oof.groupby(['Split','Model']):
        summ.append({'Split':split,'Model':model,'Compounds':len(z),'Scaffolds':d.Scaffold.nunique(),**score(z.Observed_pIC50,z.Predicted_pIC50)})
    pd.DataFrame(summ).to_csv(out/'summary_metrics.csv',index=False)
    bins=[0,.2,.4,.6,.8,1.00001]
    oof['abs_error']=(oof.Observed_pIC50-oof.Predicted_pIC50).abs()
    oof['Similarity_bin']=pd.cut(oof.Max_train_Tanimoto,bins,include_lowest=True).astype(str)
    oof.groupby(['Split','Model','Similarity_bin'],observed=True).agg(n=('abs_error','size'),MAE_pIC50=('abs_error','mean')).reset_index().to_csv(out/'error_by_similarity_bin.csv',index=False)
    plot_outputs(oof,out)
    raw=Path(args.csv).read_bytes(); meta={'UTC_time':datetime.now(timezone.utc).isoformat(),'target':args.target,
      'input_file_sha256':hashlib.sha256(raw).hexdigest(),'source_release':'TO BE RECORDED by investigator',
      'python':sys.version,'platform':platform.platform(),'rdkit':rdBase.rdkitVersion,
      'scikit_learn':sklearn.__version__,'numpy':np.__version__,'pandas':pd.__version__,
      'arguments':vars(args),'curation_counts':cnt,
      'design':'Models use predeclared hyperparameters, not optimized on held-out folds. No random fallback. ' +
        'For claims of optimized models use nested CV. A natural-product library without labels is not an external accuracy test.'}
    (out/'provenance.json').write_text(json.dumps(meta,indent=2))
    print(pd.DataFrame(summ).to_string(index=False));print('Artifacts:',out.resolve())
if __name__=='__main__':main()
