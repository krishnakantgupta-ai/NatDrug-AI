#!/usr/bin/env python3
"""NatDrug-AI v2: source-aware research modules for EGFR (CHEMBL203).

All activity targets must be experimental. Never label IMPPAT source-only rows as
independent external validation without independently confirmed assay provenance.
This module complements the pre-existing NatDrug-AI benchmark scripts.
"""
from __future__ import annotations
import argparse, hashlib, json, os, re, sys, time, warnings
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger, rdBase
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold, GroupShuffleSplit
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from benchmark_natdrug_ai import standardize, features, read_curate
RDLogger.DisableLog('rdApp.warning')

SOURCE_INFO = {"target_chembl":"CHEMBL203", "human_uniprot":"P00533", "source":"IMPPAT 3.0", "release":"2026-09-30", "url":"https://cb.imsc.res.in/imppat/download"}

def save_json(data,path):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(data,indent=2,default=lambda o:float(o) if isinstance(o,np.floating) else str(o)),encoding='utf8')

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as fh:
        for block in iter(lambda:fh.read(1<<20),b''): h.update(block)
    return h.hexdigest()

def get_data(train_csv,curated_csv=None):
    if curated_csv and Path(curated_csv).exists():
        d=pd.read_csv(curated_csv).fillna({'Scaffold':''})
        if len(d)<100:raise RuntimeError('Curated training data unexpectedly small')
        return d
    return read_curate(train_csv,'parent','pooled')[0]

def metrics(y,p):
    y=np.asarray(y,float);p=np.asarray(p,float)
    d={'N':int(len(y)),'MAE_pIC50':float(mean_absolute_error(y,p)),
       'RMSE_pIC50':float(np.sqrt(mean_squared_error(y,p)))}
    d['R2']=float(r2_score(y,p)) if len(y)>1 else None
    return d

def file_info(root):
    root=Path(root)
    names=['data/CHEMBL203_IC50.csv','sources/imppat_chemical_normalized.csv','sources/imppat_bioactivity.tsv']
    return {str(x):{'sha256':digest(root/x),'bytes':(root/x).stat().st_size} for x in names if (root/x).exists()}

def make_provenance_audit(root):
    root=Path(root);out=root/'results/v2/provenance';out.mkdir(parents=True,exist_ok=True)
    p=root/'results/legacy_imppat_audit/external_egfr_candidates_NOT_YET_VALIDATED.csv'
    if not p.exists():raise FileNotFoundError('Run the IMPPAT preparation stage first: '+str(p))
    d=pd.read_csv(p)
    info=root/'sources/imppat_bioactivity.tsv'
    raw=pd.read_csv(info,sep='\t',dtype=str,low_memory=False)
    human=raw[(raw.get('UniProt_identifier',pd.Series('',index=raw.index)).fillna('').eq('P00533')) &
              (raw.get('Target_organism',pd.Series('',index=raw.index)).fillna('').str.contains('Homo sapiens|human',case=False,na=False))].copy()
    cols=['Source_ID','Name','Canonical_SMILES','IC50_nM','Original_Source','Raw_IC50','Source_Review_Status']
    d=d[cols].copy().rename(columns={'Source_ID':'IMPPAT_ID'})
    d['Reported_Database_Not_Primary_Reference']=d.Original_Source
    d['Original_DOI']='';d['Original_PMID']='';d['Original_Assay_ID']=''
    d['Verified_Independent']=False;d['Validation_Status']='UNVERIFIED_ASSAY_OR_PUBLICATION_PROVENANCE'
    d.to_csv(out/'natural_product_candidates_NEED_VERIFICATION.csv',index=False)
    mandatory=['IMPPAT_ID','Original_DOI','Original_PMID','Original_Assay_ID','Original_Measurement_ID',
               'Assay_Type','Endpoint_Compatible','Verified_Not_In_Training_Publications','Verification_Evidence_URL']
    template=d[['IMPPAT_ID','Name']].drop_duplicates().copy()
    for c in mandatory[1:]:template[c]=''
    template.to_csv(out/'provenance_review_TEMPLATE.csv',index=False)
    summary={'status':'SOURCE_CANDIDATES_ONLY_NOT_EXTERNAL_VALIDATION','candidate_rows':len(d),
        'candidate_unique_compounds':int(d.Canonical_SMILES.nunique()),
        'all_bioactivity_rows':len(raw),'human_EGFR_source_rows':len(human),
        'missing_original_DOI_PMID_assay_ids':True,
        'manual_review_required':['original per-measurement DOI/PMID and assay ID','activity assay comparability',
                                  'shared assay or paper with CHEMBL203 training','duplicate values across NPASS/BindingDB'],
        'input_files':file_info(root),'annotation':'Exact compound separation is NOT assay or publication independence.'}
    save_json(summary,out/'source_audit_status.json')
    return summary

def fetch_chembl_assay_metadata(root,limit=None,timeout=12):
    """Optional official ChEMBL API fetch, with persistent caching; does not verify IMPPAT refs."""
    import requests
    root=Path(root);out=root/'results/v2/provenance';out.mkdir(parents=True,exist_ok=True)
    train=pd.read_csv(root/'data/CHEMBL203_IC50.csv')
    assays=sorted(set(train.Assay_ID.dropna().astype(str)))
    if limit is not None:assays=assays[:int(limit)]
    cache=out/'chembl_assays_cache';cache.mkdir(exist_ok=True)
    result=[]
    session=requests.Session();session.headers.update({'User-Agent':'NatDrug-AI-v2-reproducibility/0.2'})
    for i,assay in enumerate(assays,1):
        p=cache/(assay+'.json')
        if p.exists():
            item=json.loads(p.read_text())
        else:
            try:
                resp=session.get('https://www.ebi.ac.uk/chembl/api/data/assay/'+assay+'.json',timeout=timeout)
                resp.raise_for_status();item=resp.json();save_json(item,p)
            except Exception as ex:
                item={'assay_chembl_id':assay,'fetch_error':repr(ex)}
                # Do not cache transient API failures permanently.
        result.append({'Assay_ID':assay,'Doc_ChEMBL_ID':item.get('doc_chembl_id',''),
                       'Assay_Type':item.get('assay_type',''), 'Description':item.get('description',''),
                       'Fetch_Error':item.get('fetch_error','')})
        if i%25==0:print('ChEMBL assay metadata:',i,'/',len(assays),flush=True)
    pd.DataFrame(result).to_csv(out/'training_assay_metadata.csv',index=False)
    summ={'queried':len(assays),'complete_total_assays':len(set(train.Assay_ID.dropna().astype(str))),
          'failed_queries':sum(bool(r['Fetch_Error']) for r in result),
          'status':'ASSAY_ID_AND_DOC_ID_ONLY; retrieve publication metadata and check candidate provenance before independent validation'}
    save_json(summ,out/'chembl_assay_fetch_status.json');return summ

def compare_additional_models(root,force=False):
    root=Path(root);out=root/'results/v2/model_comparison';out.mkdir(parents=True,exist_ok=True)
    csv=out/'additional_model_metrics.csv'
    if csv.exists() and not force:return {'status':'CACHED','file':str(csv)}
    d=get_data(root/'data/CHEMBL203_IC50.csv',root/'results/legacy_benchmark/curated_activities.csv')
    x,_=features(d.Canonical_SMILES)
    f=x['fp_desc'];y=d.pIC50.to_numpy(float);groups=d.Scaffold.fillna('').to_numpy()
    if len(np.unique(groups))<5:raise RuntimeError('Insufficient scaffold groups')
    rows=[];oof=[]
    for name in ('ExtraTrees_fp_desc',):
        pred=np.full(len(d),np.nan)
        for fold,(tr,te) in enumerate(GroupKFold(n_splits=5).split(f,y,groups),1):
            model=(ExtraTreesRegressor(n_estimators=80,min_samples_leaf=2,max_features=.5,random_state=42,n_jobs=-1)
                   if name.startswith('ExtraTrees') else
                   RandomForestRegressor(n_estimators=150,min_samples_leaf=2,max_features=.5,random_state=52,n_jobs=-1))
            model.fit(f[tr],y[tr]);p=model.predict(f[te]);pred[te]=p
            for idx,v in zip(te,p):oof.append({'Model':name,'Fold':fold,'SMILES':d.iloc[idx].Canonical_SMILES,
                                 'Observed_pIC50':y[idx],'Predicted_pIC50':v})
            print('Model comparison:',name,'fold',fold,'of 5',flush=True)
        rows.append({'Model':name,'Split':'scaffold',**metrics(y,pred)})
        pd.DataFrame(rows).to_csv(csv,index=False)
        pd.DataFrame(oof).to_csv(out/'additional_model_oof.csv',index=False)
    pd.DataFrame(rows).to_csv(csv,index=False);pd.DataFrame(oof).to_csv(out/'additional_model_oof.csv',index=False)
    save_json({'status':'COMPLETE','cross_validation':'5-fold scaffold GroupKFold',
               'model_selection_note':'Comparing models on CV folds does not establish independent performance. No external candidates used to tune the models.',
               'training_input_sha256':digest(root/'data/CHEMBL203_IC50.csv')},out/'model_comparison_manifest.json')
    return {'status':'COMPLETE','metrics':rows}

def reliability_screen(root,alpha=.1,force=False):
    """Use disjoint scaffold training/calibration/test; output empirical test coverage.
    Screening is with the same training-only model used to form calibration residuals.
    Conformal intervals are not guaranteed to cover out-of-distribution natural products.
    """
    root=Path(root);out=root/'results/v2/reliability_screening';out.mkdir(parents=True,exist_ok=True)
    ranked=out/'natural_products_with_intervals.csv'
    if ranked.exists() and not force:return {'status':'CACHED','ranked_file':str(ranked)}
    d=get_data(root/'data/CHEMBL203_IC50.csv',root/'results/legacy_benchmark/curated_activities.csv')
    libpath=root/'results/legacy_imppat_audit/natural_product_screening_library.csv'
    if not libpath.exists():raise FileNotFoundError('IMPPAT library missing: '+str(libpath))
    lib=pd.read_csv(libpath,low_memory=False)
    cached=root/'results/legacy_screening/natural_products_ranked_all.csv'
    if cached.exists():
        # Reuse molecule-by-molecule structures from the prior actual IMPPAT screening,
        # only if compound IDs match the original 18,030-compound source export.
        prior=pd.read_csv(cached,low_memory=False)
        if {'Source_ID','Canonical_SMILES'}.issubset(prior.columns) and set(prior.Source_ID.astype(str))==set(lib.Source_ID.astype(str)):
            lib=lib.merge(prior[['Source_ID','Canonical_SMILES']],on='Source_ID',how='left',validate='one_to_one')
            print('Reused canonical structures from audited previous IMPPAT screening',flush=True)
        else:lib['Canonical_SMILES']=lib.SMILES.map(lambda s:standardize(s,'parent'))
    else:lib['Canonical_SMILES']=lib.SMILES.map(lambda s:standardize(s,'parent'))
    lib=lib.dropna(subset=['Canonical_SMILES']).drop_duplicates('Canonical_SMILES').copy()
    lib=lib[~lib.Canonical_SMILES.isin(set(d.Canonical_SMILES))].copy()
    lib['Scaffold']=lib.Canonical_SMILES.map(lambda s:MurckoScaffold.MurckoScaffoldSmiles(smiles=s))
    if len(lib)<1:raise RuntimeError('Empty library after overlap exclusion')
    f,fps=features(d.Canonical_SMILES)
    y=d.pIC50.to_numpy(float);g=d.Scaffold.fillna('').to_numpy()
    # First reserve untouched scaffold-separated test set; then independent calibration scaffolds.
    partition=GroupShuffleSplit(n_splits=1,test_size=.15,random_state=2026)
    develop,test=next(partition.split(f['fp_desc'],y,g))
    sub=GroupShuffleSplit(n_splits=1,test_size=.18,random_state=2027)
    reltr,relcal=next(sub.split(f['fp_desc'][develop],y[develop],g[develop]))
    train,cal=develop[reltr],develop[relcal]
    trsc=set(g[train]);calsc=set(g[cal]);testsc=set(g[test])
    assert not (trsc&calsc or trsc&testsc or calsc&testsc)
    mdl=RandomForestRegressor(n_estimators=250,min_samples_leaf=2,max_features=.5,random_state=42,n_jobs=-1)
    print('Fitting reliability model on scaffold-disjoint training partition:',len(train),'compounds',flush=True)
    mdl.fit(f['fp_desc'][train],y[train]);pcal=mdl.predict(f['fp_desc'][cal]);ptest=mdl.predict(f['fp_desc'][test])
    resid=np.abs(y[cal]-pcal);n=len(resid)
    level=min(1.,np.ceil((n+1)*(1-alpha))/n)
    q=float(np.quantile(resid,level,method='higher'))
    coverage=float(np.mean((y[test]>=ptest-q)&(y[test]<=ptest+q)))
    e=metrics(y[test],ptest)
    tst=pd.DataFrame({'SMILES':d.iloc[test].Canonical_SMILES.to_numpy(),'Observed_pIC50':y[test],
      'Predicted_pIC50':ptest,'Lower_pIC50':ptest-q,'Upper_pIC50':ptest+q,'Scaffold':g[test]})
    tst.to_csv(out/'heldout_scaffold_interval_evaluation.csv',index=False)
    xlib,libfps=features(lib.Canonical_SMILES)
    pred=mdl.predict(xlib['fp_desc'])
    trainfps=[fps[i] for i in train]
    similarities=np.array([max(DataStructs.BulkTanimotoSimilarity(p,trainfps)) for p in libfps])
    lib['Predicted_pIC50']=pred;lib['Predicted_IC50_nM']=10**np.clip(9-pred,-15,15)
    lib['Lower_pIC50']=pred-q;lib['Upper_pIC50']=pred+q
    lib['Max_train_Tanimoto']=similarities
    lib['Flag_under_0p4_exploratory']=similarities<.4
    # Acyclic compounds have an empty Murcko scaffold, which is not evidence
    # that they share a ring framework with training compounds.
    lib['Scaffold_seen_in_train']=(lib.Scaffold.ne('') & lib.Scaffold.notna() & lib.Scaffold.isin(trsc))
    lib['Scaffold_status_model_training_subset']=np.where(
        lib.Scaffold.isna() | lib.Scaffold.eq(''), 'NO_MURCKO_FRAMEWORK',
        np.where(lib.Scaffold_seen_in_train,'SHARED_FRAMEWORK','NOVEL_FRAMEWORK'))
    # Evidence-based error table from scaffold-held-out predictions, not a validated AD classifier.
    tst_features=[fps[i] for i in test]
    tst['Max_train_Tanimoto']=[max(DataStructs.BulkTanimotoSimilarity(fp,trainfps)) for fp in tst_features]
    tst['Abs_error']=abs(tst.Observed_pIC50-tst.Predicted_pIC50)
    bins=[0,.2,.4,.6,.8,1.00001]
    tst['Similarity_bin']=pd.cut(tst.Max_train_Tanimoto,bins=bins,include_lowest=True)
    (tst.groupby('Similarity_bin',observed=True).agg(N=('Abs_error','size'),MAE=('Abs_error','mean'),
       Median_Error=('Abs_error','median')).reset_index()).to_csv(out/'test_error_vs_similarity.csv',index=False)
    lib=lib.sort_values('Predicted_pIC50',ascending=False)
    lib.to_csv(ranked,index=False)
    summary={'status':'COMPLETED_PROSPECTIVE_IN_SILICO_SCREENING_ONLY',
        'splits':{'training_n':len(train),'calibration_n':len(cal),'heldout_test_n':len(test),
                  'training_scaffolds':len(trsc),'calibration_scaffolds':len(calsc),'heldout_scaffolds':len(testsc)},
        'calibration_n':n,'nominal_coverage':1-alpha,'conformal_q_pIC50':q,
        'heldout_scaffold_empirical_coverage':coverage,'heldout_scaffold_metrics':e,
        'library_n':len(lib),'library_scaffolds':int(lib.Scaffold.replace('',np.nan).nunique()),
        'library_scaffold_status_model_training_subset_counts':lib.Scaffold_status_model_training_subset.value_counts().to_dict(),
        'library_fraction_similarity_lt_0p4':float(np.mean(similarities<.4)),
        'important_warning':'Nominal split-conformal coverage assumes exchangeability, which may fail for OOD natural products. Heldout coverage is from ChEMBL scaffolds, NOT external natural-product validation. 0.4 is an exploratory cutoff, not a calibrated confidence classification.',
        'input_checksums':file_info(root),'model':'RF 250 trees, 8 descriptors + Morgan fingerprints, training-only split; no natural-product data for tuning'}
    save_json(summary,out/'reliability_manifest.json')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6,4));ax.hist(lib.Predicted_pIC50,bins=45)
    ax.set(xlabel='Predicted pIC50 (in silico)',ylabel='Natural products',title='IMPPAT natural-product screening')
    fig.tight_layout();fig.savefig(out/'screening_pIC50_distribution.png',dpi=200);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4));ax.scatter(tst.Max_train_Tanimoto,tst.Abs_error,s=9,alpha=.4)
    ax.set(xlabel='Maximum training Tanimoto',ylabel='Absolute pIC50 error',title='Held-out ChEMBL scaffold test')
    fig.tight_layout();fig.savefig(out/'similarity_error_test.png',dpi=200);plt.close(fig)
    print('Completed library screening for',len(lib),'source structures; intervals from separate scaffold calibration; held-out empirical coverage',round(coverage,3),flush=True)
    return summary

def external_validation_gate(root,verified_csv=None):
    """Require original provenance, compatibility and clear study-level non-leakage.
    No automatic assertion of independence is made from IMPPAT database tags alone.
    """
    root=Path(root);out=root/'results/v2/external_validation';out.mkdir(parents=True,exist_ok=True)
    cp=root/'results/v2/provenance/natural_product_candidates_NEED_VERIFICATION.csv'
    if not cp.exists():make_provenance_audit(root)
    candidates=pd.read_csv(cp)
    status={'candidate_compounds':int(candidates.Canonical_SMILES.nunique()),
            'status':'BLOCKED_UNVERIFIED_SOURCE_PUBLICATION_ASSAY_METADATA',
            'reason':'IMPPAT Source names are not primary assay/publication IDs; ChEMBL training dataset has only assay IDs and no primary document IDs.',
            'external_validation_metrics':None}
    if verified_csv and Path(verified_csv).exists():
        v=pd.read_csv(verified_csv).fillna('');required=['IMPPAT_ID','Original_Assay_ID','Original_Measurement_ID','Assay_Type',
            'Endpoint_Compatible','Verified_Not_In_Training_Publications','Verification_Evidence_URL']
        missing=[k for k in required if k not in v]
        if missing:status['reason']='Verification CSV lacks mandatory fields: '+','.join(missing)
        else:
            v['source_verified']=v[['Original_Assay_ID','Original_Measurement_ID','Assay_Type','Verification_Evidence_URL']].astype(str).apply(lambda col:col.str.strip().ne('')).all(axis=1)
            v['external_claim_eligible']=v.source_verified & v.Endpoint_Compatible.astype(str).str.lower().isin(['true','yes']) & v.Verified_Not_In_Training_Publications.astype(str).str.lower().isin(['true','yes']) & ((v.get('Original_DOI',pd.Series('',index=v.index)).astype(str).str.strip().ne('')) | (v.get('Original_PMID',pd.Series('',index=v.index)).astype(str).str.strip().ne('')))
            joined=candidates.merge(v,on='IMPPAT_ID',how='left',suffixes=('','_verified'))
            joined.to_csv(out/'eligibility_review.csv',index=False)
            eligible=joined[joined.external_claim_eligible.fillna(False)]
            status['reviewed_rows']=len(v);status['eligible_rows_after_user_attestation']=len(eligible)
            # Avoid false 'fully verified' because the training assay metadata is incomplete.
            trainingmeta=root/'results/v2/provenance/chembl_assay_fetch_status.json'
            if trainingmeta.exists():
                md=json.loads(trainingmeta.read_text());status['training_assays_metadata_complete']=(md['queried']==md['complete_total_assays'] and md['failed_queries']==0)
            else:status['training_assays_metadata_complete']=False
            status['reason']='User-provided verification requires independent checking of primary papers and cross-source assay/document identifiers. No external metrics computed automatically without complete assay/publication verification.'
    save_json(status,out/'external_validation_status.json');return status

def build_report(root):
    root=Path(root);out=root/'results/v2/report';out.mkdir(parents=True,exist_ok=True)
    b=root/'results/legacy_benchmark/summary_metrics.csv';s=root/'results/legacy_robustness/mlp_3seed_summary.csv'
    lines=['# NatDrug-AI v2 — Evidence-first computational study report','',
       'Report generated from supplied measured datasets and machine-generated outputs. **Not a manuscript-ready external validation claim.**','',
       '## Source and reproducibility','',f'- ChEMBL source file SHA256: `{digest(root/"data/CHEMBL203_IC50.csv")}`',
       f'- IMPPAT bioactivity table SHA256: `{digest(root/"sources/imppat_bioactivity.tsv")}`',
       '- Human EGFR: ChEMBL target CHEMBL203, UniProt P00533.','- IMPPAT 3.0 release: 2026-09-30.','',
       '## Experiment 1: measured ChEMBL assay benchmark','']
    if b.exists():
        d=pd.read_csv(b);lines.append(d.round(4).to_markdown(index=False));lines.append('')
        lines.append('Numbers above are cached actual previous-Mac-run results, NOT newly computed by this report builder.');lines.append('')
    else:lines.append('No benchmark file present.');lines.append('')
    lines+=['## Experiment 4: MLP stability','']
    if s.exists():lines += [pd.read_csv(s).round(4).to_markdown(index=False),'']
    else:lines+=['Full MLP study pending.','']
    for h,p in [('Alternative models','results/v2/model_comparison/additional_model_metrics.csv')]:
        lines+=['## '+h,''];q=root/p
        lines+=[pd.read_csv(q).round(4).to_markdown(index=False),''] if q.exists() else ['Not rerun in this package.','']
    reli=root/'results/v2/reliability_screening/reliability_manifest.json'
    lines+=['## Experiment 3: natural-product virtual screening + uncertainty','']
    if reli.exists():
        x=json.loads(reli.read_text());lines+= [f'- In-silico-screened IMPPAT structures: {x["library_n"]:,}',
            f'- ChEMBL scaffold-held-out interval coverage (NOT natural-product external coverage): {x["heldout_scaffold_empirical_coverage"]:.3f}',
            f'- Interval half-width: {x["conformal_q_pIC50"]:.3f} pIC50',
            f'- Screening structures below heuristic Tanimoto 0.4: {x["library_fraction_similarity_lt_0p4"]:.1%}',
            '- Intervals do not have guaranteed coverage on novel natural-product chemistry.','']
    else:lines+=['Calibrated screening pending; previous complete screening remains available as legacy data.','']
    ex=root/'results/v2/external_validation/external_validation_status.json'
    lines+=['## Experiment 2: externally measured natural-product EGFR IC50','']
    if ex.exists():
        x=json.loads(ex.read_text());lines += [f'- Status: **{x["status"]}**',f'- Candidate structures: {x["candidate_compounds"]}',f'- Reason: {x["reason"]}','']
    else:lines+=['Provenance gate not yet evaluated.','']
    lines+=['## Interpretation and submission blockers','',
       '1. Exact SMILES non-overlap does not rule out shared assay or publication sources.',
       '2. IMPPAT NPASS / BindingDB pairs may contain the same reported experimental measurements.',
       '3. Independent EGFR natural-product performance remains **unverified**; do not call 13-row diagnostic metrics an external validation.',
       '4. Prediction intervals calibrated on ChEMBL scaffolds may under-cover structurally novel natural products.',
       '5. Freeze all datasets, seed/split assignments, source release and environment before submission.',
       '6. A negative natural-product diagnostic result must be transparently reported; model improvements cannot be chosen using an untouched external test.','']
    (out/'RESEARCH_SUMMARY.md').write_text('\n'.join(lines),encoding='utf8')
    save_json({'UTC_time':datetime.now(timezone.utc).isoformat(),'inputs':file_info(root),'python':sys.version,'rdkit':rdBase.rdkitVersion,
               'claim_policy':'No validated natural-product external accuracy without primary-paper and assay matching.'},out/'reproducibility_manifest.json')
    return str(out/'RESEARCH_SUMMARY.md')

def main():
    pa=argparse.ArgumentParser();pa.add_argument('action',choices=['audit','fetch_assays','compare','reliability','external','report','all'])
    pa.add_argument('--root',default=str(Path(__file__).resolve().parents[1]));pa.add_argument('--verified-csv')
    pa.add_argument('--limit',type=int);pa.add_argument('--force',action='store_true')
    a=pa.parse_args();root=Path(a.root)
    jobs={'audit':lambda:make_provenance_audit(root),'fetch_assays':lambda:fetch_chembl_assay_metadata(root,a.limit),
          'compare':lambda:compare_additional_models(root,force=a.force),
          'reliability':lambda:reliability_screen(root,force=a.force),
          'external':lambda:external_validation_gate(root,a.verified_csv),'report':lambda:build_report(root)}
    if a.action=='all':
        for job in ['audit','compare','reliability','external','report']:
            print('RUNNING',job,flush=True);print(jobs[job](),flush=True)
    else:print(jobs[a.action]())
if __name__=='__main__':main()
