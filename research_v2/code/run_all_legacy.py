"""Single entrypoint for reproducible EGFR publication experiments.
Runs main CV, natural-product extraction, screening, externally audited validation,
MLP stability and curation sensitivity. Never substitutes synthetic source data.
"""
import argparse,hashlib,json,os,subprocess,sys,traceback
from datetime import datetime,timezone
from pathlib import Path
from fetch_imppat import download

def execute(argv,logfile):
    print('RUN', ' '.join(map(str,argv)),flush=True)
    with Path(logfile).open('w') as fh:
        p=subprocess.run(list(map(str,argv)),stdout=fh,stderr=subprocess.STDOUT)
    if p.returncode:print('FAILED:',str(logfile),'exit code',p.returncode,flush=True)
    else:print('COMPLETED:',str(logfile),flush=True)
    return p.returncode==0

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',default=str(Path(__file__).resolve().parent.parent))
    p.add_argument('--chemical-file',help='Official unmodified IMPPAT chemical download, if auto-download fails')
    p.add_argument('--bioactivity-file',help='Official unmodified IMPPAT bioactivity download, if auto-download fails')
    p.add_argument('--no-download',action='store_true',help='Use provided files, no requests')
    p.add_argument('--external-reviewed-file',help='Optional reviewed source-verified human EGFR IC50 dataset CSV (SMILES,IC50_nM,Source_Reference)')
    p.add_argument('--skip-heavy',action='store_true',help='Skip MLP/curation; used for smoke tests only')
    p.add_argument('--scaffold-only',action='store_true',help='Skip random CV when rerunning benchmark')
    args=p.parse_args()
    project=Path(args.project).expanduser().resolve(); code=project/'code';results=project/'results';sources=project/'sources';data=project/'data'
    for d in [results,sources]:d.mkdir(exist_ok=True,parents=True)
    train=data/'CHEMBL203_IC50.csv';assert train.exists(),f'Experimental CHEMBL203 CSV missing: {train}'
    report={'run_started_utc':datetime.now(timezone.utc).isoformat(),'train_sha256':hashlib.sha256(train.read_bytes()).hexdigest(),
      'target':'EGFR CHEMBL203 / UniProt P00533','stages':{},'scientific_limitations':[]}
    python=sys.executable
    cmd=[python,code/'benchmark_natdrug_ai.py','--csv',train,'--out',results/'experiment1','--target','CHEMBL203']
    if args.scaffold_only:cmd.append('--only-scaffold')
    report['stages']['experiment1_scaffold_benchmark']='COMPLETED' if execute(cmd,results/'experiment1.log') else 'FAILED'
    chem=Path(args.chemical_file).expanduser() if args.chemical_file else None
    bio=Path(args.bioactivity_file).expanduser() if args.bioactivity_file else None
    if (not chem or not bio) and not args.no_download:
        try:
            manifest=download(sources)
            if not chem and 'chemical' in manifest['files']:chem=Path(manifest['files']['chemical']['path'])
            if not bio and 'bioactivity' in manifest['files']:bio=Path(manifest['files']['bioactivity']['path'])
            report['stages']['source_download']='PARTIAL' if manifest['warnings'] else 'COMPLETED'
            report['scientific_limitations'].extend(manifest['warnings'])
        except Exception as e:
            report['stages']['source_download']='FAILED';report['scientific_limitations'].append('Official download failure: '+str(e))
    # Include manually downloaded official tables with stable names.
    if not chem:
        files=sorted(sources.glob('imppat_chemical.*'));chem=files[0] if files else None
    if not bio:
        files=sorted(sources.glob('imppat_bioactivity.*'));bio=files[0] if files else None
    if chem and bio and chem.exists() and bio.exists():
        out=results/'imppat_source_audit'
        ok=execute([python,code/'prepare_imppat.py','--chemical',chem,'--bioactivity',bio,'--train',train,'--out',out],results/'imppat_preparation.log')
        report['stages']['experiment2_egfr_measured_candidates']='EXTRACTED_NOT_INDEPENDENCE_VERIFIED' if ok else 'FAILED'
        if ok:
            # An external source library of IMPPAT phytochemicals is a valid in-silico screening set.
            library=out/'natural_product_screening_library.csv'
            if library.exists():
                report['stages']['experiment3_natural_product_screening']='COMPLETED_PREDICTIONS_ONLY' if execute([python,code/'screen_natural_products.py','--train-csv',train,'--library-csv',library,
                   '--library-provenance','IMPPAT 3.0, release 2026-09-30, https://cb.imsc.res.in/imppat/download; see preparation_manifest.json',
                   '--model','RF_fingerprint_descriptors','--out',results/'experiment3_screening'],results/'experiment3.log') else 'FAILED'
            report['scientific_limitations'].append('IMPPAT reuses published/BindingDB/NPASS measurements. Exact source publication and assay independence is not proven by excluding structural matches alone.')
    else:
        report['stages']['experiment2_egfr_measured_candidates']='BLOCKED_SOURCE_FILES_UNAVAILABLE'
        report['stages']['experiment3_natural_product_screening']='BLOCKED_SOURCE_FILES_UNAVAILABLE'
        report['scientific_limitations'].append('IMPPAT official chemical and bioactivity bulk files not available. No artificial natural-product values generated.')
    if args.external_reviewed_file:
        verified=Path(args.external_reviewed_file)
        if verified.exists():
            report['stages']['experiment2_reviewed_external_metrics']='COMPLETED_AUTHOR_REVIEW_REQUIRED' if execute([python,code/'external_validate.py','--train-csv',train,'--test-csv',verified,'--test-provenance',f'Independence of original publication and assay records verified by submitting author: {verified.name}',
              '--model','RF_fingerprint_descriptors','--out',results/'experiment2_reviewed_external'],results/'experiment2.log') else 'FAILED'
        else:report['stages']['experiment2_reviewed_external_metrics']='FAILED_MISSING_REVIEWED_FILE'
    else:
        report['stages']['experiment2_reviewed_external_metrics']='NOT_PERFORMED_NO_SOURCE_VERIFIED_EXTERNAL_SET'
    if args.skip_heavy:
        report['stages']['experiment4_robustness']='SKIPPED_BY_FLAG'
    else:
        report['stages']['experiment4_robustness']='COMPLETED' if execute([python,code/'robustness_studies.py','--csv',train,'--out',results/'experiment4_robustness','--mode','all','--seeds','42,52,62'],results/'experiment4.log') else 'FAILED'
    (project/'STUDY_STATUS.json').write_text(json.dumps(report,indent=2))
    md=['# NatDrug-AI automatic study status','',f"Run: {report['run_started_utc']}",'','## Research stages','']
    md.extend(f'- **{k}:** {v}' for k,v in report['stages'].items())
    md.extend(['','## Qualifications']+['- '+q for q in report['scientific_limitations']])
    (project/'STUDY_STATUS.md').write_text('\n'.join(md)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
