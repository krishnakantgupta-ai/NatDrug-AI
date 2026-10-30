"""Source-preserving human-EGFR IC50 preparation plus independent NP screening input.
Outputs *candidate* external data. Publication-grade assay/source independence requires manual verification.
"""
import argparse, io,re,zipfile,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from rdkit import Chem, RDLogger
from benchmark_natdrug_ai import read_curate,standardize
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog('rdApp.warning')

def read_any(file):
    file=Path(file);blob=file.read_bytes();name=file.name
    if file.suffix.lower()=='.zip':
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            names=[n for n in z.namelist() if n.lower().endswith(('.csv','.tsv','.txt','.xlsx','.xls')) and '__MACOSX' not in n]
            if len(names)!=1:raise ValueError(f'ZIP must contain one tabular file; found {names[:6]}')
            name=names[0];blob=z.read(name)
    suffix=Path(name).suffix.lower()
    if suffix in ('.xlsx','.xls'):return pd.read_excel(io.BytesIO(blob),dtype=str)
    for e in ['utf-8-sig','latin-1']:
        try:return pd.read_csv(io.BytesIO(blob),sep=None,engine='python',dtype=str,encoding=e,on_bad_lines='warn')
        except UnicodeDecodeError:continue
    raise ValueError('Cannot read table '+str(file))

def pick(df,names,optional=False):
    normalize=lambda a:re.sub('[^a-z0-9]','',str(a).lower())
    lookup={normalize(c):c for c in df.columns}
    for x in names:
        if normalize(x) in lookup:return lookup[normalize(x)]
    if optional:return None
    raise ValueError(f'Cannot match required column amongst {list(df.columns)}; expected one of {names}')

def create(chemfile,biofile,trainfile,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    chem=read_any(chemfile);bio=read_any(biofile)
    ids=['Source_ID','IMPPAT_Phytochemical_identifier','IMPPAT Phytochemical identifier','IMPPAT phytochemical ID','Phytochemical identifier','Phytochemical ID','IMPPAT ID','IMPPAT3 PHYID','IMPPAT3_PHYID','Phytochemical_ID','PHYID','Compound ID','ID']
    ci=pick(chem,ids);cs=pick(chem,['SMILES','Canonical SMILES','Canonical_SMILES','Isomeric SMILES','Standardized SMILES','SMILES String'])
    cn=pick(chem,['Phytochemical name_standardised','Phytochemical Name Standardised','Phytochemical Name','Compound Name','Chemical Name','Name'],optional=True)
    bi=pick(bio,ids); bt=pick(bio,['Target Name','Protein Name','Target'],optional=True)
    bu=pick(bio,['UniProt identifier','Uniprot ID','Uniprot_identifier','UniProt','UniProt ID','Target Uniprot identifier','UniProt accession'],optional=True)
    bo=pick(bio,['Target organism','Target Organism','Organism','Species'],optional=True)
    bv=pick(bio,['IC50','IC50 Value','IC50_values','IC50 (nM)','IC50_nM','IC 50'])
    bs=pick(bio,['Source','Data source','Activity Source','Activity database'],optional=True)
    bd=pick(bio,['DOI','Reference DOI','Publication DOI','Reference','PMID','PubMed ID','Article DOI'],optional=True)
    if not (bt or bu):raise ValueError('No experimentally target-identifying column')
    valid=chem[[ci,cs]+([cn] if cn else [])].copy();valid.columns=['Source_ID','SMILES']+(['Name'] if cn else [])
    if cn is None:valid['Name']=valid.Source_ID
    valid=valid.dropna(subset=['SMILES']).copy();valid['Source_ID']=valid.Source_ID.astype(str).str.strip()
    valid['SMILES']=valid.SMILES.astype(str).str.strip();valid=valid[valid.SMILES.ne('')].drop_duplicates('Source_ID')
    # Explicit database tag is essential; chemical structures alone cannot prove origin.
    valid['Source_Database']='IMPPAT 3.0 (2026-09-30)'
    valid.to_csv(out/'imppat_all_phytochemicals_raw.csv',index=False)
    matches=[];count={'bio_rows':len(bio),'human_egfr_rows':0,'candidate_activity_tokens':0,'exact_nM_values':0,'censored_or_non_numeric':0,'non_molar_units':0}
    regex=re.compile(r'^\s*(?P<relation><=|>=|<|>|=|≤|≥)?\s*(?P<value>(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(?P<unit>pM|nM|uM|µM|μM|mM|M)?\s*$',re.I)
    factors={'pm':0.001,'nm':1.,'um':1000.,'mm':1e6,'m':1e9}
    explicit_nm_in_header=bool(re.search(r'nm\b',bv,re.I))
    for _,row in bio.iterrows():
        uni=str(row[bu]).strip().upper() if bu else ''
        target=str(row[bt]).strip().lower() if bt else ''
        if uni and uni not in ('NAN','NONE',''):
            match=bool(re.search(r'(?<![A-Z0-9])P00533(?![A-Z0-9])',uni))
        else:match=target in {'egfr','epidermal growth factor receptor','epidermal growth factor receptor erbb1','epidermal growth factor receptor erb-b1'}
        if not match:continue
        org=str(row[bo]).lower() if bo else ''
        if bo and not ('homo sapiens' in org or 'human' in org):continue
        count['human_egfr_rows']+=1
        for token in re.split('[|;]',str(row[bv])):
            token=token.strip()
            if not token or token.lower() in ('nan','none'):continue
            count['candidate_activity_tokens']+=1
            m=regex.fullmatch(token)
            if m is None:
                count['censored_or_non_numeric']+=1;continue
            if m.group('relation') in ['<','>','<=','>=','≤','≥']:
                count['censored_or_non_numeric']+=1;continue
            unit=m.group('unit')
            if unit is None and not explicit_nm_in_header:
                count['non_molar_units']+=1;continue
            unit=(unit or 'nM').lower().replace('µ','u').replace('μ','u')
            value=float(m.group('value'))*factors[unit]
            if not np.isfinite(value) or value<=0:continue
            matches.append({'Source_ID':str(row[bi]).strip(),'IC50_nM':value,
              'Raw_IC50':token,'Original_Source':str(row[bs]) if bs else 'not stated',
              'Original_Reference':str(row[bd]) if bd else 'not stated',
              'Target_UniProt':'P00533','Target_Organism':'Homo sapiens'})
            count['exact_nM_values']+=1
    train,traincnt=read_curate(trainfile,'parent','pooled');trains=set(train.Canonical_SMILES)
    source=valid.copy();source['Canonical_SMILES']=source.SMILES.map(lambda x:standardize(x,'parent'))
    source=source.dropna(subset=['Canonical_SMILES']).drop_duplicates('Canonical_SMILES')
    source['Training_Compound_Overlap']=source.Canonical_SMILES.isin(trains)
    screened=source.loc[~source.Training_Compound_Overlap].copy()
    # Retain source original name and IDs to enable external audit and plant-origin traceability.
    screened[['Name','SMILES','Source_ID','Source_Database']].to_csv(out/'natural_product_screening_library.csv',index=False)
    report={'chemical_file':str(chemfile),'bioactivity_file':str(biofile),'chemical_sha256':hashlib.sha256(Path(chemfile).read_bytes()).hexdigest(),
       'bioactivity_sha256':hashlib.sha256(Path(biofile).read_bytes()).hexdigest(),
       'chemical_records':len(chem),'unique_valid_phytochemicals':len(source),'screening_after_overlap':len(screened),
       'training_compound_overlap_excluded':int(source.Training_Compound_Overlap.sum()),'egfr_curation':count,'detected_columns':{'chemical_id':ci,'chemical_smiles':cs,'chemical_name':cn,'bioactivity_id':bi,'bioactivity_target':bt,'bioactivity_uniprot':bu,'bioactivity_organism':bo,'bioactivity_ic50':bv,'bioactivity_source':bs,'bioactivity_reference':bd},
       'warning':'Only structure-disjoint candidates and in-silico predictions are produced. Independent EGFR IC50 evaluation remains provisional until original publications and assay provenance are checked.'}
    if matches:
        d=pd.DataFrame(matches).merge(source[['Source_ID','Name','SMILES','Canonical_SMILES']],on='Source_ID',how='left')
        d['Training_Compound_Overlap']=d.Canonical_SMILES.isin(trains)
        d.to_csv(out/'all_egfr_measurements_overlap_audit.csv',index=False)
        cand=d.loc[~d.Training_Compound_Overlap & d.Canonical_SMILES.notna()].copy()
        cand['Scaffold']=cand.Canonical_SMILES.map(lambda x:MurckoScaffold.MurckoScaffoldSmiles(smiles=x))
        cand['Source_Review_Status']='UNVERIFIED_original_assay_or_publication_overlap'
        cand.to_csv(out/'external_egfr_candidates_NOT_YET_VALIDATED.csv',index=False)
        report['candidate_measured_egfr_records']=len(cand)
        report['candidate_measured_egfr_unique_compounds']=cand.Canonical_SMILES.nunique()
    else:report['candidate_measured_egfr_records']=0
    (out/'preparation_manifest.json').write_text(json.dumps(report,indent=2))
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--chemical',required=True);p.add_argument('--bioactivity',required=True);p.add_argument('--train',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();print(json.dumps(create(a.chemical,a.bioactivity,a.train,a.out),indent=2))
