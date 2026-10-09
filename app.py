import io, json, hashlib, time
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import requests
import joblib
import streamlit as st
import plotly.express as px
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED, rdMolDescriptors, rdFingerprintGenerator, Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold, KFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'; MODELS=ROOT/'models'
DATA.mkdir(exist_ok=True); MODELS.mkdir(exist_ok=True)
BASE='https://www.ebi.ac.uk/chembl/api/data'
FP=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
TARGET_PRESETS={'Oncology':['EGFR','ERBB2','BRAF','MAP2K1','PARP1','CDK4','CDK6','PIK3CA','ALK','MET','KIT','FLT3','JAK2','HDAC1','ESR1','AR','KRAS','BTK','VEGFR2'], 'Metabolic':['DPP4','PPARG','SGLT2','GLP1R','HMGCR'], 'Neurology':['ACHE','MAOB','DRD2','SLC6A4','BACE1'], 'Inflammation':['PTGS2','ALOX5','JAK1','TNF'], 'Infectious disease':['HIV1 integrase','HIV1 protease','InhA','SARS-CoV-2 main protease'], 'Cardiovascular':['ACE','AGTR1','ADRB1','F2','F10']}
st.set_page_config(page_title='NatDrug-AI Research',page_icon='🧬',layout='wide')
st.markdown('''<style>.stApp{background:linear-gradient(135deg,#091526,#11243b);color:#e8f3fc}.block-container{padding-top:1.5rem}.stButton>button[kind="primary"]{background:#0fae9e;color:white;border-radius:12px}.stMetric{background:#182c44;border-radius:12px;padding:14px}h1,h2,h3{color:#70e2d1!important}div[data-testid="stTabs"] button{font-weight:650}</style>''',unsafe_allow_html=True)
st.title('🧬 NatDrug-AI | Research Studio')
st.caption('Live ChEMBL evidence • RDKit • Random Forest • PyTorch MLP • scaffold validation • natural-product screening')
st.warning('Research prototype: neither ChEMBL target membership nor predicted pIC50 establishes clinical efficacy. The model is trained on experimental IC50, not on approved-drug labels.')

def get_json(endpoint,params=None):
    r=requests.get(f'{BASE}/{endpoint}',params=params,timeout=40,headers={'Accept':'application/json','User-Agent':'NatDrugAI-research/1.0'})
    r.raise_for_status(); return r.json()

def pages(endpoint,params=None,max_records=10000,key=None):
    params=dict(params or {}); params['limit']=min(1000,max_records);params['offset']=0
    out=[]
    while len(out)<max_records:
        payload=get_json(endpoint,params)
        chunk=payload.get(key or endpoint.split('/')[0],[])
        if not chunk:break
        out.extend(chunk[:max_records-len(out)])
        if not payload.get('page_meta',{}).get('next'):break
        params['offset']+=len(chunk)
    return out

@st.cache_data(ttl=86400,show_spinner=False)
def clinical_catalog(max_records=15000):
    # Evidence definition: approved-drug mechanism records, not all clinical-trial targets.
    mechanisms=pages('mechanism.json',{'max_phase':4},max_records,'mechanisms')
    rows=[]
    for m in mechanisms:
        tid=m.get('target_chembl_id')
        if tid and tid!='CHEMBL0':
            rows.append({'Target_ID':tid,'Drug_ID':m.get('molecule_chembl_id'),'Mechanism':m.get('mechanism_of_action'),'Action':m.get('action_type'),'Evidence':'ChEMBL approved-drug mechanism'})
    return pd.DataFrame(rows).drop_duplicates() if rows else pd.DataFrame(columns=['Target_ID','Drug_ID','Mechanism','Action','Evidence'])

@st.cache_data(ttl=86400,show_spinner=False)
def search_targets(term):
    d=get_json('target/search.json',{'q':term,'limit':100})
    return [{'ID':x['target_chembl_id'],'Name':x.get('pref_name',''),'Organism':x.get('organism',''),'Type':x.get('target_type','')} for x in d.get('targets',[])]

@st.cache_data(ttl=86400,show_spinner=False)
def target_details(tid):
    return get_json(f'target/{tid}.json')

@st.cache_data(ttl=86400,show_spinner=False)
def fetch_activities(tid,max_records):
    # Exact nM IC50, high target confidence, binding assays.
    rows=pages('activity.json',{'target_chembl_id':tid,'standard_type':'IC50','standard_units':'nM','standard_relation':'=','assay_type':'B','confidence_score__gte':8},max_records,'activities')
    out=[]
    for x in rows:
        try:
            v=float(x.get('standard_value'))
            if not np.isfinite(v) or v<=0 or not x.get('canonical_smiles'):continue
            out.append({'SMILES':x['canonical_smiles'],'IC50_nM':v,'Assay_ID':x.get('assay_chembl_id'),'Molecule_ID':x.get('molecule_chembl_id'),'Target_ID':tid})
        except (TypeError,ValueError):pass
    return pd.DataFrame(out,columns=['SMILES','IC50_nM','Assay_ID','Molecule_ID','Target_ID'])

def mol(s):
    if not isinstance(s,str) or not s.strip():return None
    try:
        m=Chem.MolFromSmiles(s.strip());return m if m is not None and m.GetNumHeavyAtoms()>0 else None
    except Exception:return None

def canon(s):
    m=mol(s);return Chem.MolToSmiles(m) if m else None

def features(s):
    m=mol(s)
    if m is None:return None
    a=np.zeros((2048,),dtype=np.float32)
    DataStructs.ConvertToNumpyArray(FP.GetFingerprint(m),a)
    return a

def properties(s):
    m=mol(s)
    if m is None:raise ValueError('Invalid SMILES')
    mw=Descriptors.MolWt(m);logp=Crippen.MolLogP(m);hbd=Lipinski.NumHDonors(m);hba=Lipinski.NumHAcceptors(m);tpsa=rdMolDescriptors.CalcTPSA(m);rot=Lipinski.NumRotatableBonds(m)
    return {'SMILES':Chem.MolToSmiles(m),'MW':round(mw,2),'LogP':round(logp,2),'HBD':hbd,'HBA':hba,'TPSA':round(tpsa,2),'Rotatable_Bonds':rot,'QED':round(QED.qed(m),3),'Lipinski_Violations':sum([mw>500,logp>5,hbd>5,hba>10]),'Veber_Pass':bool(rot<=10 and tpsa<=140)}

def read_molecules(file):
    raw=file.getvalue();name=file.name.lower()
    if name.endswith('.sdf'):
        supplier=Chem.ForwardSDMolSupplier(io.BytesIO(raw),removeHs=False)
        return pd.DataFrame([{'Name':m.GetProp('_Name') if m.HasProp('_Name') else f'Compound_{i}','SMILES':Chem.MolToSmiles(Chem.RemoveHs(m))} for i,m in enumerate(supplier) if m is not None])
    df=pd.read_csv(io.BytesIO(raw))
    if 'SMILES' not in df.columns:raise ValueError('CSV requires SMILES column')
    return df

def curate(df):
    if not {'SMILES','IC50_nM'}.issubset(df.columns):raise ValueError('Training CSV requires SMILES and IC50_nM')
    d=df.copy();d['IC50_nM']=pd.to_numeric(d['IC50_nM'],errors='coerce');d=d[np.isfinite(d['IC50_nM'])&(d['IC50_nM']>0)].copy()
    d['SMILES']=d.SMILES.map(canon);d=d.dropna(subset=['SMILES'])
    d['pIC50']=9-np.log10(d.IC50_nM)
    d=d.groupby('SMILES',as_index=False).agg(pIC50=('pIC50','median'),Replicates=('pIC50','size'),Assay_SD=('pIC50','std'))
    d['Scaffold']=d.SMILES.map(lambda s:MurckoScaffold.MurckoScaffoldSmiles(smiles=s))
    return d

class TorchMLP:
    def __init__(self,epochs=120):self.epochs=epochs
    def fit(self,X,y):
        import torch
        torch.manual_seed(42)
        self.mean=float(np.mean(y));self.sd=max(float(np.std(y)),1e-6)
        xx=torch.tensor(X,dtype=torch.float32);yy=torch.tensor((y-self.mean)/self.sd,dtype=torch.float32).reshape(-1,1)
        self.net=torch.nn.Sequential(torch.nn.Linear(X.shape[1],256),torch.nn.ReLU(),torch.nn.Dropout(0.2),torch.nn.Linear(256,64),torch.nn.ReLU(),torch.nn.Linear(64,1))
        opt=torch.optim.AdamW(self.net.parameters(),lr=0.001,weight_decay=0.01)
        self.net.train()
        for _ in range(self.epochs):
            opt.zero_grad();loss=torch.nn.functional.mse_loss(self.net(xx),yy);loss.backward();opt.step()
        self.net.eval();return self
    def predict(self,X):
        import torch
        with torch.no_grad():return self.net(torch.tensor(X,dtype=torch.float32)).numpy().ravel()*self.sd+self.mean
    def state(self):return {k:v.detach().cpu().numpy() for k,v in self.net.state_dict().items()}
    def load_state(self,state):
        import torch
        self.net=torch.nn.Sequential(torch.nn.Linear(2048,256),torch.nn.ReLU(),torch.nn.Dropout(0.2),torch.nn.Linear(256,64),torch.nn.ReLU(),torch.nn.Linear(64,1))
        self.net.load_state_dict({k:torch.tensor(v) for k,v in state.items()});self.net.eval()

def make_model(kind):return RandomForestRegressor(n_estimators=250,min_samples_leaf=2,max_features='sqrt',random_state=42,n_jobs=-1) if kind=='Random Forest' else TorchMLP()

def train(d,kind):
    clean=curate(d);n=len(clean)
    if n<20:raise ValueError(f'{n} unique compounds: minimum 20 required for exploratory training')
    X=np.stack([features(s) for s in clean.SMILES]);y=clean.pIC50.to_numpy();groups=clean.Scaffold.to_numpy();unique=len(set(groups))
    if unique>=3:
        splits=list(GroupKFold(n_splits=min(5,unique)).split(X,y,groups));method='Scaffold GroupKFold'
    else:
        splits=list(KFold(n_splits=min(5,n),shuffle=True,random_state=42).split(X,y));method='Random KFold: insufficient distinct scaffolds'
    pred=np.full(n,np.nan)
    for tr,te in splits:
        m=make_model(kind);m.fit(X[tr],y[tr]);pred[te]=m.predict(X[te])
    metrics={'MAE':float(mean_absolute_error(y,pred)),'RMSE':float(np.sqrt(mean_squared_error(y,pred))),'R2':float(r2_score(y,pred)),'N':n,'Scaffolds':unique,'Validation':method,'Model':kind,'Exploratory':n<50}
    final=make_model(kind);final.fit(X,y)
    comparison=clean[['SMILES','pIC50','Scaffold']].copy();comparison['Predicted_pIC50']=pred
    return final,clean,metrics,comparison

def save_model(model,clean,metrics,target):
    token=''.join(c if c.isalnum() or c in '_-' else '_' for c in target)[:64]
    stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    key=f'{token}_{metrics["Model"].replace(" ","_")}_{stamp}'
    payload={'metrics':metrics,'target':target,'smiles':clean.SMILES.tolist(),'created':stamp,'kind':metrics['Model'],'fingerprint':{'radius':2,'bits':2048}}
    if metrics['Model']=='Random Forest':joblib.dump(model,MODELS/f'{key}.joblib')
    else:
        import torch
        torch.save({'state':{k:torch.tensor(v) for k,v in model.state().items()},'mean':model.mean,'sd':model.sd},MODELS/f'{key}.pt')
    (MODELS/f'{key}.json').write_text(json.dumps(payload,indent=2));return key

def load_model(key):
    meta=json.loads((MODELS/f'{key}.json').read_text())
    if meta['kind']=='Random Forest':model=joblib.load(MODELS/f'{key}.joblib')
    else:
        import torch
        obj=torch.load(MODELS/f'{key}.pt',map_location='cpu',weights_only=True)
        model=TorchMLP();model.mean=obj['mean'];model.sd=obj['sd'];model.load_state(obj['state'])
    return model,meta

def predict(s,model,meta):
    p=properties(s);x=features(s).reshape(1,-1);score=float(model.predict(x)[0]);q=FP.GetFingerprint(mol(s));sims=[DataStructs.TanimotoSimilarity(q,FP.GetFingerprint(mol(t))) for t in meta['smiles']]
    p.update({'Target':meta['target'],'Predicted_pIC50':round(score,3),'Estimated_IC50_nM':round(float(10**np.clip(9-score,-15,15)),3),'Max_Tanimoto':round(max(sims),3),'Model':meta['kind'],'Exploratory_Model':meta['metrics']['Exploratory']})
    return p

with st.sidebar:
    st.header('Research workspace')
    st.write('Local models:',len(list(MODELS.glob('*.json'))))
    st.caption('Approved-target catalog is assembled from ChEMBL mechanism records. Verify drug approval and target assignment individually.')
    st.markdown('[ChEMBL database](https://www.ebi.ac.uk/chembl/)')

tabs=st.tabs(['🎯 Clinical targets','🧪 Train AI','🔬 Compound explorer','📦 Batch screening','📈 Model reports'])
with tabs[0]:
    st.subheader('Evidence-linked drug-target discovery')
    c1,c2=st.columns(2)
    with c1:
        area=st.selectbox('Therapeutic area',list(TARGET_PRESETS));preset=st.selectbox('Example target search',TARGET_PRESETS[area]);term=st.text_input('Search ChEMBL target',value=preset)
        if st.button('Search targets',type='primary'):
            try:st.session_state['hits']=search_targets(term)
            except Exception as e:st.error(f'ChEMBL search failed: {e}')
        hits=st.session_state.get('hits',[])
        if hits:st.dataframe(pd.DataFrame(hits),use_container_width=True)
    with c2:
        if st.button('Build approved-drug mechanism catalog'):
            try:
                with st.spinner('Retrieving approved-drug mechanism records...'):
                    cat=clinical_catalog();st.session_state['catalog']=cat
                st.success(f'{cat.Target_ID.nunique()} distinct target identifiers in retrieved records')
            except Exception as e:st.error(f'Catalog retrieval failed: {e}')
        cat=st.session_state.get('catalog')
        if cat is not None:
            st.dataframe(cat,use_container_width=True,height=330)
            st.download_button('Download evidence catalog CSV',cat.to_csv(index=False),'approved_mechanism_targets.csv')
    st.info('Clinical target status: ChEMBL approved-drug mechanism associations are evidence for approved-drug targeting, not a comprehensive registry of all trial-stage targets. Complexes and nonprotein targets may be included. The catalog is fetched live, not hard-coded.')
with tabs[1]:
    st.subheader('Experimental bioactivity → reproducible training')
    mode=st.radio('Dataset source',['Live ChEMBL IC50','Upload measured IC50 CSV'],horizontal=True)
    if mode=='Live ChEMBL IC50':
        tid=st.text_input('Exact ChEMBL target ID (from target search)',value='CHEMBL203')
        limit=st.number_input('Maximum records',min_value=20,max_value=50000,value=5000,step=500)
        if st.button('Fetch experimental IC50 records',type='primary'):
            try:
                with st.spinner('Retrieving ChEMBL records...'):df=fetch_activities(tid.strip(),int(limit))
                st.session_state['train_data']=df;st.session_state['train_target']=tid.strip()
                if len(df):df.to_csv(DATA/f'{tid.strip()}_IC50.csv',index=False)
                st.success(f'{len(df)} exact IC50 measurements retrieved')
            except Exception as e:st.error(f'Download error: {e}')
    else:
        up=st.file_uploader('Upload measured IC50 CSV (SMILES, IC50_nM)',type=['csv'],key='training')
        label=st.text_input('Dataset target identifier',value='Custom_Target')
        if up is not None:
            st.session_state['train_data']=pd.read_csv(up);st.session_state['train_target']=label
    df=st.session_state.get('train_data')
    if df is not None:
        st.write('Measurements:',len(df));st.dataframe(df.head(30),use_container_width=True)
        try:
            curated=curate(df);st.metric('Unique compounds after curation',len(curated));st.metric('Unique scaffolds',curated.Scaffold.nunique())
            st.download_button('Download curated training CSV',curated.to_csv(index=False),'curated_training.csv')
        except Exception as e:st.error(f'Curation error: {e}')
        kind=st.selectbox('Model architecture',['Random Forest','PyTorch MLP (deep learning)'])
        if st.button('Train, cross-validate and save',type='primary'):
            try:
                with st.spinner('Cross-validating and training; deep learning may take several minutes...'):
                    m,clean,metrics,comparison=train(df,kind);key=save_model(m,clean,metrics,st.session_state.get('train_target','Unknown'))
                st.session_state['latest_metrics']=metrics;st.session_state['latest_cv']=comparison
                st.success(f'Saved model: {key}')
            except Exception as e:st.exception(e)
    metrics=st.session_state.get('latest_metrics')
    if metrics:
        a,b,c=st.columns(3);a.metric('CV MAE',f'{metrics["MAE"]:.3f}');b.metric('CV RMSE',f'{metrics["RMSE"]:.3f}');c.metric('CV R²',f'{metrics["R2"]:.3f}')
        st.caption(f'{metrics["Validation"]} • {metrics["N"]} compounds • {metrics["Scaffolds"]} scaffolds')
        if metrics['Exploratory']:st.warning('Small dataset: exploratory model, not publication-ready.')
        comp=st.session_state['latest_cv'];fig=px.scatter(comp,x='pIC50',y='Predicted_pIC50',title='Cross-validated actual vs predicted pIC50',template='plotly_dark');st.plotly_chart(fig,use_container_width=True)
        st.download_button('Download out-of-fold predictions',comp.to_csv(index=False),'cv_predictions.csv')

model_keys=sorted([p.stem for p in MODELS.glob('*.json')],reverse=True)
selected=st.sidebar.selectbox('Saved prediction model',model_keys) if model_keys else None
loaded=None
if selected:
    try:loaded=load_model(selected)
    except Exception as e:st.sidebar.error(f'Model load error: {e}')
with tabs[2]:
    st.subheader('Natural compound explorer')
    s=st.text_input('Enter SMILES',value='Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1')
    if st.button('Analyze molecule',type='primary'):
        try:
            p=properties(s);a,b=st.columns([1,2]);a.image(Draw.MolToImage(mol(s),size=(400,330)));b.dataframe(pd.DataFrame([p]),use_container_width=True)
            if loaded:
                m,meta=loaded;pred=predict(s,m,meta);st.metric('Predicted pIC50',pred['Predicted_pIC50']);st.metric('Nearest training fingerprint similarity',pred['Max_Tanimoto']);st.dataframe(pd.DataFrame([pred]),use_container_width=True)
                if pred['Max_Tanimoto']<0.4:st.warning('Low structural similarity to training compounds; prediction may extrapolate.')
            else:st.info('Train or select a saved model for target-specific predictions.')
        except Exception as e:st.error(str(e))
with tabs[3]:
    st.subheader('Upload your own compound library')
    st.caption('CSV: Name,SMILES • SDF: molecular structures. Invalid structures are reported rather than silently removed.')
    up=st.file_uploader('Choose CSV or SDF',type=['csv','sdf'],key='screen')
    if up is not None and st.button('Screen uploaded compounds',type='primary'):
        try:
            d=read_molecules(up);rows=[]
            for i,row in d.iterrows():
                try:
                    p=predict(row.SMILES,*loaded) if loaded else properties(row.SMILES)
                    p['Status']='Valid'
                except Exception as e:p={'SMILES':str(row.SMILES),'Status':str(e)}
                p['Name']=row.get('Name',f'Molecule_{i+1}');rows.append(p)
            result=pd.DataFrame(rows);st.dataframe(result,use_container_width=True)
            st.download_button('Download screening results',result.to_csv(index=False),'NatDrug_AI_predictions.csv')
        except Exception as e:st.error(f'Upload failed: {e}')
with tabs[4]:
    st.subheader('Reproducibility and model provenance')
    if loaded:
        _,meta=loaded;st.json({k:v for k,v in meta.items() if k!='smiles'})
        st.caption('Model persistence is local. Only load models you trained or trust. Results require external experimental validation.')
    else:st.info('No trained model selected.')
    st.markdown('**Publication checklist:** external scaffold holdout; target and assay harmonization; natural-product external set; matched baselines; uncertainty calibration; reproducible versions; ablation studies; experimental confirmation.')
