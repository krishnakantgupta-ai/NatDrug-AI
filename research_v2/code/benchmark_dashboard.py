"""NatDrug-AI Research Evidence Dashboard: view *real* exported benchmark outputs.

Use `streamlit run code/benchmark_dashboard.py` after running benchmark scripts.
The dashboard never trains a model or makes unsupported biological claims.
"""
from pathlib import Path
import pandas as pd
import streamlit as st
import plotly.express as px

st.set_page_config(page_title='NatDrug-AI | Validation Evidence',page_icon='🧬',layout='wide')
st.title('🧬 NatDrug-AI | Validation Evidence')
st.caption('Scaffold-aware out-of-fold results • error versus similarity • external measured validation • independent natural-product screening')
st.warning('The dashboard shows exported computations, not experimental confirmation. An unlabeled library is not external predictive validation.')
root=Path(st.sidebar.text_input('Results directory (local path)',value='results/MolInf_benchmark')).expanduser()
st.sidebar.caption('Run the benchmark CLI first. This view reads files only and cannot fabricate missing results.')

@st.cache_data(show_spinner=False)
def load_csv(path):return pd.read_csv(path)

def needed(path):
 if not path.exists():
  st.error(f'Missing: {path}. Run the associated CLI script first.');return None
 return load_csv(str(path.resolve()))

summary=needed(root/'summary_metrics.csv')
if summary is not None:
 st.subheader('Fixed-comparator cross-validation')
 st.dataframe(summary,use_container_width=True,hide_index=True)
 st.caption('Scaffold CV estimates unseen-framework generalization; random CV is an easier comparison and must be labeled as such.')
 if {'Split','Model','MAE_pIC50'}.issubset(summary):
  st.plotly_chart(px.bar(summary,x='Model',y='MAE_pIC50',color='Split',barmode='group',title='Out-of-fold MAE by model and split'),use_container_width=True)

oof=needed(root/'oof_predictions.csv')
if oof is not None:
 st.subheader('Compound-level out-of-fold evidence')
 a,b=st.columns(2)
 split=a.selectbox('Validation split',sorted(oof.Split.unique()))
 model=b.selectbox('Model',sorted(oof.Model.unique()))
 part=oof[(oof.Split==split)&(oof.Model==model)].copy()
 part['Absolute_error_pIC50']=(part.Observed_pIC50-part.Predicted_pIC50).abs()
 a,b=st.columns(2)
 a.plotly_chart(px.scatter(part,x='Observed_pIC50',y='Predicted_pIC50',hover_data=['SMILES','Fold'],title='Held-out measured versus predicted'),use_container_width=True)
 b.plotly_chart(px.scatter(part,x='Max_train_Tanimoto',y='Absolute_error_pIC50',hover_data=['SMILES','Fold'],title='Held-out error versus train-only similarity'),use_container_width=True)
 st.download_button('Download selected out-of-fold records',part.to_csv(index=False).encode(),'selected_oof.csv','text/csv')

curve=root/'mlp_learning_curves.csv'
if curve.exists():
 st.subheader('MLP training stability')
 x=load_csv(str(curve.resolve()))
 sel=x[(x.Split=='scaffold')&(x.Seed==x.Seed.unique()[0])&(x.Fold==x.Fold.unique()[0])]
 st.plotly_chart(px.line(sel,x='epoch',y=['train_loss','val_loss'],title='Representative scaffold-held-out early-stopping learning curve'),use_container_width=True)

st.subheader('Measured external-test results')
ext_dir=Path(st.text_input('External result directory',value='results/external')).expanduser()
if (ext_dir/'external_predictions.csv').exists():
 ext=load_csv(str((ext_dir/'external_predictions.csv').resolve()))
 st.dataframe(ext[['Canonical_SMILES','pIC50','Prediction_pIC50','Max_train_Tanimoto','Scaffold_in_training']].head(40),hide_index=True)
 st.plotly_chart(px.scatter(ext,x='pIC50',y='Prediction_pIC50',color='Scaffold_in_training',title='Measured independent external test'),use_container_width=True)
else:st.info('No independently measured external dataset output has been loaded.')

st.subheader('Natural-product library deployment')
npdir=Path(st.text_input('Natural-product screening result directory',value='results/NP_screening')).expanduser()
if (npdir/'natural_products_ranked_all.csv').exists():
 np=load_csv(str((npdir/'natural_products_ranked_all.csv').resolve()))
 st.metric('Independent valid structures screened',len(np))
 st.plotly_chart(px.histogram(np,x='Predicted_pIC50',nbins=40,title='Unvalidated prediction distribution'),use_container_width=True)
 st.plotly_chart(px.histogram(np,x='Max_Train_Tanimoto',nbins=20,title='Training-set similarity distribution'),use_container_width=True)
 st.dataframe(np.head(30),hide_index=True,use_container_width=True)
 st.download_button('Download ranked natural products',np.to_csv(index=False).encode(),'NP_screening_ranked.csv','text/csv')
else:st.info('No external natural-product screening dataset output has been loaded.')
