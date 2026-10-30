"""Download official IMPPAT 3.0 chemical and bioactivity exports with provenance.
Fails closed if download links are unavailable; never substitutes invented records.
"""
import argparse, hashlib,json,re,sys,urllib.parse
from pathlib import Path
import requests
from bs4 import BeautifulSoup
URL='https://cb.imsc.res.in/imppat/download'
LABELS={'chemical':'Chemical information on IMPPAT phytochemicals','bioactivity':'Bioactivity information for IMPPAT Phytochemicals'}

def get_links(html):
    soup=BeautifulSoup(html,'html.parser')
    links={}
    for typ,label in LABELS.items():
        rows=[tr for tr in soup.select('tr') if label.lower() in tr.get_text(' ',strip=True).lower()]
        if not rows:continue
        row=rows[0]; candidates=[]
        for a in row.select('a, button, [onclick], [data-url]'):
            for attr in ['href','data-url','data-href','onclick','download']:
                v=a.get(attr,'')
                if not v: continue
                if attr=='onclick':
                    matches=re.findall(r'''["']([^"']+\.(?:csv|tsv|zip|xlsx|xls|txt|sdf)(?:\?[^"']*)?)["']''',v,re.I)
                    candidates.extend(matches)
                elif re.search(r'\.(?:csv|tsv|zip|xlsx|xls|txt|sdf)(?:[?#]|$)|/download/|/static/|/media/|/data/',v,re.I):
                    candidates.append(v)
        for cand in candidates:
            u=urllib.parse.urljoin(URL,cand)
            if u.startswith('https://cb.imsc.res.in/') or u.startswith('https://www.cb.imsc.res.in/'):
                links[typ]=u;break
    return links

def download(out,timeout=45):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    headers={'User-Agent':'NatDrugAI-academic-reproducibility/1.0'}
    sess=requests.Session();resp=sess.get(URL,timeout=timeout,headers=headers);resp.raise_for_status()
    (out/'imppat_download_page.html').write_text(resp.text,encoding='utf-8')
    links=get_links(resp.text)
    info={'source':URL,'retrieval_utc':__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),'extracted_links':links,'files':{},'warnings':[]}
    for typ in LABELS:
        if typ not in links:
            info['warnings'].append(f'No direct file link found for {typ}; download manually from official page and place in sources/')
            continue
        url=links[typ]
        try:
            r=sess.get(url,timeout=timeout,headers=headers,stream=True);r.raise_for_status()
            content_type=r.headers.get('content-type','')
            ext=Path(urllib.parse.urlparse(url).path).suffix.lower()
            if ext not in ['.csv','.tsv','.zip','.xlsx','.xls','.txt','.sdf']:ext='.dat'
            p=out/(f'imppat_{typ}'+ext)
            digest=hashlib.sha256();n=0
            with p.open('wb') as fh:
                for part in r.iter_content(1024*512):
                    if part:fh.write(part);digest.update(part);n+=len(part)
            if n<1000 or ('text/html' in content_type.lower()):
                p.unlink(missing_ok=True);raise ValueError(f'File was empty or HTML: {content_type} {n} bytes')
            info['files'][typ]={'path':str(p),'url':url,'bytes':n,'sha256':digest.hexdigest(),'content_type':content_type}
        except Exception as e:info['warnings'].append(f'Could not automatically download {typ}: {e}')
    (out/'download_manifest.json').write_text(json.dumps(info,indent=2))
    return info

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',default='sources');a=p.parse_args(); print(json.dumps(download(a.out),indent=2))
