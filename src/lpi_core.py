from __future__ import annotations
import csv, re, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

AXIS = ['GHR','IGF1','IGF1R','IRS1','IRS2','PIK3CA','PIK3CB','PDPK1','AKT1','AKT2','MTOR','RPTOR','RPS6KB1','EIF4EBP1']
FOXO = ['BCL2L11','CDKN1B','GADD45A','SESN3','SOD2','CAT']
PROLIF = ['MKI67','PCNA','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','CDK1','CCNB1','TOP2A','UBE2C','CDC20']
SOLID_LINEAGES = ['Lung','CNS/Brain','Skin','Esophagus/Stomach','Bowel','Head and Neck','Ovary/Fallopian Tube','Breast','Pancreas','Soft Tissue','Biliary Tract','Bladder/Urinary Tract','Peripheral Nervous System','Uterus','Bone','Kidney','Liver','Pleura','Cervix','Eye']

def symbol(label: str) -> str:
    return re.sub(r' \(\d+\)$','',str(label).strip('"'))

def bh_adjust(p):
    p=np.asarray(p,float)
    n=len(p)
    order=np.argsort(p)
    ranked=p[order]*n/np.arange(1,n+1)
    ranked=np.minimum.accumulate(ranked[::-1])[::-1]
    q=np.empty(n,float); q[order]=np.minimum(ranked,1.0)
    return q

def zscore_cols(a):
    a=np.asarray(a,float)
    mu=np.nanmean(a,axis=0)
    sd=np.nanstd(a,axis=0,ddof=1)
    sd=np.where(sd>0,sd,np.nan)
    return (a-mu)/sd

def z_within_groups(v, groups):
    v=np.asarray(v,float); groups=np.asarray(groups,object)
    out=np.full(len(v),np.nan,float)
    for g in np.unique(groups):
        ix=np.flatnonzero(groups==g)
        x=v[ix]; m=np.nanmean(x); s=np.nanstd(x,ddof=1)
        out[ix]=(x-m)/s if np.isfinite(s) and s>0 else 0.0
    return out

def load_metadata(model_csv: Path):
    with open(model_csv,newline='',encoding='utf-8') as f:
        return {r['ModelID']:r for r in csv.DictReader(f)}

def expression_rows(expression_csv: Path):
    with open(expression_csv,newline='',encoding='utf-8') as f:
        rr=csv.reader(f); header=next(rr); rows=list(rr)
    return header, rows

def determine_pan_ids(data_dir: Path):
    model_csv=data_dir/'Model(2).csv'
    expr_csv=data_dir/'OmicsExpressionProteinCodingGenesTPMLogp1(2).csv'
    crispr_csv=data_dir/'CRISPRGeneEffect (1).csv'
    if not crispr_csv.exists(): crispr_csv=data_dir/'CRISPRGeneEffect(2).csv'
    meta=load_metadata(model_csv)
    with open(crispr_csv,encoding='utf-8',errors='replace') as f:
        header=f.readline().rstrip('\r\n').split(',')
        cids=[line.split(',',1)[0].strip('"') for line in f if line.strip()]
    with open(expr_csv,newline='',encoding='utf-8') as f:
        rr=csv.reader(f); expr_header=next(rr); eids={r[0] for r in rr}
    ids=[mid for mid in cids if mid in eids and mid in meta and meta[mid].get('OncotreeLineage') in SOLID_LINEAGES]
    return meta, header, expr_header, ids, crispr_csv

def load_lpi_expression(data_dir: Path, ids, meta):
    expr_csv=data_dir/'OmicsExpressionProteinCodingGenesTPMLogp1(2).csv'
    genes=AXIS+FOXO+PROLIF
    with open(expr_csv,newline='',encoding='utf-8') as f:
        rr=csv.reader(f); header=next(rr)
        gene_col={symbol(c):i for i,c in enumerate(header[1:],start=1)}
        missing=[g for g in genes if g not in gene_col]
        if missing: raise RuntimeError(f'Missing expression genes: {missing}')
        idset=set(ids); rows={}
        for r in rr:
            if r and r[0] in idset:
                rows[r[0]]=np.array([float(r[gene_col[g]]) for g in genes],dtype=np.float64)
    X=np.vstack([rows[mid] for mid in ids])
    Z=zscore_cols(X[:,:20])
    raw=.5*(-np.nanmean(Z[:,:len(AXIS)],axis=1))+.5*np.nanmean(Z[:,len(AXIS):],axis=1)
    lineages=np.array([meta[mid]['OncotreeLineage'] for mid in ids],dtype=object)
    lpi=z_within_groups(raw,lineages)
    Zp=zscore_cols(X[:,20:])
    pro=np.nanmean(Zp,axis=1)
    pro=(pro-np.nanmean(pro))/np.nanstd(pro,ddof=1)
    return X, lineages, lpi, pro

def parse_numeric_csv_rest(rest: str, n: int):
    s=rest.rstrip('\r\n')
    if s.startswith(','): s='nan'+s
    if s.endswith(','): s=s+'nan'
    while ',,' in s: s=s.replace(',,',',nan,')
    arr=np.fromstring(s,sep=',',dtype=np.float32)
    if len(arr)!=n:
        raise RuntimeError(f'Expected {n} numeric fields, parsed {len(arr)}')
    return arr

def load_crispr_matrix(data_dir: Path, ids, header, crispr_csv: Path, coverage_fraction=0.90):
    raw_symbols=[symbol(x) for x in header[1:]]
    idset=set(ids)
    nonmiss=np.zeros(len(raw_symbols),dtype=np.int32)
    with open(crispr_csv,encoding='utf-8',errors='replace') as f:
        next(f)
        for line in f:
            mid,rest=line.split(',',1); mid=mid.strip('"')
            if mid not in idset: continue
            fields=rest.rstrip('\r\n').split(',')
            nonmiss += np.fromiter((x!='' and x not in ('NA','NaN','nan') for x in fields),dtype=np.bool_,count=len(fields))
    threshold=math.ceil(coverage_fraction*len(ids))
    sel=np.flatnonzero(nonmiss>=threshold)
    genes=[raw_symbols[i] for i in sel]
    idpos={mid:i for i,mid in enumerate(ids)}
    Y=np.full((len(ids),len(sel)),np.nan,dtype=np.float32)
    with open(crispr_csv,encoding='utf-8',errors='replace') as f:
        next(f)
        for line in f:
            mid,rest=line.split(',',1); mid=mid.strip('"')
            j=idpos.get(mid)
            if j is None: continue
            Y[j,:]=parse_numeric_csv_rest(rest,len(raw_symbols))[sel]
    return genes,Y,nonmiss[sel],threshold

def design_covariates(proliferation, lineage):
    D=pd.get_dummies(np.asarray(lineage,object),drop_first=True,dtype=float).values
    return np.column_stack([np.ones(len(lineage)),np.asarray(proliferation,float),D])

def fit_matrix(Y, x, proliferation, lineage):
    Y=np.asarray(Y)
    n0,g=Y.shape
    C0=design_covariates(proliferation,lineage)
    groups={}
    for j in range(g):
        key=np.isfinite(Y[:,j]).tobytes(); groups.setdefault(key,[]).append(j)
    out=np.empty((g,6),float)
    for key,cols in groups.items():
        ok=np.frombuffer(key,dtype=np.bool_,count=n0)
        yy=Y[ok][:,cols].astype(float); xx=np.asarray(x)[ok]; cc=C0[ok]
        xr=xx-cc@np.linalg.lstsq(cc,xx,rcond=None)[0]
        yr=yy-cc@np.linalg.lstsq(cc,yy,rcond=None)[0]
        denom=float(xr@xr)
        beta=(xr[:,None]*yr).sum(axis=0)/denom
        resid=yr-xr[:,None]*beta
        rank=np.linalg.matrix_rank(np.column_stack([cc,xx])); df=yy.shape[0]-rank
        se=np.sqrt(((resid*resid).sum(axis=0)/df)/denom)
        t=beta/se; p=2*stats.t.sf(np.abs(t),df); partial=t*t/(t*t+df)
        out[cols,:]=np.column_stack([beta,se,t,p,np.repeat(yy.shape[0],len(cols)),partial])
    return out

def fit_response(y, x, proliferation, lineage, robust=False):
    y=np.asarray(y,float); x=np.asarray(x,float); proliferation=np.asarray(proliferation,float); lineage=np.asarray(lineage,object)
    ok=np.isfinite(y)&np.isfinite(x)&np.isfinite(proliferation)
    y=y[ok]; x=x[ok]; proliferation=proliferation[ok]; lineage=lineage[ok]
    C=design_covariates(proliferation,lineage)
    X=np.column_stack([C,x])
    if robust:
        import statsmodels.api as sm
        fit=sm.OLS(y,X).fit(cov_type='HC3')
        b=float(fit.params[-1]); se=float(fit.bse[-1]); t=b/se; p=float(fit.pvalues[-1])
        rank=np.linalg.matrix_rank(X); df=len(y)-rank; partial=t*t/(t*t+df)
    else:
        xr=x-C@np.linalg.lstsq(C,x,rcond=None)[0]
        yr=y-C@np.linalg.lstsq(C,y,rcond=None)[0]
        b=float(xr@yr/(xr@xr)); resid=yr-b*xr
        rank=np.linalg.matrix_rank(X); df=len(y)-rank
        se=float(np.sqrt(((resid@resid)/df)/(xr@xr))); t=b/se; p=float(2*stats.t.sf(abs(t),df)); partial=t*t/(t*t+df)
    return {'n':int(len(y)),'n_lineages':int(len(set(lineage))),'beta':b,'se':se,'t':t,'p':p,'partial_r2':float(partial)}

def discovery_mask(ids,meta):
    return np.array([(meta[mid].get('OncotreePrimaryDisease')=='Invasive Breast Carcinoma') or (meta[mid].get('OncotreeCode')=='LUAD') or (meta[mid].get('OncotreePrimaryDisease')=='Colorectal Adenocarcinoma') for mid in ids],dtype=bool)

def heldout_mask(ids,meta):
    return np.array([meta[mid].get('OncotreeLineage') not in {'Breast','Lung','Bowel'} for mid in ids],dtype=bool)
