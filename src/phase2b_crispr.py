from __future__ import annotations
import argparse, csv, json, math, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from lpi_core import (AXIS,FOXO,PROLIF,SOLID_LINEAGES,symbol,bh_adjust,zscore_cols,z_within_groups,
                      determine_pan_ids,load_lpi_expression,load_crispr_matrix,fit_matrix,
                      design_covariates,discovery_mask,heldout_mask)

SEED=20260915

def frame_from_fit(genes, arr):
    q=bh_adjust(arr[:,3])
    return pd.DataFrame({'gene':genes,'beta':arr[:,0],'se':arr[:,1],'t':arr[:,2],'p':arr[:,3],
                         'q':q,'n':arr[:,4].astype(int),'partial_r2':arr[:,5]})

def lpi_from_20(X20,lineages):
    Z=zscore_cols(X20)
    raw=.5*(-np.nanmean(Z[:,:len(AXIS)],axis=1))+.5*np.nanmean(Z[:,len(AXIS):],axis=1)
    return z_within_groups(raw,lineages)

def stability_table(X20,lineages,base):
    rows=[]
    def add(name,s):
        rows.append({'variant':name,'pearson':float(stats.pearsonr(base,s).statistic),'spearman':float(stats.spearmanr(base,s).statistic)})
    Z=zscore_cols(X20)
    for i,g in enumerate(AXIS+FOXO):
        if i < len(AXIS):
            a=[j for j in range(len(AXIS)) if j!=i]; f=list(range(len(AXIS),len(AXIS)+len(FOXO)))
        else:
            a=list(range(len(AXIS))); f=[j for j in range(len(AXIS),len(AXIS)+len(FOXO)) if j!=i]
        raw=.5*(-np.nanmean(Z[:,a],axis=1))+.5*np.nanmean(Z[:,f],axis=1)
        add('LOO_'+g,z_within_groups(raw,lineages))
    lo=np.nanpercentile(X20,1,axis=0); hi=np.nanpercentile(X20,99,axis=0)
    add('winsor',lpi_from_20(np.clip(X20,lo,hi),lineages))
    med=np.nanmedian(X20,axis=0); mad=np.nanmedian(np.abs(X20-med),axis=0)*1.4826
    mad=np.where(mad>0,mad,np.nan); R=(X20-med)/mad
    raw=.5*(-np.nanmean(R[:,:len(AXIS)],axis=1))+.5*np.nanmean(R[:,len(AXIS):],axis=1)
    add('robust',z_within_groups(raw,lineages))
    Zw=np.full_like(X20,np.nan,dtype=float)
    for g in np.unique(lineages):
        ix=np.flatnonzero(lineages==g); Zw[ix]=zscore_cols(X20[ix])
    raw=.5*(-np.nanmean(Zw[:,:len(AXIS)],axis=1))+.5*np.nanmean(Zw[:,len(AXIS):],axis=1)
    add('within_gene',z_within_groups(raw,lineages))
    add('axis_only',z_within_groups(-np.nanmean(Z[:,:len(AXIS)],axis=1),lineages))
    add('foxo_only',z_within_groups(np.nanmean(Z[:,len(AXIS):],axis=1),lineages))
    return pd.DataFrame(rows)

class MaxTComputer:
    def __init__(self,Y,proliferation,lineages):
        self.n,self.g=Y.shape
        self.C=design_covariates(proliferation,lineages)
        groups={}
        for j in range(self.g):
            key=np.isfinite(Y[:,j]).tobytes(); groups.setdefault(key,[]).append(j)
        self.groups=[]
        for key,cols in groups.items():
            ok=np.frombuffer(key,dtype=np.bool_,count=self.n)
            yy=Y[ok][:,cols].astype(float); cc=self.C[ok]
            yr=yy-cc@np.linalg.lstsq(cc,yy,rcond=None)[0]
            yss=(yr*yr).sum(axis=0)
            rankc=np.linalg.matrix_rank(cc); df=yy.shape[0]-rankc-1
            self.groups.append((ok,cc,yr,yss,df))
    def batch(self,X):
        X=np.asarray(X,float)
        if X.ndim==1: X=X[:,None]
        P=X.shape[1]; maxima=np.zeros(P,float)
        for ok,cc,yr,yss,df in self.groups:
            xx=X[ok,:]
            xr=xx-cc@np.linalg.lstsq(cc,xx,rcond=None)[0]
            xss=(xr*xr).sum(axis=0)
            cross=xr.T@yr
            sse=yss[None,:]-(cross*cross)/xss[:,None]
            sse=np.maximum(sse,np.finfo(float).tiny)
            t=cross/np.sqrt(xss[:,None]*sse/df)
            maxima=np.maximum(maxima,np.nanmax(np.abs(t),axis=1))
        return maxima

def permuted_lpi_matrix(lpi,lineages,n_perm,rng):
    X=np.empty((len(lpi),n_perm),float)
    groups=[np.flatnonzero(lineages==g) for g in np.unique(lineages)]
    for p in range(n_perm):
        v=lpi.copy()
        for ix in groups: v[ix]=v[ix][rng.permutation(len(ix))]
        X[:,p]=v
    return X

def load_all_expression(data_dir,ids):
    path=data_dir/'OmicsExpressionProteinCodingGenesTPMLogp1(2).csv'
    idpos={m:i for i,m in enumerate(ids)}
    with open(path,encoding='utf-8',errors='replace') as f:
        header=f.readline().rstrip('\r\n').split(',')
        genes=[symbol(x) for x in header[1:]]
        E=np.full((len(ids),len(genes)),np.nan,dtype=np.float32)
        for line in f:
            mid,rest=line.split(',',1)
            mid=mid.strip('"')
            j=idpos.get(mid)
            if j is None: continue
            arr=np.fromstring(rest.rstrip('\r\n'),sep=',',dtype=np.float32)
            if len(arr)!=len(genes): raise RuntimeError(f'Expression parse failed for {mid}: {len(arr)} != {len(genes)}')
            E[j]=arr
    return genes,E

def matched_random_scores(data_dir,ids,lineages,n_sig,rng):
    genes,E=load_all_expression(data_dir,ids)
    means=np.nanmean(E,axis=0); sds=np.nanstd(E,axis=0,ddof=1)
    stat=pd.DataFrame({'gene':genes,'mean':means,'sd':sds,'idx':np.arange(len(genes))})
    stat=stat.replace([np.inf,-np.inf],np.nan).dropna().query('sd > 0').copy()
    stat['mean_bin']=pd.qcut(stat['mean'],10,labels=False,duplicates='drop')
    stat['sd_bin']=pd.qcut(stat['sd'],10,labels=False,duplicates='drop')
    lookup=stat.set_index('gene')[['mean_bin','sd_bin']].to_dict('index')
    pools={(a,b):g.idx.to_numpy() for (a,b),g in stat.groupby(['mean_bin','sd_bin'])}
    exclude=set(AXIS+FOXO)
    symbol_by_idx=np.array(genes,object)
    scores=[]; chosen_names=[]
    for _ in range(n_sig):
        chosen=[]
        for target in AXIS+FOXO:
            key=(lookup[target]['mean_bin'],lookup[target]['sd_bin'])
            cand=[int(i) for i in pools[key] if symbol_by_idx[int(i)] not in exclude and int(i) not in chosen]
            if not cand:
                cand=[int(i) for i in stat.idx if symbol_by_idx[int(i)] not in exclude and int(i) not in chosen]
            chosen.append(int(rng.choice(cand)))
        Z=zscore_cols(E[:,chosen])
        raw=.5*(-np.nanmean(Z[:,:len(AXIS)],axis=1))+.5*np.nanmean(Z[:,len(AXIS):],axis=1)
        scores.append(z_within_groups(raw,lineages))
        chosen_names.append('|'.join(symbol_by_idx[chosen]))
    return np.column_stack(scores),chosen_names

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--data-dir',type=Path,required=True); ap.add_argument('--out-dir',type=Path,required=True)
    ap.add_argument('--permutations',type=int,default=300); ap.add_argument('--random-signatures',type=int,default=100); ap.add_argument('--seed',type=int,default=SEED)
    args=ap.parse_args(); args.out_dir.mkdir(parents=True,exist_ok=True)
    meta,header,expr_header,ids,crispr_path=determine_pan_ids(args.data_dir)
    Xsmall,lineages,lpi,pro=load_lpi_expression(args.data_dir,ids,meta)
    genes,Y,coverage,threshold=load_crispr_matrix(args.data_dir,ids,header,crispr_path,0.90)
    disc=discovery_mask(ids,meta); val=~disc; strict=heldout_mask(ids,meta)
    cohorts={'pan':np.ones(len(ids),bool),'discovery':disc,'validation_nonoverlap':val,'validation_lineage_heldout':strict}
    results={}
    for name,mask in cohorts.items():
        arr=fit_matrix(Y[mask],lpi[mask],pro[mask],lineages[mask])
        df=frame_from_fit(genes,arr); df.to_csv(args.out_dir/f'{name}_continuous.csv',index=False); results[name]=df
    score_df=pd.DataFrame({'ModelID':ids,'lineage':lineages,'LPI':lpi,'proliferation':pro,
                           'discovery':disc,'lineage_heldout':strict})
    score_df.to_csv(args.out_dir/'cell_line_scores.csv',index=False)
    stab=stability_table(Xsmall[:,:20],lineages,lpi); stab.to_csv(args.out_dir/'lpi_stability.csv',index=False)

    rng=np.random.default_rng(args.seed)
    maxcomp=MaxTComputer(Y,pro,lineages)
    P=permuted_lpi_matrix(lpi,lineages,args.permutations,rng)
    perm_max=maxcomp.batch(P)
    R,names=matched_random_scores(args.data_dir,ids,lineages,args.random_signatures,rng)
    random_max=maxcomp.batch(R)
    obs=float(np.nanmax(np.abs(results['pan']['t'])))
    perm_p=float((1+np.sum(perm_max>=obs))/(1+len(perm_max)))
    random_p=float((1+np.sum(random_max>=obs))/(1+len(random_max)))
    ctrl=pd.concat([
        pd.DataFrame({'control':'within_lineage_LPI_permutation','iteration':np.arange(1,len(perm_max)+1),'max_abs_t':perm_max}),
        pd.DataFrame({'control':'expression_matched_random_signature','iteration':np.arange(1,len(random_max)+1),'max_abs_t':random_max})
    ],ignore_index=True)
    ctrl.to_csv(args.out_dir/'permutation_controls.csv',index=False)
    pd.DataFrame({'iteration':np.arange(1,len(names)+1),'genes':names,'max_abs_t':random_max}).to_csv(args.out_dir/'random_signature_details.csv',index=False)

    d=results['discovery'].sort_values('p').head(100).set_index('gene')
    v=results['validation_nonoverlap'].set_index('gene'); s=results['validation_lineage_heldout'].set_index('gene'); p=results['pan'].set_index('gene')
    gate=[]
    for g,r in d.iterrows():
        fwer=(1+np.sum(perm_max>=abs(float(p.loc[g,'t']))))/(1+len(perm_max))
        same=(np.sign(r.beta)==np.sign(v.loc[g,'beta'])==np.sign(s.loc[g,'beta']))
        pre=(r.q<.05 and v.loc[g,'q']<.05 and s.loc[g,'q']<.05 and same)
        gate.append({'gene':g,'disc_beta':r.beta,'disc_p':r.p,'disc_q':r.q,'val_beta':v.loc[g,'beta'],'val_q':v.loc[g,'q'],
                     'strict_beta':s.loc[g,'beta'],'strict_q':s.loc[g,'q'],'pan_beta':p.loc[g,'beta'],'pan_q':p.loc[g,'q'],
                     'perm_fwer_p':fwer,'same_direction_both':bool(same),'pre_gate':bool(pre),'final_gate':bool(pre and fwer<.05)})
    pd.DataFrame(gate).to_csv(args.out_dir/'validation_gate.csv',index=False)

    summary={
      'seed':args.seed,'pan_solid_n':len(ids),'lineages':int(len(set(lineages))),'discovery_n':int(disc.sum()),
      'validation_nonoverlap_n':int(val.sum()),'validation_lineage_heldout_n':int(strict.sum()),'genes_tested':len(genes),
      'coverage_threshold_n':int(threshold),'pan_q_lt_0_05':int((results['pan'].q<.05).sum()),
      'discovery_q_lt_0_05':int((results['discovery'].q<.05).sum()),'validation_q_lt_0_05':int((results['validation_nonoverlap'].q<.05).sum()),
      'heldout_q_lt_0_05':int((results['validation_lineage_heldout'].q<.05).sum()),'final_gate_count':int(sum(x['final_gate'] for x in gate)),
      'observed_max_abs_t':obs,'permutation_maxT_empirical_p':perm_p,'random_signature_maxT_empirical_p':random_p,
      'n_permutations':args.permutations,'n_random_signatures':args.random_signatures,
      'ncoa4':results['pan'].set_index('gene').loc['NCOA4',['beta','p','q','partial_r2']].to_dict()
    }
    with open(args.out_dir/'summary.json','w') as f: json.dump(summary,f,indent=2)
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
