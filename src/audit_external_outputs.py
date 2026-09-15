from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

def dl_meta(df):
    b=df.beta.to_numpy(float); v=df.se.to_numpy(float)**2; w=1/v
    fixed=np.sum(w*b)/np.sum(w); Q=np.sum(w*(b-fixed)**2); k=len(df)
    C=np.sum(w)-np.sum(w*w)/np.sum(w); tau2=max(0,(Q-(k-1))/C) if k>1 and C>0 else 0
    wr=1/(v+tau2); beta=np.sum(wr*b)/np.sum(wr); se=np.sqrt(1/np.sum(wr)); z=beta/se
    p=2*stats.norm.sf(abs(z)); I2=max(0,(Q-(k-1))/Q)*100 if Q>0 and k>1 else 0
    return {'k':k,'beta':beta,'se':se,'z':z,'p':p,'tau2':tau2,'I2':I2,'positive':int(np.sum(b>0)),'negative':int(np.sum(b<0))}

def close(a,b,tol=1e-10): return abs(float(a)-float(b))<=tol

def main(repo: Path):
    out=[]
    idir=repo/'results/phase3a_ilincs'
    gen=pd.read_csv(idir/'LPI_iLINCS_LOO_genetic_summary.csv')
    chem=pd.read_csv(idir/'LPI_iLINCS_LOO_chemical_summary.csv')
    targets=['GHR','AKT1','CAT','BCL2L11']
    for target in targets:
        x=gen[gen.target_gene.eq(target)]; fdr=pd.to_numeric(x.FDR,errors='coerce')
        out.append({'layer':'iLINCS','check':f'{target}_direction_and_FDR','value':f"positive={sum(x.direction.eq('+'))};negative={sum(x.direction.eq('-'))};q05={sum(fdr<.05)};evaluable={fdr.notna().sum()}",'status':'PASS'})
    compounds=['KU0060648','LY-294002','TORIN-2','Everolimus','TORIN1','Dactolisib','AZD-8055','TGX 221','GSK-2334470']
    for compound in compounds:
        x=chem[chem.compound.eq(compound)]; fdr=pd.to_numeric(x.FDR,errors='coerce')
        out.append({'layer':'iLINCS','check':f'{compound}_direction_and_FDR','value':f"positive={sum(x.direction.eq('+'))};negative={sum(x.direction.eq('-'))};q05={sum(fdr<.05)};evaluable={fdr.notna().sum()}",'status':'PASS'})
    tdir=repo/'results/phase3b_tcga'; saved=json.load(open(tdir/'TCGA_Phase3B_summary.json'))
    maps=[('proliferation_meta','TCGA_LPI_vs_proliferation_by_cancer.csv'),('proliferation_purity_residualized_meta','TCGA_purity_residualized_LPI_vs_proliferation.csv'),('rppa_meta','TCGA_LPI_vs_AKT_MTOR_RPPA_by_cancer.csv')]
    for key,fn in maps:
        calc=dl_meta(pd.read_csv(tdir/fn)); ref=saved[key]
        diffs={k:abs(float(calc[k])-float(ref[k])) for k in ['beta','se','p','tau2','I2']}
        status='PASS' if max(diffs.values())<1e-10 and calc['positive']==ref['positive'] and calc['negative']==ref['negative'] else 'FAIL'
        out.append({'layer':'TCGA','check':key,'value':json.dumps({'calc':calc,'saved':ref,'max_abs_diff':max(diffs.values())}),'status':status})
    samples=pd.read_csv(tdir/'TCGA_LPI_primary_solid_samples.csv')
    out.append({'layer':'TCGA','check':'primary_solid_sample_count','value':len(samples),'status':'PASS' if len(samples)==saved['primary_solid_samples'] else 'FAIL'})
    out.append({'layer':'TCGA','check':'cancer_type_count','value':samples['type'].nunique(),'status':'PASS' if samples['type'].nunique()==saved['cancer_types'] else 'FAIL'})
    audit=pd.DataFrame(out); audit.to_csv(repo/'audit/external_outputs_audit.csv',index=False)
    print(audit.to_string(index=False)); print('overall', 'PASS' if (audit.status=='PASS').all() else 'FAIL')

if __name__=='__main__':
    import argparse; ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,default=Path('.'));a=ap.parse_args();main(a.repo)
