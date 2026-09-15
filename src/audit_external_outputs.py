from __future__ import annotations
import json
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
    return {'k':k,'beta':beta,'se':se,'z':z,'p':p,'tau2':tau2,'I2':I2,
            'positive':int(np.sum(b>0)),'negative':int(np.sum(b<0))}


def main(repo: Path):
    out=[]

    # GitHub stores a compact, paper-facing iLINCS robustness table. The complete
    # query/LOO outputs are retained in the archived v1.0.0 package for exact provenance.
    idir=repo/'results/phase3a_ilincs'
    key=pd.read_csv(idir/'key_robustness_summary.csv')
    expected={
      'GHR':('genetic',9,9,0,7), 'AKT1':('genetic',7,7,0,5),
      'CAT':('genetic',8,0,8,8), 'BCL2L11':('genetic',3,0,3,2),
      'KU0060648':('chemical',9,9,0,9), 'LY-294002':('chemical',9,9,0,9),
      'TORIN-2':('chemical',8,8,0,8), 'Everolimus':('chemical',9,8,1,8),
      'TORIN1':('chemical',6,6,0,6), 'Dactolisib':('chemical',7,7,0,5),
      'AZD-8055':('chemical',4,4,0,4), 'TGX 221':('chemical',9,0,9,9),
      'GSK-2334470':('chemical',5,0,5,5)
    }
    for pert,(kind,evaluable,pos,neg,q05) in expected.items():
        x=key[(key['class']==kind)&(key.perturbagen==pert)]
        ok=(len(x)==1 and int(x.evaluable.iloc[0])==evaluable and int(x.positive.iloc[0])==pos and
            int(x.negative.iloc[0])==neg and int(x.fdr_lt_0_05.iloc[0])==q05)
        value='missing' if x.empty else f"evaluable={int(x.evaluable.iloc[0])};positive={int(x.positive.iloc[0])};negative={int(x.negative.iloc[0])};q05={int(x.fdr_lt_0_05.iloc[0])}"
        out.append({'layer':'iLINCS','check':f'{pert}_LOO_robustness','value':value,'status':'PASS' if ok else 'FAIL'})

    # Recompute the three paper-level TCGA random-effects meta-analyses directly
    # from the versioned per-cancer source tables.
    tdir=repo/'results/phase3b_tcga'; saved=json.load(open(tdir/'TCGA_Phase3B_summary.json'))
    maps=[('proliferation_meta','TCGA_LPI_vs_proliferation_by_cancer.csv'),
          ('proliferation_purity_residualized_meta','TCGA_purity_residualized_LPI_vs_proliferation.csv'),
          ('rppa_meta','TCGA_LPI_vs_AKT_MTOR_RPPA_by_cancer.csv')]
    for key_name,fn in maps:
        calc=dl_meta(pd.read_csv(tdir/fn)); ref=saved[key_name]
        diffs={k:abs(float(calc[k])-float(ref[k])) for k in ['beta','se','p','tau2','I2']}
        ok=max(diffs.values())<1e-10 and calc['positive']==ref['positive'] and calc['negative']==ref['negative']
        out.append({'layer':'TCGA','check':key_name,
                    'value':json.dumps({'calc':calc,'saved':ref,'max_abs_diff':max(diffs.values())}),
                    'status':'PASS' if ok else 'FAIL'})

    # The 9,632-sample table is in the archived full release; the lightweight
    # GitHub checkout records its audited count in the frozen summary JSON.
    out.append({'layer':'TCGA','check':'primary_solid_sample_count_recorded',
                'value':saved['primary_solid_samples'],
                'status':'PASS' if saved['primary_solid_samples']==9632 else 'FAIL'})
    out.append({'layer':'TCGA','check':'cancer_type_count_recorded',
                'value':saved['cancer_types'],
                'status':'PASS' if saved['cancer_types']==31 else 'FAIL'})

    audit=pd.DataFrame(out); (repo/'audit').mkdir(exist_ok=True)
    audit.to_csv(repo/'audit/external_outputs_audit.csv',index=False)
    print(audit.to_string(index=False))
    passed=bool((audit.status=='PASS').all())
    print('\nOVERALL', 'PASS' if passed else 'FAIL')
    return 0 if passed else 1


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('--repo',type=Path,default=Path('.')); a=ap.parse_args()
    raise SystemExit(main(a.repo))
