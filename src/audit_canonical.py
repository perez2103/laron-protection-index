from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd


def main(repo: Path) -> int:
    rows=[]
    p2=repo/'results/phase2b'
    s2=json.load(open(p2/'summary.json'))
    checks={
        'pan_solid_n':(s2.get('pan_solid_n'),914),
        'discovery_n':(s2.get('discovery_n'),152),
        'validation_nonoverlap_n':(s2.get('validation_nonoverlap_n'),762),
        'validation_lineage_heldout_n':(s2.get('validation_lineage_heldout_n'),688),
        'genes_tested':(s2.get('genes_tested'),17787),
        'final_gate_count':(s2.get('final_gate_count'),0),
        'seed':(s2.get('seed'),20260915),
    }
    ok=all(a==b for a,b in checks.values()) and s2['permutation_maxT_empirical_p']>0.05 and s2['random_signature_maxT_empirical_p']>0.05
    rows.append({'layer':'Phase2B','check':'canonical summary','status':'PASS' if ok else 'FAIL','detail':json.dumps({'checks':checks,'maxT_p':s2['permutation_maxT_empirical_p'],'random_signature_p':s2['random_signature_maxT_empirical_p']})})

    p3=repo/'results/phase3a_pharmacogenomics'
    s3=json.load(open(p3/'summary.json'))
    expect={'prism_models':637,'prism_compounds':6790,'prism_global_q05':0,'prism_family_q05':0,
            'gdsc1_models_overlap':448,'gdsc1_drugs':402,'gdsc1_global_q05':0,'gdsc1_family_q05':1,
            'gdsc2_models_overlap':447,'gdsc2_drugs':295,'gdsc2_global_q05':0,'gdsc2_family_q05':0}
    ok=all(s3.get(k)==v for k,v in expect.items())
    rows.append({'layer':'Phase3A','check':'canonical pharmacogenomic summary','status':'PASS' if ok else 'FAIL','detail':json.dumps({k:s3.get(k) for k in expect})})

    pr=pd.read_csv(p3/'PRISM_results.csv'); g1=pd.read_csv(p3/'GDSC1_results.csv'); g2=pd.read_csv(p3/'GDSC2_results.csv')
    observed={'PRISM_global':int((pr.q_global<.05).sum()),'PRISM_family':int((pr.q_family<.05).sum()),
              'GDSC1_global':int((g1.q_global<.05).sum()),'GDSC1_family':int((g1.q_family<.05).sum()),
              'GDSC2_global':int((g2.q_global<.05).sum()),'GDSC2_family':int((g2.q_family<.05).sum())}
    expected={'PRISM_global':0,'PRISM_family':0,'GDSC1_global':0,'GDSC1_family':1,'GDSC2_global':0,'GDSC2_family':0}
    rows.append({'layer':'Phase3A','check':'canonical FDR gate counts','status':'PASS' if observed==expected else 'FAIL','detail':json.dumps(observed)})

    ext_path=repo/'audit/external_outputs_audit.csv'
    ext=pd.read_csv(ext_path)
    rows.append({'layer':'External','check':'iLINCS/TCGA processed-output audit','status':'PASS' if (ext.status=='PASS').all() else 'FAIL','detail':f'{len(ext)} checks; failures={(ext.status!="PASS").sum()}'})

    audit=pd.DataFrame(rows)
    out=repo/'audit/canonical_release_audit.csv'
    audit.to_csv(out,index=False)
    print(audit.to_string(index=False))
    passed=bool((audit.status=='PASS').all())
    print('\nOVERALL', 'PASS' if passed else 'FAIL')
    return 0 if passed else 1


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo',type=Path,default=Path('.'))
    args=ap.parse_args()
    raise SystemExit(main(args.repo))
