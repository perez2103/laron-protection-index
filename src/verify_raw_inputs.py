from __future__ import annotations
import argparse, csv, hashlib
from pathlib import Path


def sha256(path: Path, block=8*1024*1024):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        while True:
            b=f.read(block)
            if not b: break
            h.update(b)
    return h.hexdigest()


def main(repo: Path, data_dir: Path) -> int:
    manifest=repo/'data/raw_manifest.csv'
    rows=[]; failed=0
    with open(manifest,newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            p=data_dir/r['filename']
            if not p.exists():
                status='MISSING'; actual=''; failed+=1
            else:
                actual=sha256(p)
                status='PASS' if actual.lower()==r['sha256'].lower() else 'FAIL'
                if status!='PASS': failed+=1
            rows.append((r['source'],r['filename'],status,actual))
    for source,fn,status,actual in rows:
        print(f'{status:7s} {source:8s} {fn}' + (f'  {actual}' if actual else ''))
    print(f'\nVerified {len(rows)-failed}/{len(rows)} manifest inputs.')
    return 0 if failed==0 else 1


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo',type=Path,default=Path('.'))
    ap.add_argument('--data-dir',type=Path,required=True)
    a=ap.parse_args()
    raise SystemExit(main(a.repo,a.data_dir))
