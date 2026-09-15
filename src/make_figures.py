from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures'; OUT.mkdir(exist_ok=True)

def save(fig,name):
    fig.savefig(OUT/f'{name}.png',dpi=300,bbox_inches='tight')
    fig.savefig(OUT/f'{name}.svg',bbox_inches='tight')
    plt.close(fig)

def dl_meta(df):
    b=df.beta.to_numpy(float);v=df.se.to_numpy(float)**2;w=1/v;fixed=(w*b).sum()/w.sum();Q=(w*(b-fixed)**2).sum();k=len(df);C=w.sum()-(w*w).sum()/w.sum();tau=max(0,(Q-(k-1))/C) if C>0 else 0;wr=1/(v+tau);beta=(wr*b).sum()/wr.sum();se=(1/wr.sum())**.5;p=2*norm.sf(abs(beta/se));I2=max(0,(Q-(k-1))/Q)*100 if Q>0 else 0;return beta,se,p,I2

def forest(ax,df,title,xlabel):
    d=df.sort_values('beta').reset_index(drop=True); y=np.arange(len(d));lo=d.beta-1.96*d.se;hi=d.beta+1.96*d.se
    ax.errorbar(d.beta,y,xerr=[d.beta-lo,hi-d.beta],fmt='o',markersize=3,capsize=1.5,linewidth=.8);ax.axvline(0,linewidth=.8);ax.set_yticks(y,d.type,fontsize=7);ax.set_xlabel(xlabel,fontsize=8);ax.set_title(title,fontsize=10);m=dl_meta(d);ax.text(.02,.01,f'RE meta beta={m[0]:.3f}; P={m[2]:.3g}; I2={m[3]:.1f}%',transform=ax.transAxes,fontsize=7)

def main():
    p2=ROOT/'results/phase2b';p3=ROOT/'results/phase3a_pharmacogenomics';il=ROOT/'results/phase3a_ilincs';tc=ROOT/'results/phase3b_tcga'
    stab=pd.read_csv(p2/'lpi_stability.csv')
    fig,ax=plt.subplots(figsize=(8,4.8));ax.bar(np.arange(len(stab)),stab.pearson);ax.set_xticks(np.arange(len(stab)),stab.variant,rotation=60,ha='right',fontsize=7);ax.set_ylim(0,1.05);ax.set_ylabel('Pearson r with frozen LPI');ax.set_title('Figure 1 | LPI stability');fig.tight_layout();save(fig,'Figure1_LPI_stability')
    disc=pd.read_csv(p2/'discovery_continuous.csv');val=pd.read_csv(p2/'validation_nonoverlap_continuous.csv');pan=pd.read_csv(p2/'pan_continuous.csv');m=disc[['gene','beta','q']].rename(columns={'beta':'disc_beta','q':'disc_q'}).merge(val[['gene','beta','q']].rename(columns={'beta':'val_beta','q':'val_q'}),on='gene').merge(pan[['gene','q']].rename(columns={'q':'pan_q'}),on='gene')
    fig,axs=plt.subplots(1,2,figsize=(12,5));axs[0].scatter(m.disc_beta,m.val_beta,s=7,alpha=.25);lim=np.nanmax(np.abs(np.r_[m.disc_beta,m.val_beta]));axs[0].plot([-lim,lim],[-lim,lim],linewidth=.8);axs[0].axhline(0,linewidth=.7);axs[0].axvline(0,linewidth=.7);axs[0].set_xlabel('Discovery beta (n=152)');axs[0].set_ylabel('Non-overlap validation beta (n=762)');axs[0].set_title('a  Discovery vs validation',loc='left',fontweight='bold');axs[0].text(.02,.98,'Discovery FDR<0.05: 0\nFinal gate: 0',transform=axs[0].transAxes,va='top',fontsize=8)
    top=pan.sort_values('q').head(12).sort_values('q',ascending=False);axs[1].barh(np.arange(len(top)),-np.log10(top.q));axs[1].set_yticks(np.arange(len(top)),top.gene,fontsize=8);axs[1].axvline(-np.log10(.05),linewidth=.8);axs[1].set_xlabel('-log10(FDR)');axs[1].set_title('b  Pan-solid signals',loc='left',fontweight='bold');fig.tight_layout();save(fig,'Figure2_CRISPR_no_validated_dependency')
    cross=pd.read_csv(p3/'cross_platform_family.csv');fig,ax=plt.subplots(figsize=(6.5,5.5));
    for platform,g in cross.groupby('validation_platform'):ax.scatter(g.prism_beta,g.gdsc_beta_median,s=30,alpha=.7,label=platform)
    ax.axhline(0,linewidth=.7);ax.axvline(0,linewidth=.7);ax.set_xlabel('PRISM beta');ax.set_ylabel('GDSC median beta');ax.legend(frameon=False);ax.set_title('Figure 3 | No replicated pharmacogenomic liability');fig.tight_layout();save(fig,'Figure3_pharmacogenomics')
    gs=pd.read_csv(il/'LPI_iLINCS_LOO_genetic_summary.csv');cs=pd.read_csv(il/'LPI_iLINCS_LOO_chemical_summary.csv')
    gt=[]
    for target in ['GHR','AKT1','CAT','BCL2L11']:
        x=gs[gs.target_gene.eq(target)];f=pd.to_numeric(x.FDR,errors='coerce');pos=sum(x.direction.eq('+'));neg=sum(x.direction.eq('-'));gt.append((target,(1 if pos>=neg else -1)*sum(f<.05)))
    ct=[]
    for c in ['KU0060648','LY-294002','TORIN-2','Everolimus','TORIN1','Dactolisib','AZD-8055','TGX 221','GSK-2334470']:
        x=cs[cs.compound.eq(c)];f=pd.to_numeric(x.FDR,errors='coerce');pos=sum(x.direction.eq('+'));neg=sum(x.direction.eq('-'));ct.append((c,(1 if pos>=neg else -1)*sum(f<.05)))
    fig,axs=plt.subplots(1,2,figsize=(12,5));gt=sorted(gt,key=lambda z:z[1]);axs[0].barh(range(len(gt)),[x[1] for x in gt]);axs[0].set_yticks(range(len(gt)),[x[0] for x in gt]);axs[0].axvline(0,linewidth=.7);axs[0].set_title('a  Genetic perturbations',loc='left',fontweight='bold');ct=sorted(ct,key=lambda z:z[1]);axs[1].barh(range(len(ct)),[x[1] for x in ct]);axs[1].set_yticks(range(len(ct)),[x[0] for x in ct],fontsize=8);axs[1].axvline(0,linewidth=.7);axs[1].set_title('b  Chemical perturbations',loc='left',fontweight='bold');fig.suptitle('Figure 4 | iLINCS leave-one-landmark-out robustness');fig.tight_layout();save(fig,'Figure4_iLINCS_robustness')
    prol=pd.read_csv(tc/'TCGA_LPI_vs_proliferation_by_cancer.csv');rppa=pd.read_csv(tc/'TCGA_LPI_vs_AKT_MTOR_RPPA_by_cancer.csv');fig,axs=plt.subplots(1,2,figsize=(12,9));forest(axs[0],prol,'a  Proliferation','Purity-adjusted beta');forest(axs[1],rppa,'b  AKT-mTOR phosphosignaling','Purity-adjusted beta');fig.suptitle('Figure 5 | TCGA separates phosphosignaling from proliferation');fig.tight_layout(rect=[0,0,1,.97]);save(fig,'Figure5_TCGA_molecular_validation')
    fig,ax=plt.subplots(figsize=(10,5.5));ax.axis('off');ax.text(.5,.9,'Integrated evidence',ha='center',fontsize=16,fontweight='bold');items=[('CRISPR','0 validated targets'),('PRISM/GDSC','0 replicated sensitivities'),('iLINCS','perturbational convergence'),('TCGA','lower AKT-mTOR RPPA\nno pan-cancer proliferation effect')];xs=[.12,.37,.62,.87]
    for x,(h,b) in zip(xs,items):ax.text(x,.62,h+'\n'+b,ha='center',va='center',fontsize=9,bbox=dict(boxstyle='round',facecolor='none'))
    for a,b in zip(xs[:-1],xs[1:]):ax.annotate('',xy=(b-.08,.62),xytext=(a+.08,.62),arrowprops=dict(arrowstyle='->'))
    ax.text(.5,.25,'Reproducible molecular state != targetable cancer vulnerability',ha='center',fontsize=13,fontweight='bold');save(fig,'Figure6_integrated_model')

if __name__=='__main__':main()
