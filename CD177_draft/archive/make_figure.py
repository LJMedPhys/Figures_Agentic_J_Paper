"""Draft Fig. CD177 validation (183 mm wide). Reads only existing run outputs + mask_agreement_per_cell.csv."""
import os, numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],'font.size':6,
    'axes.linewidth':0.5,'xtick.major.width':0.5,'ytick.major.width':0.5,'xtick.major.size':2,'ytick.major.size':2,
    'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
MM=1/25.4; W,H=183,222
fig=plt.figure(figsize=(W*MM,H*MM))
def ax_mm(x,y,w,h,**k): return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H],**k)
def letter(x,y,s): fig.text(x/W,1-y/H,s,fontsize=8,fontweight='bold',va='top')
COL={'published':'#8e3fb0','automated':'#2e9e45','hil':'#1f9fb0','codex':'#7a7a7a'}
FILL={'manual':'#e5c4f5','hf':'#bfeff3','auto':'#c6f6cf','other':'#e6e6e6'}
LAB={'automated':'Agent only','hil':'Agent + human','codex':'Codex CLI'}
M=['automated','hil','codex']

# ---------- a: intent & data
letter(0,0,'a'); a=ax_mm(0,0,52,58); a.axis('off'); a.set_xlim(0,1); a.set_ylim(0,1)
a.text(0.04,0.95,'Do CD177$^+$ and CD177$^-$ neutrophils\ndiffer in morphology?',fontsize=6.5,fontweight='bold',va='top')
a.text(0.04,0.72,'Task',fontweight='bold',va='top')
a.text(0.04,0.64,'Segment cells and multi-lobed nuclei,\nextract morphology, test CD177$^+$ vs CD177$^-$',va='top')
a.text(0.04,0.44,'Data',fontweight='bold',va='top')
a.text(0.04,0.36,'• 177 bright-field images, May–Grünwald stain\n• 10 donors, CD177$^+$ / CD177$^-$ sorted (82 / 95)\n• 2 slide batches with different staining (Top/Bottom)\n• Reference: published analysis$^{\\mathrm{ref}}$, ~1 month analyst',va='top',linespacing=1.5)

# ---------- b: pipelines
letter(55,0,'b'); b=ax_mm(55,0,128,58); b.axis('off'); b.set_xlim(0,128); b.set_ylim(58,0)
steps=['Cell\nsegmentation','Nucleus\nsegmentation','Post-processing\n& QC','Feature\nextraction','Statistics']
feat_n={'automated':len(pd.read_csv('/auto_stats/cd177_feature_statistics.csv')),
        'codex':len(pd.read_csv('/codex_stats/cd177_feature_statistics.csv')),
        'hil':pd.read_csv('/hil_stats/Paired_Donor_Results.csv').metric.nunique()}
rows=[('Published\nreference','published',[('Cellpose (default)','manual'),('ilastik pixel classifier','manual'),('dilation, crop,\nsize filter','manual'),('25 hand-picked','manual'),('cell-level\nMann–Whitney','manual')]),
      ('Agent +\nhuman','hil',[('micro-sam\nvit_b_lm','auto'),('micro-sam fine-tuned\non expert tiles','hf'),('border & lobe QC','auto'),(f"{feat_n['hil']} features",'auto'),('donor-level\nWilcoxon + FDR','auto')]),
      ('Agent\nonly','automated',[('Cellpose cyto3\non optical density','auto'),('threshold\nwithin cell','auto'),('size/shape QC','auto'),(f"{feat_n['automated']} features",'auto'),('donor-level\npermutation + FDR','auto')]),
      ('Codex CLI','codex',[('Cellpose cyto3\n(half resolution)','other'),('erosion/threshold\nwithin cell','other'),('shape QC','other'),(f"{feat_n['codex']} features",'other'),('donor-level\nWilcoxon + FDR','other')])]
x0,cw,gap=19,20.6,1.2
for i,s in enumerate(steps):
    b.add_patch(FancyBboxPatch((x0+i*(cw+gap),2),cw,6.5,boxstyle='round,pad=0,rounding_size=1',fc='white',ec='k',lw=0.4))
    b.text(x0+i*(cw+gap)+cw/2,5.25,s,ha='center',va='center',fontsize=5.5)
for r,(name,key,cells) in enumerate(rows):
    y=11+r*10.3
    b.text(9,y+4.3,name,ha='center',va='center',fontsize=6,fontweight='bold',color=COL[key])
    for i,(t,kind) in enumerate(cells):
        b.add_patch(FancyBboxPatch((x0+i*(cw+gap),y),cw,8.6,boxstyle='round,pad=0,rounding_size=1',fc=FILL[kind],ec='k',lw=0.4))
        b.text(x0+i*(cw+gap)+cw/2,y+4.3,t,ha='center',va='center',fontsize=5.2)
for i,(k,t) in enumerate([('manual','Manual'),('hf','Human feedback'),('auto','Automatic (Agentic-J)'),('other','Automatic (other agent)')]):
    b.add_patch(FancyBboxPatch((19+i*27,53.3),3,2.6,boxstyle='round,pad=0,rounding_size=0.5',fc=FILL[k],ec='k',lw=0.3))
    b.text(23+i*27,54.6,t,va='center',fontsize=5.5)

# ---------- c: gallery (collaborator's existing panel)
letter(0,62,'c'); g=plt.imread('/work/gallery.png'); c=ax_mm(0,63,118,118*g.shape[0]/g.shape[1]); c.imshow(g); c.axis('off')

# ---------- d: agreement with published masks
ag=pd.read_csv('/work/mask_agreement_per_cell.csv'); ag['batch']=np.where(ag.image.str.contains('Top'),'Top','Bottom'); ag['matched']=ag.cell_iou>=0.5
letter(122,62,'d')
mets=[('matched','Reference cells\nrecovered (IoU≥0.5)',lambda d:d.matched.mean()),
      ('cell_iou','Cell IoU\n(matched)',lambda d:d.cell_iou[d.matched].median()),
      ('nuc_dice','Nucleus Dice\n(matched)',lambda d:d.nuc_dice[d.matched].median())]
for j,(k,t,f) in enumerate(mets):
    d=ax_mm(134,70+j*24,48,18)
    for bi,(bt,mk) in enumerate([('Bottom','o'),('Top','s')]):
        for mi,m in enumerate(M):
            sub=ag[(ag.method==m)&(ag.batch==bt)]
            if k!='matched':
                q=sub[k][sub.matched]; lo,hi=q.quantile([.25,.75])
                d.plot([mi+bi*0.25-0.125]*2,[lo,hi],color=COL[m],lw=0.8)
            d.plot(mi+bi*0.25-0.125,f(sub),mk,color=COL[m],ms=3.2,mfc=COL[m] if bt=='Bottom' else 'white',mew=0.7)
    d.set_xticks(range(3)); d.set_xticklabels([LAB[m] for m in M] if j==2 else [],fontsize=5.5,rotation=20,ha='right')
    d.set_ylim(0.5,1.0); d.set_ylabel(t,fontsize=5.5); d.set_xlim(-0.5,2.5); d.grid(axis='y',lw=0.3,color='#dddddd')
d0=fig.axes[-3]
d0.plot([],[],'o',color='k',ms=3,label='Bottom'); d0.plot([],[],'s',color='k',mfc='white',ms=3,label='Top')
d0.legend(frameon=False,fontsize=5,loc='lower right',handletextpad=0.2,borderaxespad=0)
fig.text(126/W,1-63/H,f"Agreement with published masks\n(n = {int((ag.method=='hil').sum())} reference cells; median ± IQR)",fontsize=5.5,va='top')

# ---------- e: per-cell feature agreement
yE=152; letter(0,yE,'e')
mt=ag[ag.matched]
for r,(pk,rk,lab,lim) in enumerate([('pred_cell_px','ref_cell_px','Cell area (px$^2$)',(0,9000)),('pred_nuc_px','ref_nuc_px','Nucleus area (px$^2$)',(0,4000))]):
    for mi,m in enumerate(M):
        e=ax_mm(9+mi*22,yE+4+r*24,18,18); s=mt[mt.method==m]; x,y=s[rk].values,s[pk].values
        e.hexbin(x,y,gridsize=30,extent=lim*2,cmap=mpl.colors.LinearSegmentedColormap.from_list('c',['#ffffff',COL[m]]),bins='log',mincnt=1,lw=0)
        e.plot(lim,lim,'k--',lw=0.4)
        ccc=2*np.cov(x,y)[0,1]/(x.var()+y.var()+(x.mean()-y.mean())**2)
        e.text(0.04,0.96,f'CCC {ccc:.2f}',transform=e.transAxes,va='top',fontsize=5)
        e.set_xlim(lim); e.set_ylim(lim); e.set_xticks([0,lim[1]/2,lim[1]]); e.set_yticks([0,lim[1]/2,lim[1]]); e.tick_params(labelsize=4.5)
        if r==0: e.set_title(LAB[m],fontsize=5.5,color=COL[m],pad=2)
        if mi==0: e.set_ylabel(lab.replace('area','area,\nagent'),fontsize=5)
        else: e.set_yticklabels([])
        if r==1: e.set_xlabel('published',fontsize=5)
        else: e.set_xticklabels([])

# ---------- f: donor-level CD177+ vs CD177- differences
def pct(pos,neg): return 100*(pos-neg)/neg
auto=pd.read_csv('/auto_stats/donor_contrasts.csv')
auto['d']=pct(auto.donor_pos_summary,auto.donor_neg_summary)
hil=pd.read_csv('/hil_stats/Donor_Level_Metrics.csv'); hil=hil[hil.session_subset=='all']
hp=hil.pivot_table(index=['metric','donor'],columns='cd177_status',values='donor_median_of_image_medians').reset_index(); hp['d']=pct(hp.pos,hp.neg)
cx=pd.read_csv('/codex_cells/cd177_cell_features.csv'); cx=cx[~cx.magnification.astype(str).str.contains('60')]
FEATS=[('Cell area',{'automated':'cell_area','hil':'cell_area_um2','codex':'cell_area'},'p = 3.5×10$^{-52}$, r = 0.27'),
       ('Nucleus area',{'automated':'nucleus_area','hil':'nucleus_area_um2','codex':'nucleus_area'},'p = 0.2, r = 0.02'),
       ('Nucleus / cell area',{'automated':'nuclear_fraction_of_cell','hil':'nucleus_cell_area_ratio','codex':'nucleus_cell_area_ratio'},'p = 2.4×10$^{-11}$, r = 0.12')]
qa=pd.read_csv('/auto_stats/cd177_feature_statistics.csv').set_index('feature').p_adjusted
qh=pd.read_csv('/hil_stats/Paired_Donor_Results.csv'); qh=qh[qh.get('session_subset','all')=='all'] if 'session_subset' in qh else qh; qh=qh.drop_duplicates('metric').set_index('metric').p_fdr_bh
qc=pd.read_csv('/codex_stats/cd177_feature_statistics.csv').set_index('feature').q_value_bh
letter(75,yE,'f'); rng=np.random.default_rng(0); out=[]
for fi,(name,fm,pub) in enumerate(FEATS):
    f=ax_mm(88+fi*32,yE+6,27,36)
    for mi,m in enumerate(M):
        if m=='automated': v=auto[auto.feature==fm[m]].d.values; q=qa[fm[m]]
        elif m=='hil': v=hp[hp.metric==fm[m]].d.values; q=qh[fm[m]]
        else:
            g2=cx.groupby(['donor','slide','condition'])[fm[m]].median().unstack('condition')
            pos=[c for c in g2.columns if 'pos' in str(c).lower()][0]; neg=[c for c in g2.columns if 'neg' in str(c).lower()][0]
            v=pct(g2[pos],g2[neg]).dropna().values; q=qc[fm[m]]
        f.scatter(mi+rng.uniform(-.15,.15,len(v)),v,s=5,color=COL[m],alpha=.8,lw=0)
        f.plot([mi-.25,mi+.25],[np.median(v)]*2,color='k',lw=0.8)
        f.text(mi,1.02,f'q={q:.2f}',transform=f.get_xaxis_transform(),ha='center',fontsize=4.8)
        out.append(dict(feature=name,method=m,n_units=len(v),median_pct_diff=np.median(v),q=q))
    f.axhline(0,color='#999999',lw=0.4,ls='--'); f.set_xlim(-.6,2.6)
    f.set_xticks(range(3)); f.set_xticklabels([LAB[m] for m in M],rotation=30,ha='right',fontsize=5)
    f.set_title(name,fontsize=6,pad=9); f.text(0.5,-0.36,'published (cell-level):\n'+pub,transform=f.transAxes,ha='center',va='top',fontsize=5,color=COL['published'])
    if fi==0: f.set_ylabel('CD177$^+$ vs CD177$^-$ (% difference\nper donor or donor·slide)',fontsize=5.5)
pd.DataFrame(out).to_csv('/work/panel_f_values.csv',index=False)
for ext in ['pdf','svg','png']: fig.savefig(f'/work/Fig_CD177_draft.{ext}',dpi=450 if ext=='png' else None)
print(pd.DataFrame(out).round(3).to_string())
