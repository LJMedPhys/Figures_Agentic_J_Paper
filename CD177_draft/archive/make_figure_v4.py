"""Half-page CD177 validation figure, v4 (adds fast mode; step roles from INTERVENTIONS_HIL1.md).
Inputs: run outputs, collaborator gallery PNG, mask_agreement_per_cell.csv, nuclei_vs_expert.csv."""
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],'font.size':5,
    'axes.linewidth':0.4,'xtick.major.width':0.4,'ytick.major.width':0.4,'xtick.major.size':1.5,'ytick.major.size':1.5,
    'xtick.labelsize':4.5,'ytick.labelsize':4.5,'axes.labelsize':5,'axes.spines.top':False,'axes.spines.right':False,
    'svg.fonttype':'none','pdf.fonttype':42,'xtick.major.pad':1.2,'ytick.major.pad':1.2,'axes.labelpad':1.5})
MM=1/25.4; W,H=183,129
fig=plt.figure(figsize=(W*MM,H*MM))
def ax_mm(x,y,w,h,**k): return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H],**k)
def letter(x,y,s): fig.text(x/W,1-y/H,s,fontsize=8,fontweight='bold',va='top')
# method colours (data marks) and step-role fills (panel a), chosen to be mutually distinct
COL={'published':'#7b3294','hil':'#0571b0','automated':'#1a9641','fast':'#a6611a','codex':'#636363','hil_stock':'#0571b0'}
FILL={'manual':'#f7c873','hf':'#8ecdf2','ai':'#f4a3c0','auto':'#a6dba0','other':'#d9d9d9'}
LAB={'published':'Published','automated':'Agent only (advanced)','fast':'Agent only (fast)','hil':'Agent + human','codex':'Codex CLI','hil_stock':'before fine-tuning'}
SH={'automated':'Advanced','fast':'Fast','hil':'Agent +\nhuman','codex':'Codex'}
M=['automated','fast','hil','codex']

# ---------------- a
letter(0,0,'a'); a=ax_mm(0,0,183,49); a.axis('off'); a.set_xlim(0,183); a.set_ylim(49,0)
steps=['Cell\nsegmentation','Nucleus\nsegmentation','QC &\npost-processing','Features','Statistics']
rows=[('Published','published',[('Cellpose (default)','manual'),('ilastik pixel\nclassifier','manual'),('dilation, crop,\nsize filter','manual'),('25 hand-picked','manual'),('cell-level\nMann–Whitney','manual')]),
      ('Agent + human','hil',[('micro-sam vit_b_lm\n(settings by AI assistant)','ai'),(('15 tiles corrected by user','hf'),('retrained by AI assistant','ai')),('debris & border QC\n(set by AI assistant)','ai'),('split-invariant metrics\n(set by AI assistant)','ai'),('paired Wilcoxon\n(set by AI assistant)','ai')]),
      ('Agent only\n(advanced)','automated',[('Cellpose cyto3\non optical density','auto'),('threshold\nwithin cell','auto'),('size/shape QC','auto'),('32 features','auto'),('donor-level\npermutation + FDR','auto')]),
      ('Agent only\n(fast)','fast',[('Cellpose cyto3\n(diameter 50 px)','auto'),('Otsu threshold\nwithin cell','auto'),('border, size &\nsolidity filter','auto'),('25 features','auto'),('donor-level\nWilcoxon + FDR','auto')]),
      ('Codex CLI','codex',[('Cellpose cyto3\n(half resolution)','other'),('erosion/threshold\nwithin cell','other'),('shape QC','other'),('30 features','other'),('donor-level\nWilcoxon + FDR','other')])]
TIME={'published':(3*7*24*60,'weeks'),'hil':(345,'5.8 h*'),'automated':(67,'67 min'),'fast':(35,'35 min'),'codex':(13,'13 min')}
VERD={'published':('ok','no meaningful difference\n(small effect sizes)'),
      'hil':('ok','0/11 features significant\n(FDR); cell-area trend'),
      'automated':('ok','0/32 features significant (FDR);\nsame in 3/3 replicates'),
      'fast':('x','5/25 significant (cell shape);\ndifferences in 3/3 replicates'),
      'codex':('ok','0/30 features significant\n(FDR)')}
x0,cw,gap,rh,y0=17,18.4,0.5,6.2,6
for i,s_ in enumerate(steps): a.text(x0+i*(cw+gap)+cw/2,3,s_,ha='center',va='center',fontsize=5,fontweight='bold',linespacing=1.0)
a.text(126,3,'Time to result',ha='center',va='center',fontsize=5,fontweight='bold'); a.text(163,3,'Scientific verdict',ha='center',va='center',fontsize=5,fontweight='bold')
for r,(name,key,cells) in enumerate(rows):
    y=y0+r*(rh+0.6)
    a.text(0.5,y+rh/2,name,va='center',fontsize=5.2,fontweight='bold',color=COL[key],linespacing=1.0)
    for i,c in enumerate(cells):
        xx=x0+i*(cw+gap)
        if isinstance(c[0],tuple):
            for j,(tt,kind) in enumerate(c):
                a.add_patch(Rectangle((xx,y+j*rh/2),cw,rh/2,fc=FILL[kind],ec='k',lw=0.35))
                a.text(xx+cw/2,y+j*rh/2+rh/4,tt,ha='center',va='center',fontsize=4.0,linespacing=0.9)
        else:
            a.add_patch(Rectangle((xx,y),cw,rh,fc=FILL[c[1]],ec='k',lw=0.35))
            a.text(xx+cw/2,y+rh/2,c[0],ha='center',va='center',fontsize=4.4,linespacing=1.0)
    ok,txt=VERD[key]
    a.text(149,y+rh/2,'\u2713' if ok=='ok' else '\u2717',ha='center',va='center',fontsize=7,fontweight='bold',
           color=COL[key] if ok=='ok' else '#c0392b',fontfamily='DejaVu Sans')
    a.text(152,y+rh/2,txt,va='center',fontsize=4.4,linespacing=1.05)
ta=ax_mm(112,y0,28,5*(rh+0.6)-0.6)
for r,(name,key,_) in enumerate(rows):
    v,lab=TIME[key]; ta.barh(r,v,height=0.55,color=COL[key],lw=0); ta.text(v*1.25,r,lab,va='center',fontsize=4.4)
ta.set_xscale('log'); ta.set_xlim(5,3e5); ta.set_ylim(4.5,-0.5); ta.set_yticks([]); ta.minorticks_off()
ta.set_xticks([10,60,1440,10080]); ta.set_xticklabels(['10 min','1 h','1 day','1 week'],fontsize=4.2); ta.spines['left'].set_visible(False)
xl=x0
for k,t in [('manual','Manual (analyst)'),('hf','Human feedback (user)'),('ai','AI assistant (Claude Code)'),('auto','Automatic (Agentic-J)'),('other','Automatic (other agent)')]:
    a.add_patch(Rectangle((xl,43.3),2.2,2,fc=FILL[k],ec='k',lw=0.3)); a.text(xl+3,44.3,t,va='center',fontsize=4.5); xl+=3.2+len(t)*0.88+3
a.text(x0,47.6,'All agent runs started from an expert-written task description. AI assistant = Claude Code, sending instructions through the user\'s chat (HIL1 intervention log). *wall-clock incl. human time.',fontsize=4.2,color='#555555',va='center')

# ---------------- b gallery (all four rows)
yb=51; letter(0,yb,'b'); g=plt.imread('/work/gallery.png'); sc=g.shape[1]/2000
cols=[212,512,813,1114,1414,1715]; rws=[100,421,743,1064]; T=264
clab=['raw','Published','Agent only\n(advanced)','Agent only\n(fast)','Agent +\nhuman','Codex CLI']
ccol=['k',COL['published'],COL['automated'],COL['fast'],COL['hil'],COL['codex']]
rlab=['V8 neg\nBottom','V10 pos\nBottom','V8 neg\nTop','V9 pos\nTop']
notes={(0,5):'missed',(1,1):'debris as cells',(2,1):'nucleus only',(2,2):'missed',(3,3):'wrong split'}
tw,tg,rp=12.4,0.45,14.4
for ri,yy in enumerate(rws):
    for ci,xx in enumerate(cols):
        t=g[int(yy*sc):int((yy+T)*sc),int(xx*sc):int((xx+T)*sc)]
        ax=ax_mm(8+ci*(tw+tg),yb+5.5+ri*rp,tw,tw); ax.imshow(t); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(0.3)
        if ri==0: ax.set_title(clab[ci],fontsize=4.5,color=ccol[ci],pad=1.2,linespacing=0.95)
        if ci==0: ax.set_ylabel(rlab[ri],fontsize=4.4,rotation=0,ha='right',va='center',labelpad=1.5)
        if (ri,ci) in notes: ax.text(0.5,-0.02,notes[(ri,ci)],transform=ax.transAxes,ha='center',va='top',fontsize=4,color='#b0172a')
lx=ax_mm(8,yb+64.5,80,3); lx.axis('off'); lx.set_xlim(0,80); lx.set_ylim(0,1)
for x,c,ls,t in [(8,'#f2c200','-','cell'),(22,'#27c7f0','-','nucleus'),(38,'#999999','--','segmented, not measured')]:
    lx.plot([x,x+4],[0.5,0.5],color=c,ls=ls,lw=1); lx.text(x+5,0.5,t,va='center',fontsize=4.6)

# ---------------- c nucleus accuracy vs expert
letter(88,yb,'c'); nx=pd.read_csv('/work/nuclei_vs_expert.csv'); th=np.arange(0.5,0.85,0.05)
n_nuc=int(nx[nx.method=='hil'].n_true.sum()); n_t=nx.tile.nunique()
fig.text(96/W,1-(yb+1.2)/H,f'Nuclei vs expert annotation\n({n_nuc} nuclei, {n_t} held-out tiles)',fontsize=4.8,va='top',linespacing=1.0)
c=ax_mm(96,yb+7,29,22)
for m,ls,lw in [('hil','-',1.1),('codex','-',0.8),('published','-',1.1),('fast','-',0.8),('automated','-',0.8),('hil_stock',':',0.8)]:
    s=nx[nx.method==m]; c.plot(th,[s[f'acc_{t:.2f}'].mean() for t in th],ls,color=COL[m],lw=lw,label=f"{LAB[m]}  {s.msa.mean():.2f}")
c.set_xlim(0.5,0.8); c.set_ylim(0,0.6); c.set_xticks([0.5,0.6,0.7,0.8]); c.set_yticks([0,0.2,0.4,0.6])
c.set_xlabel('IoU threshold (stricter →)'); c.set_ylabel('Nucleus accuracy')
c.legend(frameon=False,fontsize=3.7,loc='upper right',bbox_to_anchor=(1.04,1.06),handlelength=1.0,borderaxespad=0,labelspacing=0.18,title='mSA',title_fontsize=3.7,alignment='right')

# ---------------- d agreement with published masks
letter(128,yb,'d'); ag=pd.read_csv('/work/mask_agreement_per_cell.csv'); ag['batch']=np.where(ag.image.str.contains('Top'),'Top','Bottom'); ag['matched']=ag.cell_iou>=0.5
fig.text(136/W,1-(yb+1.2)/H,f'Agreement with published masks\n(n = {int((ag.method=="hil").sum())} cells)',fontsize=4.8,va='top',linespacing=1.0)
fig.legend(handles=[Line2D([],[],marker='o',ls='',color='k',ms=2.6,label='Bottom slides'),Line2D([],[],marker='s',ls='',color='k',mfc='white',ms=2.6,label='Top slides')],
           frameon=False,fontsize=4.3,loc='upper right',bbox_to_anchor=(182/W,1-(yb+0.8)/H),handletextpad=0.1,labelspacing=0.2,borderaxespad=0)
mets=[('matched','Published cells\nfound (IoU ≥ 0.5)'),('cell_iou','Cell overlap\n(IoU)'),('nuc_dice','Nucleus overlap\n(Dice)')]
for j,(k,t) in enumerate(mets):
    d=ax_mm(137+j*16,yb+10,13.5,17)
    for mi,m in enumerate(M):
        for bi,(bt,mk) in enumerate([('Bottom','o'),('Top','s')]):
            s=ag[(ag.method==m)&(ag.batch==bt)]; xp=mi+(bi-0.5)*0.38
            if k=='matched': val=s.matched.mean()
            else:
                q=s[k][s.matched]; val=q.median(); d.plot([xp,xp],q.quantile([.25,.75]),color=COL[m],lw=0.6)
            d.plot(xp,val,mk,color=COL[m],ms=2.3,mfc=COL[m] if bt=='Bottom' else 'white',mew=0.6)
    d.set_ylim(0.6,1.0); d.set_xlim(-0.6,3.6); d.set_title(t,fontsize=4.5,pad=1.5,linespacing=0.95)
    d.set_xticks(range(4)); d.set_xticklabels([SH[m].replace('\n',' ') for m in M],rotation=45,ha='right',fontsize=4.1)
    d.grid(axis='y',lw=0.25,color='#dddddd'); d.set_yticks([0.6,0.8,1.0])
    if j: d.set_yticklabels([])

# ---------------- e per-cell areas vs published
ye=yb+38; letter(88,ye,'e'); mt=ag[ag.matched]
fig.text(96/W,1-(ye+1.2)/H,'Per-cell area, agent vs published\n(number = concordance, CCC)',fontsize=4.8,va='top',linespacing=1.0)
for r,(pk,rk,lab,lim) in enumerate([('pred_cell_px','ref_cell_px','cell area',(0,8000)),('pred_nuc_px','ref_nuc_px','nucleus area',(0,4000))]):
    for mi,m in enumerate(M):
        e=ax_mm(99+mi*9,ye+9+r*11,7.6,7.6); s=mt[mt.method==m]; x,y=s[rk].values,s[pk].values
        e.hexbin(x,y,gridsize=20,extent=lim*2,cmap=mpl.colors.LinearSegmentedColormap.from_list('c',['#ffffff',COL[m]]),bins='log',mincnt=1,lw=0)
        e.plot(lim,lim,'k--',lw=0.3); ccc=2*np.cov(x,y)[0,1]/(x.var()+y.var()+(x.mean()-y.mean())**2)
        e.text(0.03,0.97,f'{ccc:.2f}',transform=e.transAxes,va='top',fontsize=4)
        e.set_xlim(lim); e.set_ylim(lim); e.set_xticks(lim); e.set_yticks(lim); e.tick_params(labelsize=3.6)
        e.set_xticklabels(['0',f'{lim[1]//1000}k']); e.set_yticklabels(['0',f'{lim[1]//1000}k'] if mi==0 else [])
        if r==0: e.set_title(SH[m],fontsize=4,color=COL[m],pad=1,linespacing=0.9)
        if mi==0: e.set_ylabel(f'{lab}\nagent (px²)',fontsize=4.2)
        if r==1 and mi==1: e.set_xlabel('published (px²)',fontsize=4.2,x=1.1)

# ---------------- f donor-level differences
def pct(pos,neg): return 100*(pos-neg)/neg
def slide_avg(df,col,cond='condition'):
    g2=df.groupby(['donor','slide',cond])[col].median().unstack(cond)
    pos=[k for k in g2.columns if 'pos' in str(k).lower()][0]; neg=[k for k in g2.columns if 'neg' in str(k).lower()][0]
    return pct(g2[pos],g2[neg])
auto=pd.read_csv('/auto_stats/donor_contrasts.csv'); auto['d']=pct(auto.donor_pos_summary,auto.donor_neg_summary)
hil=pd.read_csv('/hil_stats/Donor_Level_Metrics.csv'); hil=hil[hil.session_subset=='all']
hp=hil.pivot_table(index=['metric','donor'],columns='cd177_status',values='donor_median_of_image_medians').reset_index(); hp['d']=pct(hp.pos,hp.neg)
cx=pd.read_csv('/codex_cells/cd177_cell_features.csv'); cx=cx[~cx.magnification.astype(str).str.contains('60')]
fc=pd.read_csv('/fast_cells/cell_features.csv'); fc=fc[~fc.image.astype(str).str.contains('60x')]
FEATS=[('Cell area',{'automated':'cell_area','hil':'cell_area_um2','codex':'cell_area','fast':'cell_area'},'3.5×10$^{-52}$, r = 0.27'),
       ('Nucleus area',{'automated':'nucleus_area','hil':'nucleus_area_um2','codex':'nucleus_area','fast':'nucleus_area'},'0.2, r = 0.02'),
       ('Nucleus / cell',{'automated':'nuclear_fraction_of_cell','hil':'nucleus_cell_area_ratio','codex':'nucleus_cell_area_ratio','fast':'nuclear_fraction'},'2.4×10$^{-11}$, r = 0.12')]
qa=pd.read_csv('/auto_stats/cd177_feature_statistics.csv').set_index('feature').p_adjusted
qh=pd.read_csv('/hil_stats/Paired_Donor_Results.csv').drop_duplicates('metric').set_index('metric').p_fdr_bh
qc=pd.read_csv('/codex_stats/cd177_feature_statistics.csv').set_index('feature').q_value_bh
qf=pd.read_csv('/fast_stats/feature_statistics.csv').set_index('feature').p_value_fdr
letter(134,ye,'f'); fig.text(141/W,1-(ye+1.2)/H,'Donor-level CD177$^+$ vs CD177$^-$\n(numbers = FDR q)',fontsize=4.8,va='top',linespacing=1.0)
rng=np.random.default_rng(0)
for fi,(name,fm,pub) in enumerate(FEATS):
    f=ax_mm(143+fi*13.6,ye+11,12.4,14.5)
    for mi,m in enumerate(M):
        if m=='automated': v=auto[auto.feature==fm[m]].d.values; q=qa[fm[m]]
        elif m=='hil': v=hp[hp.metric==fm[m]].d.values; q=qh[fm[m]]
        elif m=='fast': v=slide_avg(fc,fm[m]).groupby('donor').mean().dropna().values; q=qf.get(fm[m],np.nan)
        else: v=slide_avg(cx,fm[m]).dropna().values; q=qc[fm[m]]
        f.scatter(mi+rng.uniform(-.17,.17,len(v)),np.clip(v,-60,150),s=1.8,color=COL[m],alpha=.85,lw=0)
        f.plot([mi-.3,mi+.3],[np.median(v)]*2,color='k',lw=0.6)
        f.text(mi,1.01,(f'{q:.2f}'.lstrip('0') if q<1 else '1.0') if q==q else 'n/a',transform=f.get_xaxis_transform(),ha='center',va='bottom',fontsize=3.3,
               color='#c0392b' if q<0.05 else 'k')
    f.axhline(0,color='#999999',lw=0.3,ls='--'); f.set_xlim(-.6,3.6); f.set_ylim(-60,150); f.set_yticks([-50,0,50,100,150])
    f.set_xticks(range(4)); f.set_xticklabels([SH[m].replace('\n',' ') for m in M],rotation=45,ha='right',fontsize=3.9)
    f.set_title(name,fontsize=4.7,pad=5)
    if fi==0: f.set_ylabel('CD177$^+$ vs CD177$^-$\n(% difference per donor)',fontsize=4.4)
    else: f.set_yticklabels([])
fig.text(141/W,1-(ye+36.3)/H,'Published (cell-level): cell area p = 3.5×10$^{-52}$, r = 0.27;\nnucleus area p = 0.2, r = 0.02; nucleus/cell p = 2.4×10$^{-11}$,\nr = 0.12. Red: q < 0.05.',fontsize=3.6,color=COL['published'],va='top',linespacing=1.1)
for ext in ['pdf','svg','png']: fig.savefig(f'/work/Fig_CD177_halfpage_v4.{ext}',dpi=450 if ext=='png' else None)
