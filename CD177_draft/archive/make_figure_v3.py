"""Half-page CD177 validation figure, v3. Inputs: run outputs, collaborator gallery PNG,
mask_agreement_per_cell.csv (agent vs published masks), nuclei_vs_expert.csv (all vs expert nuclei)."""
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],'font.size':5,
    'axes.linewidth':0.4,'xtick.major.width':0.4,'ytick.major.width':0.4,'xtick.major.size':1.5,'ytick.major.size':1.5,
    'xtick.labelsize':4.5,'ytick.labelsize':4.5,'axes.labelsize':5,'axes.spines.top':False,'axes.spines.right':False,
    'svg.fonttype':'none','pdf.fonttype':42,'xtick.major.pad':1.2,'ytick.major.pad':1.2,'axes.labelpad':1.5})
MM=1/25.4; W,H=183,122
fig=plt.figure(figsize=(W*MM,H*MM))
def ax_mm(x,y,w,h,**k): return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H],**k)
def letter(x,y,s): fig.text(x/W,1-y/H,s,fontsize=8,fontweight='bold',va='top')
COL={'published':'#8e3fb0','automated':'#2e9e45','hil':'#1f9fb0','codex':'#6e6e6e','hil_stock':'#1f9fb0'}
FILL={'manual':'#e5c4f5','hf':'#bfeff3','auto':'#c6f6cf','other':'#e6e6e6'}
LAB={'published':'Published','automated':'Agent only','hil':'Agent + human','codex':'Codex CLI','hil_stock':'micro-sam before fine-tuning'}
M=['automated','hil','codex']

# ---------------- a: pipelines, time, verdict
letter(0,0,'a'); a=ax_mm(0,0,183,42); a.axis('off'); a.set_xlim(0,183); a.set_ylim(42,0)
steps=['Cell\nsegmentation','Nucleus\nsegmentation','QC &\npost-processing','Features','Statistics']
rows=[('Published','published',[('Cellpose (default)','manual'),('ilastik pixel\nclassifier','manual'),('dilation, crop,\nsize filter','manual'),('25 hand-picked','manual'),('cell-level\nMann–Whitney','manual')]),
      ('Agent + human','hil',[('micro-sam vit_b_lm\n(model set by expert)','hf'),(('15 tiles corrected\nin agent-built tool','hf'),('model retrained by expert','manual')),('debris & border QC\n(expert-specified)','hf'),('split-invariant metrics\n(expert-specified)','hf'),('paired Wilcoxon\n(expert-specified)','hf')]),
      ('Agent only','automated',[('Cellpose cyto3\non optical density','auto'),('threshold\nwithin cell','auto'),('size/shape QC','auto'),('32 features','auto'),('donor-level\npermutation + FDR','auto')]),
      ('Codex CLI','codex',[('Cellpose cyto3','other'),('erosion/threshold\nwithin cell','other'),('shape QC','other'),('30 features','other'),('donor-level\nWilcoxon + FDR','other')])]
TIME={'published':(3*7*24*60,'weeks'),'hil':(345,'5.8 h*'),'automated':(67,'67 min'),'codex':(13,'13 min')}
VERD={'published':'no meaningful difference\n(small effect sizes)','hil':'0/11 features significant\n(FDR); cell-area trend','automated':'0/32 features significant\n(FDR); cell-size trend','codex':'0/30 features significant\n(FDR)'}
x0,cw,gap,rh,y0=17,18.4,0.5,6.8,6
for i,s_ in enumerate(steps): a.text(x0+i*(cw+gap)+cw/2,3,s_,ha='center',va='center',fontsize=5,fontweight='bold',linespacing=1.0)
a.text(126,3,'Time to result',ha='center',va='center',fontsize=5,fontweight='bold'); a.text(163,3,'Scientific verdict',ha='center',va='center',fontsize=5,fontweight='bold')
for r,(name,key,cells) in enumerate(rows):
    y=y0+r*(rh+0.7)
    a.text(0.5,y+rh/2,name,va='center',fontsize=5.5,fontweight='bold',color=COL[key])
    for i,c in enumerate(cells):
        xx=x0+i*(cw+gap)
        if isinstance(c[0],tuple):
            for j,(tt,kind) in enumerate(c):
                a.add_patch(Rectangle((xx,y+j*rh/2),cw,rh/2,fc=FILL[kind],ec='k',lw=0.35))
                a.text(xx+cw/2,y+j*rh/2+rh/4,tt,ha='center',va='center',fontsize=4.1,linespacing=0.9)
        else:
            a.add_patch(Rectangle((xx,y),cw,rh,fc=FILL[c[1]],ec='k',lw=0.35))
            a.text(xx+cw/2,y+rh/2,c[0],ha='center',va='center',fontsize=4.5,linespacing=1.0)
    a.text(149,y+rh/2,'\u2713',ha='center',va='center',fontsize=7,color=COL[key],fontweight='bold',fontfamily='DejaVu Sans')
    a.text(152,y+rh/2,VERD[key],va='center',fontsize=4.5,linespacing=1.05)
ta=ax_mm(112,y0,28,4*(rh+0.7)-0.7)
for r,(name,key,_) in enumerate(rows):
    v,lab=TIME[key]; ta.barh(r,v,height=0.55,color=COL[key],lw=0); ta.text(v*1.25,r,lab,va='center',fontsize=4.5)
ta.set_xscale('log'); ta.set_xlim(5,3e5); ta.set_ylim(3.5,-0.5); ta.set_yticks([]); ta.minorticks_off()
ta.set_xticks([10,60,1440,10080]); ta.set_xticklabels(['10 min','1 h','1 day','1 week'],fontsize=4.2); ta.spines['left'].set_visible(False)
for i,(k,t) in enumerate([('manual','Manual'),('hf','Human feedback'),('auto','Automatic (Agentic-J)'),('other','Automatic (other agent)')]):
    a.add_patch(Rectangle((x0+i*23,36.3),2.2,2,fc=FILL[k],ec='k',lw=0.3)); a.text(x0+i*23+3,37.3,t,va='center',fontsize=4.8)
a.text(x0,40.6,'All agent runs started from an expert-written task description. *wall-clock incl. expert time.',fontsize=4.4,color='#555555',va='center')

# ---------------- b: gallery, all four rows
letter(0,44,'b'); g=plt.imread('/work/gallery.png'); sc=g.shape[1]/2000
cols=[212,512,813,1114,1414,1715]; rws=[100,421,743,1064]; T=264
clab=['raw','Published','Agent only\n(advanced)','Agent only\n(fast)','Agent +\nhuman','Codex CLI']
ccol=['k',COL['published'],COL['automated'],COL['automated'],COL['hil'],COL['codex']]
rlab=['V8 neg\nBottom','V10 pos\nBottom','V8 neg\nTop','V9 pos\nTop']
notes={(0,5):'missed',(1,1):'debris as cells',(2,1):'nucleus only',(2,2):'missed',(3,3):'wrong split'}
tw,tg,rp=12.4,0.45,14.4
for ri,yy in enumerate(rws):
    for ci,xx in enumerate(cols):
        t=g[int(yy*sc):int((yy+T)*sc),int(xx*sc):int((xx+T)*sc)]
        ax=ax_mm(8+ci*(tw+tg),49.5+ri*rp,tw,tw); ax.imshow(t); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(0.3)
        if ri==0: ax.set_title(clab[ci],fontsize=4.5,color=ccol[ci],pad=1.2,linespacing=0.95)
        if ci==0: ax.set_ylabel(rlab[ri],fontsize=4.4,rotation=0,ha='right',va='center',labelpad=1.5)
        if (ri,ci) in notes: ax.text(0.5,-0.02,notes[(ri,ci)],transform=ax.transAxes,ha='center',va='top',fontsize=4,color='#b0172a')
lx=ax_mm(8,108.5,80,3); lx.axis('off'); lx.set_xlim(0,80); lx.set_ylim(0,1)
for x,c,ls,t in [(8,'#f2c200','-','cell'),(22,'#27c7f0','-','nucleus'),(38,'#999999','--','segmented, not measured')]:
    lx.plot([x,x+4],[0.5,0.5],color=c,ls=ls,lw=1); lx.text(x+5,0.5,t,va='center',fontsize=4.6)

# ---------------- c: nucleus accuracy against the expert's annotations
letter(88,44,'c'); nx=pd.read_csv('/work/nuclei_vs_expert.csv')
th=np.arange(0.5,0.85,0.05); c=ax_mm(96,51,29,22)
for m,ls,lw in [('hil','-',1.1),('codex','-',0.8),('published','-',1.1),('automated','-',0.8),('hil_stock',':',0.8)]:
    s=nx[nx.method==m]; c.plot(th,[s[f'acc_{t:.2f}'].mean() for t in th],ls,color=COL[m],lw=lw,label=f"{LAB[m].replace('micro-sam before fine-tuning','before fine-tuning')}  {s.msa.mean():.2f}")
c.set_xlim(0.5,0.8); c.set_ylim(0,0.6); c.set_xticks([0.5,0.6,0.7,0.8]); c.set_yticks([0,0.2,0.4,0.6])
c.set_xlabel('IoU threshold (stricter →)'); c.set_ylabel('Nucleus accuracy')
n_nuc=int(nx[nx.method=='hil'].n_true.sum()); n_t=nx.tile.nunique()
fig.text(96/W,1-45.2/H,f'Nuclei vs expert annotation\n({n_nuc} nuclei, {n_t} held-out tiles)',fontsize=4.8,va='top',linespacing=1.0)
c.legend(frameon=False,fontsize=3.9,loc='upper right',bbox_to_anchor=(1.02,1.04),handlelength=1.0,borderaxespad=0,labelspacing=0.2,title='mSA',title_fontsize=3.9,alignment='right')

# ---------------- d: agreement with the published masks (clearer version of the first draft)
letter(128,44,'d'); ag=pd.read_csv('/work/mask_agreement_per_cell.csv'); ag['batch']=np.where(ag.image.str.contains('Top'),'Top','Bottom'); ag['matched']=ag.cell_iou>=0.5
fig.text(136/W,1-45.2/H,f'Agreement with published masks\n(n = {int((ag.method=="hil").sum())} cells)',fontsize=4.8,va='top',linespacing=1.0)
mets=[('matched','Published cells\nfound (IoU ≥ 0.5)'),('cell_iou','Cell overlap\n(IoU)'),('nuc_dice','Nucleus overlap\n(Dice)')]
for j,(k,t) in enumerate(mets):
    d=ax_mm(137+j*16,54,13,18)
    for mi,m in enumerate(M):
        for bi,(bt,mk) in enumerate([('Bottom','o'),('Top','s')]):
            s=ag[(ag.method==m)&(ag.batch==bt)]; xp=mi+(bi-0.5)*0.36
            if k=='matched': val=s.matched.mean()
            else:
                q=s[k][s.matched]; val=q.median(); d.plot([xp,xp],q.quantile([.25,.75]),color=COL[m],lw=0.6)
            d.plot(xp,val,mk,color=COL[m],ms=2.6,mfc=COL[m] if bt=='Bottom' else 'white',mew=0.6)
    d.set_ylim(0.6,1.0); d.set_xlim(-0.6,2.6); d.set_title(t,fontsize=4.6,pad=1.5,linespacing=0.95)
    d.set_xticks(range(3)); d.set_xticklabels([LAB[m] for m in M],rotation=40,ha='right',fontsize=4.3)
    d.grid(axis='y',lw=0.25,color='#dddddd'); d.set_yticks([0.6,0.8,1.0])
    if j: d.set_yticklabels([])
d.legend(handles=[Line2D([],[],marker='o',ls='',color='k',ms=2.6,label='Bottom slides'),Line2D([],[],marker='s',ls='',color='k',mfc='white',ms=2.6,label='Top slides')],
         frameon=False,fontsize=4.2,loc='lower right',handletextpad=0.1,borderaxespad=0,labelspacing=0.2)

# ---------------- e: per-cell areas vs published (matched cells)
letter(88,82,'e'); mt=ag[ag.matched]
for r,(pk,rk,lab,lim) in enumerate([('pred_cell_px','ref_cell_px','cell area',(0,8000)),('pred_nuc_px','ref_nuc_px','nucleus area',(0,4000))]):
    for mi,m in enumerate(M):
        e=ax_mm(99+mi*10.5,90+r*11.8,8.6,8.6); s=mt[mt.method==m]; x,y=s[rk].values,s[pk].values
        e.hexbin(x,y,gridsize=22,extent=lim*2,cmap=mpl.colors.LinearSegmentedColormap.from_list('c',['#ffffff',COL[m]]),bins='log',mincnt=1,lw=0)
        e.plot(lim,lim,'k--',lw=0.3); ccc=2*np.cov(x,y)[0,1]/(x.var()+y.var()+(x.mean()-y.mean())**2)
        e.text(0.03,0.97,f'{ccc:.2f}',transform=e.transAxes,va='top',fontsize=4.2)
        e.set_xlim(lim); e.set_ylim(lim); e.set_xticks(lim); e.set_yticks(lim); e.tick_params(labelsize=3.8)
        e.set_xticklabels(['0',f'{lim[1]//1000}k']); e.set_yticklabels(['0',f'{lim[1]//1000}k'] if mi==0 else [])
        if r==0: e.set_title(LAB[m].replace(' + ','+\n') if m=='hil' else LAB[m].replace(' ','\n'),fontsize=4.2,color=COL[m],pad=1,linespacing=0.9)
        if mi==0: e.set_ylabel(f'{lab}\nagent (px²)',fontsize=4.3)
        if r==1 and mi==1: e.set_xlabel('published (px²)',fontsize=4.3)
fig.text(94/W,1-82.6/H,'Per-cell area, agent vs published\n(number = concordance, CCC)',fontsize=4.8,va='top',linespacing=1.0)

# ---------------- f: donor-level CD177+ vs CD177- (labelling of the first draft)
def pct(pos,neg): return 100*(pos-neg)/neg
auto=pd.read_csv('/auto_stats/donor_contrasts.csv'); auto['d']=pct(auto.donor_pos_summary,auto.donor_neg_summary)
hil=pd.read_csv('/hil_stats/Donor_Level_Metrics.csv'); hil=hil[hil.session_subset=='all']
hp=hil.pivot_table(index=['metric','donor'],columns='cd177_status',values='donor_median_of_image_medians').reset_index(); hp['d']=pct(hp.pos,hp.neg)
cx=pd.read_csv('/codex_cells/cd177_cell_features.csv'); cx=cx[~cx.magnification.astype(str).str.contains('60')]
FEATS=[('Cell area',{'automated':'cell_area','hil':'cell_area_um2','codex':'cell_area'},'p = 3.5×10$^{-52}$\nr = 0.27'),
       ('Nucleus area',{'automated':'nucleus_area','hil':'nucleus_area_um2','codex':'nucleus_area'},'p = 0.2\nr = 0.02'),
       ('Nucleus / cell',{'automated':'nuclear_fraction_of_cell','hil':'nucleus_cell_area_ratio','codex':'nucleus_cell_area_ratio'},'p = 2.4×10$^{-11}$\nr = 0.12')]
qa=pd.read_csv('/auto_stats/cd177_feature_statistics.csv').set_index('feature').p_adjusted
qh=pd.read_csv('/hil_stats/Paired_Donor_Results.csv').drop_duplicates('metric').set_index('metric').p_fdr_bh
qc=pd.read_csv('/codex_stats/cd177_feature_statistics.csv').set_index('feature').q_value_bh
letter(134,82,'f'); fig.text(141/W,1-82.6/H,'Donor-level CD177$^+$ vs CD177$^-$\n(numbers = FDR q)',fontsize=4.8,va='top'); rng=np.random.default_rng(0)
for fi,(name,fm,pub) in enumerate(FEATS):
    f=ax_mm(143+fi*13.7,92,12,14.5)
    for mi,m in enumerate(M):
        if m=='automated': v=auto[auto.feature==fm[m]].d.values; q=qa[fm[m]]
        elif m=='hil': v=hp[hp.metric==fm[m]].d.values; q=qh[fm[m]]
        else:
            g2=cx.groupby(['donor','slide','condition'])[fm[m]].median().unstack('condition')
            pos=[k for k in g2.columns if 'pos' in str(k).lower()][0]; neg=[k for k in g2.columns if 'neg' in str(k).lower()][0]
            v=pct(g2[pos],g2[neg]).dropna().values; q=qc[fm[m]]
        f.scatter(mi+rng.uniform(-.17,.17,len(v)),np.clip(v,-60,150),s=2.2,color=COL[m],alpha=.85,lw=0)
        f.plot([mi-.3,mi+.3],[np.median(v)]*2,color='k',lw=0.6)
        f.text(mi,1.01,f'{q:.2f}'.lstrip('0'),transform=f.get_xaxis_transform(),ha='center',va='bottom',fontsize=3.7)
    f.axhline(0,color='#999999',lw=0.3,ls='--'); f.set_xlim(-.6,2.6); f.set_ylim(-60,150); f.set_yticks([-50,0,50,100,150])
    f.set_xticks(range(3)); f.set_xticklabels([LAB[m] for m in M],rotation=40,ha='right',fontsize=4.1)
    f.set_title(name,fontsize=4.8,pad=5)
    f.text(0.5,-0.62,'published:\n'+pub.replace('\n',', '),transform=f.transAxes,ha='center',va='top',fontsize=3.5,color=COL['published'],linespacing=1.0)
    if fi==0:
        f.set_ylabel('CD177$^+$ vs CD177$^-$\n(% difference per donor)',fontsize=4.5)
    else: f.set_yticklabels([])
for ext in ['pdf','svg','png']: fig.savefig(f'/work/Fig_CD177_halfpage_v3.{ext}',dpi=450 if ext=='png' else None)
