"""CD177 validation figure v7: Nature type sizes (min 5 pt, labels 6 pt, panel letters 8 pt bold).
Footnote text lives in the caption, not the figure."""
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
FS,FL=5,6          # minimum text size / axis-title size (pt)
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],'font.size':FS,
    'axes.linewidth':0.5,'xtick.major.width':0.5,'ytick.major.width':0.5,'xtick.major.size':2,'ytick.major.size':2,
    'xtick.labelsize':FS,'ytick.labelsize':FS,'axes.labelsize':FS,'axes.titlesize':FS,'axes.spines.top':False,'axes.spines.right':False,
    'svg.fonttype':'none','pdf.fonttype':42,'xtick.major.pad':1.2,'ytick.major.pad':1.2,'axes.labelpad':1.5})
MM=1/25.4; W,H=183,134
fig=plt.figure(figsize=(W*MM,H*MM))
def ax_mm(x,y,w,h,**k): return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H],**k)
def letter(x,y,s): fig.text(x/W,1-y/H,s,fontsize=8,fontweight='bold',va='top')
def head(x,y,s): fig.text(x/W,1-y/H,s,fontsize=FL,va='top')
COL={'published':'#7b3294','hil':'#0571b0','automated':'#1a9641','fast':'#a6611a','codex':'#636363'}
FILL={'manual':'#f7c873','hf':'#8ecdf2','ai':'#f4a3c0','auto':'#a6dba0','other':'#d9d9d9'}
SH={'automated':'Advanced','fast':'Fast','hil':'Agent + human','codex':'Codex'}
M=['automated','fast','hil','codex']

# ---------------- a pipelines / time / verdict
letter(0,0,'a'); a=ax_mm(0,0,183,48); a.axis('off'); a.set_xlim(0,183); a.set_ylim(48,0)
steps=['Cell\nsegmentation','Nucleus\nsegmentation','QC &\npost-processing','Features','Statistics']
rows=[('Published','published',[('Cellpose (default)','manual'),('ilastik pixel\nclassifier','manual'),('dilation, crop,\nsize filter','manual'),('25 hand-picked','manual'),('cell-level\nMann–Whitney','manual')]),
      ('Agent +\nhuman','hil',[('micro-sam vit_b_lm\n(AI assistant)','ai'),(('user fixed 15 tiles','hf'),('AI assistant retrained','ai')),('debris/border QC\n(AI assistant)','ai'),('lobe metrics\n(AI assistant)','ai'),('paired Wilcoxon\n(AI assistant)','ai')]),
      ('Advanced','automated',[('Cellpose cyto3\non optical density','auto'),('threshold\nwithin cell','auto'),('size/shape QC','auto'),('32 features','auto'),('donor-level\npermutation, FDR','auto')]),
      ('Fast','fast',[('Cellpose cyto3','auto'),('Otsu threshold\nwithin cell','auto'),('border, size,\nsolidity filter','auto'),('25 features','auto'),('donor-level\nWilcoxon, FDR','auto')]),
      ('Codex CLI','codex',[('Cellpose cyto3\n(half resolution)','other'),('erosion/threshold\nwithin cell','other'),('shape QC','other'),('30 features','other'),('donor-level\nWilcoxon, FDR','other')])]
TIME={'published':(3*7*24*60,'weeks'),'hil':(345,'5.8 h'),'automated':(67,'67 min'),'fast':(35,'35 min'),'codex':(13,'13 min')}
VERD={'published':('ok','no meaningful difference'),'hil':('ok','0/11 features significant'),
      'automated':('ok','0/32 significant (3/3 runs)'),'fast':('x','5/25 significant (3/3 runs)'),'codex':('ok','0/30 features significant')}
x0,cw,gap,rh,y0=13,19.6,0.5,6.6,6.2
for i,s_ in enumerate(steps): a.text(x0+i*(cw+gap)+cw/2,3.1,s_,ha='center',va='center',fontsize=FS,fontweight='bold',linespacing=1.0)
a.text(129,3.1,'Time to result',ha='center',va='center',fontsize=FS,fontweight='bold'); a.text(163,3.1,'Scientific verdict',ha='center',va='center',fontsize=FS,fontweight='bold')
for r,(name,key,cells) in enumerate(rows):
    y=y0+r*(rh+0.6)
    a.text(0.3,y+rh/2,name,va='center',fontsize=FS,fontweight='bold',color=COL[key],linespacing=1.0)
    for i,c in enumerate(cells):
        xx=x0+i*(cw+gap)
        if isinstance(c[0],tuple):
            for j,(tt,kind) in enumerate(c):
                a.add_patch(Rectangle((xx,y+j*rh/2),cw,rh/2,fc=FILL[kind],ec='k',lw=0.35))
                a.text(xx+cw/2,y+j*rh/2+rh/4,tt,ha='center',va='center',fontsize=FS)
        else:
            a.add_patch(Rectangle((xx,y),cw,rh,fc=FILL[c[1]],ec='k',lw=0.35))
            a.text(xx+cw/2,y+rh/2,c[0],ha='center',va='center',fontsize=FS,linespacing=1.0)
    ok,txt=VERD[key]
    a.text(146.5,y+rh/2,'\u2713' if ok=='ok' else '\u2717',ha='center',va='center',fontsize=7,fontweight='bold',
           color=COL[key] if ok=='ok' else '#c0392b',fontfamily='DejaVu Sans')
    a.text(149,y+rh/2,txt,va='center',fontsize=FS)
ta=ax_mm(114,y0,29,5*(rh+0.6)-0.6)
for r,(name,key,_) in enumerate(rows):
    v,lab=TIME[key]; ta.barh(r,v,height=0.55,color=COL[key],lw=0); ta.text(v*1.3,r,lab,va='center',fontsize=FS)
ta.set_xscale('log'); ta.set_xlim(5,4e5); ta.set_ylim(4.5,-0.5); ta.set_yticks([]); ta.minorticks_off()
ta.set_xticks([10,60,1440,10080]); ta.set_xticklabels(['10 min','1 h','1 d','1 wk']); ta.spines['left'].set_visible(False)
xl=x0
for k,t in [('manual','Manual (analyst)'),('hf','Human feedback (user)'),('ai','AI assistant (Claude Code)'),('auto','Automatic (Agentic-J)'),('other','Automatic (other agent)')]:
    a.add_patch(Rectangle((xl,44.8),2.4,2.2,fc=FILL[k],ec='k',lw=0.3)); a.text(xl+3.2,45.9,t,va='center',fontsize=FS); xl+=3.2+len(t)*0.98+3.5

# ---------------- b gallery (all four rows)
yb=51; letter(0,yb,'b'); g=plt.imread('/work/gallery.png'); sc=g.shape[1]/2000
cols=[212,512,813,1114,1414,1715]; rws=[100,421,743,1064]; T=264
clab=['raw','Published','Advanced','Fast','Agent +\nhuman','Codex']
ccol=['k',COL['published'],COL['automated'],COL['fast'],COL['hil'],COL['codex']]
rlab=['V8 neg\nBottom','V10 pos\nBottom','V8 neg\nTop','V9 pos\nTop']
notes={(0,5):'missed',(1,1):'debris',(2,1):'nucleus only',(2,2):'missed',(3,3):'wrong split'}
tw,tg,rp=12.5,0.45,15.2
for ri,yy in enumerate(rws):
    for ci,xx in enumerate(cols):
        t=g[int(yy*sc):int((yy+T)*sc),int(xx*sc):int((xx+T)*sc)]
        ax=ax_mm(9+ci*(tw+tg),yb+6+ri*rp,tw,tw); ax.imshow(t); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(0.3)
        if ri==0: ax.set_title(clab[ci],fontsize=FS,color=ccol[ci],pad=1.5,linespacing=0.95)
        if ci==0: ax.set_ylabel(rlab[ri],fontsize=FS,rotation=0,ha='right',va='center',labelpad=1.5)
        if (ri,ci) in notes: ax.text(0.5,-0.03,notes[(ri,ci)],transform=ax.transAxes,ha='center',va='top',fontsize=FS,color='#b0172a')
lx=ax_mm(9,yb+68,78,3); lx.axis('off'); lx.set_xlim(0,78); lx.set_ylim(0,1)
for x,c,ls,t in [(2,'#f2c200','-','cell'),(15,'#27c7f0','-','nucleus'),(31,'#999999','--','segmented, not measured')]:
    lx.plot([x,x+4],[0.5,0.5],color=c,ls=ls,lw=1); lx.text(x+5,0.5,t,va='center',fontsize=FS)

# ---------------- c agreement with published masks
letter(90,yb,'c'); ag=pd.read_csv('/work/mask_agreement_per_cell.csv'); ag['matched']=ag.cell_iou>=0.5
head(97,yb+0.6,f'Agreement with the published masks ({int((ag.method=="hil").sum())} cells)')
mets=[('matched','Published cells\nfound (IoU ≥ 0.5)'),('cell_iou','Cell overlap\n(IoU)'),('nuc_dice','Nucleus overlap\n(Dice)')]
for j,(k,t) in enumerate(mets):
    d=ax_mm(100+j*28,yb+9,24,14)
    for mi,m in enumerate(M):
        s_=ag[ag.method==m]
        if k=='matched': d.bar(mi,s_.matched.mean(),0.62,color=COL[m],lw=0)
        else:
            q=s_[k][s_.matched]; d.plot([mi,mi],q.quantile([.25,.75]),color=COL[m],lw=1.0); d.plot(mi,q.median(),'o',color=COL[m],ms=3)
    d.set_ylim(0.6,1.0); d.set_xlim(-0.6,3.6); d.set_title(t,fontsize=FS,pad=2,linespacing=0.95)
    d.set_xticks(range(4)); d.set_xticklabels([SH[m] for m in M],rotation=40,ha='right')
    d.grid(axis='y',lw=0.25,color='#dddddd'); d.set_axisbelow(True); d.set_yticks([0.6,0.8,1.0])
    if j: d.set_yticklabels([])

# ---------------- d histograms vs published (Fig. S2c style)
qa=pd.read_csv('/auto_stats/cd177_feature_statistics.csv').set_index('feature').p_adjusted
qh=pd.read_csv('/hil_stats/Paired_Donor_Results.csv').drop_duplicates('metric').set_index('metric').p_fdr_bh
qc=pd.read_csv('/codex_stats/cd177_feature_statistics.csv').set_index('feature').q_value_bh
qf=pd.read_csv('/fast_stats/feature_statistics.csv').set_index('feature').p_value_fdr
ye=yb+37.5; letter(90,ye,'d'); head(97,ye+0.6,'Per-cell distributions vs the published Fig. S2c')
PX2=0.2515**2
ref=ag[ag.method=='hil'].copy(); ref['cond']=np.where(ref.image.str.contains('pos'),'pos','neg')
ref=pd.DataFrame({'cond':ref.cond,'cell':ref.ref_cell_px*PX2,'nuc':ref.ref_nuc_px*PX2,'nc':ref.ref_nuc_px/ref.ref_cell_px})
A=pd.read_csv('/auto_cells/cd177_cell_features.csv'); A=A[A.size_analysis_eligible.astype(str)=='True']
F=pd.read_csv('/fast_cells/cell_features.csv'); F=F[~F.image.astype(str).str.contains('60x')]
Cx=pd.read_csv('/codex_cells/cd177_cell_features.csv'); Cx=Cx[~Cx.magnification.astype(str).str.contains('60')]
Hh=pd.read_csv('/hil_cells/per_cell_morphology_eligible.csv'); Hh=Hh[(Hh.eligible_cell_metrics.astype(str)=='True')&(Hh.zero_lobes.astype(str)=='False')]
def std(df,cond,c,n,r): return pd.DataFrame({'cond':np.where(df[cond].astype(str).str.lower().str.contains('pos'),'pos','neg'),'cell':df[c],'nuc':df[n],'nc':df[r]})
DATA={'published':ref,'automated':std(A,'condition','cell_area','nucleus_area','nuclear_fraction_of_cell'),
      'fast':std(F,'condition','cell_area','nucleus_area','nuclear_fraction'),'hil':std(Hh,'cd177_status','cell_area_um2','nucleus_area_um2','nucleus_cell_area_ratio'),
      'codex':std(Cx,'condition','cell_area','nucleus_area','nucleus_cell_area_ratio')}
Q={'automated':(qa,('cell_area','nucleus_area','nuclear_fraction_of_cell')),'fast':(qf,('cell_area','nucleus_area','nuclear_fraction')),
   'hil':(qh,('cell_area_um2','nucleus_area_um2','nucleus_cell_area_ratio')),'codex':(qc,('cell_area','nucleus_area','nucleus_cell_area_ratio'))}
PUB=('p = 3×10$^{-52}$','p = 0.2','p = 2×10$^{-11}$')
FEAT=[('cell','Cell area (µm²)',(0,450)),('nuc','Nucleus area (µm²)',(0,160)),('nc','Nucleus / cell area',(0,1))]
RM=['published','automated','fast','hil','codex']; RL={'published':'Published','automated':'Advanced','fast':'Fast','hil':'Agent +\nhuman','codex':'Codex'}
NEG,POS='#f08080','#7b68ee'
for ri,m in enumerate(RM):
    for fj,(k,t,lim) in enumerate(FEAT):
        h=ax_mm(104+fj*26.5,ye+7+ri*6.6,24.5,5.4); v=DATA[m]; bins=np.linspace(*lim,46)
        for cnd,colr in [('neg',NEG),('pos',POS)]:
            h.hist(v[v.cond==cnd][k].dropna().clip(*lim),bins=bins,density=True,color=colr,alpha=0.55,lw=0)
        h.set_xlim(lim); h.set_xticks({'cell':[0,200,400],'nuc':[0,50,100,150],'nc':[0,0.5,1]}[k]); h.set_yticks([]); h.spines['left'].set_visible(False)
        if m=='published': txt,tc=PUB[fj],COL['published']
        else:
            q=Q[m][0].get(Q[m][1][fj],np.nan); txt=f'q = {q:.2f}'; tc='#c0392b' if q<0.05 else '#333333'
        h.text(1.0,1.0,txt,transform=h.transAxes,ha='right',va='top',fontsize=FS,color=tc)
        if ri<len(RM)-1: h.set_xticklabels([])
        else: h.set_xlabel(t)
        if fj==0: h.text(-0.05,0.4,RL[m],transform=h.transAxes,ha='right',va='center',fontsize=FS,fontweight='bold',color=COL[m],linespacing=0.95)
fig.legend(handles=[Patch(fc=NEG,alpha=0.55,label='CD177$^-$'),Patch(fc=POS,alpha=0.55,label='CD177$^+$')],frameon=False,fontsize=FS,ncol=2,
           loc='upper right',bbox_to_anchor=(182.5/W,1-(ye+0.3)/H),handlelength=1.0,columnspacing=0.8,borderaxespad=0)
for ext in ['pdf','svg','png']: fig.savefig(f'/work/Fig_CD177_v7.{ext}',dpi=450 if ext=='png' else None)
