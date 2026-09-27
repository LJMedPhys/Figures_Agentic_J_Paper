"""Half-page draft (183 x 120 mm) of the CD177 validation figure. Reads only existing run outputs,
the collaborator's gallery PNG and mask_agreement_per_cell.csv (overlap of existing masks)."""
import json, numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],'font.size':5.5,
    'axes.linewidth':0.4,'xtick.major.width':0.4,'ytick.major.width':0.4,'xtick.major.size':1.5,'ytick.major.size':1.5,
    'xtick.labelsize':5,'ytick.labelsize':5,'axes.labelsize':5.5,'axes.spines.top':False,'axes.spines.right':False,
    'svg.fonttype':'none','pdf.fonttype':42})
MM=1/25.4; W,H=183,120
fig=plt.figure(figsize=(W*MM,H*MM))
def ax_mm(x,y,w,h,**k): return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H],**k)
def letter(x,y,s): fig.text(x/W,1-y/H,s,fontsize=8,fontweight='bold',va='top')
COL={'published':'#8e3fb0','automated':'#2e9e45','hil':'#1f9fb0','codex':'#6e6e6e'}
FILL={'manual':'#e5c4f5','hf':'#bfeff3','auto':'#c6f6cf','other':'#e6e6e6'}
M=['automated','hil','codex']; SHORT={'automated':'Agent\nonly','hil':'Agent +\nhuman','codex':'Codex'}

# ---------------- a: pipelines + time + verdict (human involvement from the HIL chat log)
letter(0,0,'a'); a=ax_mm(0,0,183,44); a.axis('off'); a.set_xlim(0,183); a.set_ylim(44,0)
steps=['Cell\nsegmentation','Nucleus\nsegmentation','QC &\npost-processing','Features','Statistics']
rows=[('Published','published',[('Cellpose (default)','manual'),('ilastik pixel\nclassifier','manual'),('dilation, crop,\nsize filter','manual'),('25 hand-picked','manual'),('cell-level\nMann–Whitney','manual')]),
      ('Agent + human','hil',[('micro-sam vit_b_lm\n(model set by expert)','hf'),(('15 tiles corrected\nin agent-built tool','hf'),('model retrained by expert','manual')),('debris & border QC\n(expert-specified)','hf'),('split-invariant metrics\n(expert-specified)','hf'),('paired Wilcoxon\n(expert-specified)','hf')]),
      ('Agent only','automated',[('Cellpose cyto3\non optical density','auto'),('threshold\nwithin cell','auto'),('size/shape QC','auto'),('32 features','auto'),('donor-level\npermutation + FDR','auto')]),
      ('Codex CLI','codex',[('Cellpose cyto3','other'),('erosion/threshold\nwithin cell','other'),('shape QC','other'),('30 features','other'),('donor-level\nWilcoxon + FDR','other')])]
TIME={'published':(3*7*24*60,'weeks'),'hil':(345,'5.8 h*'),'automated':(67,'67 min'),'codex':(13,'13 min')}
VERD={'published':'no meaningful difference\n(small effect sizes)','hil':'0/11 features significant\n(FDR); cell-area trend','automated':'0/32 features significant\n(FDR); cell-size trend','codex':'0/30 features significant\n(FDR)'}
x0,cw,gap,rh,y0=17,18.4,0.5,7.3,6.5
for i,s_ in enumerate(steps): a.text(x0+i*(cw+gap)+cw/2,3.2,s_,ha='center',va='center',fontsize=5,fontweight='bold',linespacing=1.0)
a.text(126,3.2,'Time to result',ha='center',va='center',fontsize=5,fontweight='bold')
a.text(163,3.2,'Scientific verdict',ha='center',va='center',fontsize=5,fontweight='bold')
for r,(name,key,cells) in enumerate(rows):
    y=y0+r*(rh+0.8)
    a.text(0.5,y+rh/2,name,va='center',fontsize=5.5,fontweight='bold',color=COL[key])
    for i,c in enumerate(cells):
        xx=x0+i*(cw+gap)
        if isinstance(c[0],tuple):
            for j,(tt,kind) in enumerate(c):
                a.add_patch(Rectangle((xx,y+j*rh/2),cw,rh/2,fc=FILL[kind],ec='k',lw=0.35))
                a.text(xx+cw/2,y+j*rh/2+rh/4,tt,ha='center',va='center',fontsize=4.2,linespacing=0.95)
        else:
            a.add_patch(Rectangle((xx,y),cw,rh,fc=FILL[c[1]],ec='k',lw=0.35))
            a.text(xx+cw/2,y+rh/2,c[0],ha='center',va='center',fontsize=4.5,linespacing=1.0)
    a.text(149,y+rh/2,'\u2713',ha='center',va='center',fontsize=7,color=COL[key],fontweight='bold',fontfamily='DejaVu Sans')
    a.text(152,y+rh/2,VERD[key],va='center',fontsize=4.5,linespacing=1.05)
# time: log-scale bars on a shared axis aligned with the rows
ta=ax_mm(112,y0,28,4*(rh+0.8)-0.8)
for r,(name,key,_) in enumerate(rows):
    v,lab=TIME[key]; ta.barh(r,v,height=0.55,color=COL[key],lw=0); ta.text(v*1.25,r,lab,va='center',fontsize=4.5)
ta.set_xscale('log'); ta.set_xlim(5,3e5); ta.set_ylim(3.5,-0.5); ta.set_yticks([])
ta.set_xticks([10,60,1440,10080]); ta.set_xticklabels(['10 min','1 h','1 day','1 week'],fontsize=4.2)
ta.minorticks_off(); ta.spines['left'].set_visible(False); ta.tick_params(axis='x',pad=1)
for i,(k,t) in enumerate([('manual','Manual'),('hf','Human feedback'),('auto','Automatic (Agentic-J)'),('other','Automatic (other agent)')]):
    a.add_patch(Rectangle((x0+i*23,40.2),2.2,2,fc=FILL[k],ec='k',lw=0.3)); a.text(x0+i*23+3,41.2,t,va='center',fontsize=4.8)
a.text(x0,43.6,'All agent runs started from an expert-written task description. *wall-clock incl. expert time.',fontsize=4.4,color='#555555',va='center')

# ---------------- b: gallery (tiles cropped from the collaborator's panel)
letter(0,47,'b')
g=plt.imread('/work/gallery.png'); sc=g.shape[1]/2000
cols=[212,512,813,1114,1414,1715]; rws=[421,743,1064]; T=264
clab=['raw','published','Agent only\n(advanced)','Agent only\n(fast)','Agent + human','Codex']
ccol=['k',COL['published'],COL['automated'],COL['automated'],COL['hil'],COL['codex']]
rlab=['V10 pos\nBottom','V8 neg\nTop','V9 pos\nTop']
notes={(0,1):'debris as cells',(1,1):'nucleus only',(1,2):'missed',(2,3):'wrong split'}
tw,tg=14.6,0.6
for ri,yy in enumerate(rws):
    for ci,xx in enumerate(cols):
        t=g[int(yy*sc):int((yy+T)*sc),int(xx*sc):int((xx+T)*sc)]
        ax=ax_mm(8+ci*(tw+tg),52.5+ri*(tw+2.4),tw,tw); ax.imshow(t); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(0.3)
        if ri==0: ax.set_title(clab[ci],fontsize=4.8,color=ccol[ci],pad=1.5,linespacing=1.0)
        if ci==0: ax.set_ylabel(rlab[ri],fontsize=4.6,rotation=0,ha='right',va='center',labelpad=2)
        if (ri,ci) in notes: ax.text(0.5,-0.02,notes[(ri,ci)],transform=ax.transAxes,ha='center',va='top',fontsize=4.2,color='#b0172a')
lx=ax_mm(8,104.5,92,3); lx.axis('off'); lx.set_xlim(0,97); lx.set_ylim(0,1)
for x,c,ls,t in [(20,'#f2c200','-','cell'),(38,'#27c7f0','-','nucleus'),(58,'#999999','--','segmented, not measured')]:
    lx.plot([x,x+4],[0.5,0.5],color=c,ls=ls,lw=1); lx.text(x+5,0.5,t,va='center',fontsize=4.8)

# ---------------- c: fine-tuning effect (expert's 4-fold cross-validation, rescore.json)
letter(106,47,'c'); rs=json.load(open('/rescore/rescore.json'))['vit_b_lm']
c=ax_mm(114,52,17,25); tiles=sorted(rs['stock'])
st=[rs['stock'][t]['nuc_sa50'] for t in tiles]; ft=[rs['best'][t]['nuc_sa50'] for t in tiles]
for s_,f_ in zip(st,ft): c.plot([0,1],[s_,f_],color='#bbbbbb',lw=0.4)
c.scatter(np.zeros(len(st)),st,s=5,color='#999999',zorder=3,lw=0); c.scatter(np.ones(len(ft)),ft,s=5,color=COL['hil'],zorder=3,lw=0)
for x_,v in [(0,st),(1,ft)]: c.plot([x_-.2,x_+.2],[np.mean(v)]*2,color='k',lw=0.8)
c.set_xticks([0,1]); c.set_xticklabels(['stock','fine-\ntuned']); c.set_xlim(-.5,1.5); c.set_ylim(-0.02,1)
c.set_ylabel('Nucleus detection (SA50)'); c.set_title(f'micro-sam, {len(tiles)} expert tiles\n(4-fold cross-validation)',fontsize=4.8,pad=2)

# ---------------- d: agreement with published masks
letter(138,47,'d'); ag=pd.read_csv('/work/mask_agreement_per_cell.csv'); ag['matched']=ag.cell_iou>=0.5
d=ax_mm(148,52,34,25)
for mi,m in enumerate(M):
    s=ag[ag.method==m]; rec=s.matched.mean(); nd=s.nuc_dice[s.matched]
    d.bar(mi-0.2,rec,0.36,color=COL[m],lw=0); d.bar(mi+0.2,nd.median(),0.36,color=COL[m],alpha=0.45,lw=0)
    d.errorbar(mi+0.2,nd.median(),yerr=[[nd.median()-nd.quantile(.25)],[nd.quantile(.75)-nd.median()]],color='k',lw=0.4,capsize=0)
d.set_xticks(range(3)); d.set_xticklabels([SHORT[m] for m in M],linespacing=1.0); d.set_ylim(0,1); d.set_xlim(-.6,2.6)
d.set_ylabel(f'Agreement with published\nmasks (n = {int((ag.method=="hil").sum())} cells)',fontsize=5)
from matplotlib.patches import Patch
d.legend(handles=[Patch(fc='#555555',label='cells found'),Patch(fc='#555555',alpha=0.45,label='nucleus Dice (median, IQR)')],frameon=False,fontsize=4.6,loc='upper center',bbox_to_anchor=(0.5,1.2),ncol=2,handlelength=0.8,columnspacing=0.6,handletextpad=0.3)

# ---------------- e: donor-level CD177+ vs CD177- (each run's own statistics)
def pct(pos,neg): return 100*(pos-neg)/neg
auto=pd.read_csv('/auto_stats/donor_contrasts.csv'); auto['d']=pct(auto.donor_pos_summary,auto.donor_neg_summary)
hil=pd.read_csv('/hil_stats/Donor_Level_Metrics.csv'); hil=hil[hil.session_subset=='all']
hp=hil.pivot_table(index=['metric','donor'],columns='cd177_status',values='donor_median_of_image_medians').reset_index(); hp['d']=pct(hp.pos,hp.neg)
cx=pd.read_csv('/codex_cells/cd177_cell_features.csv'); cx=cx[~cx.magnification.astype(str).str.contains('60')]
FEATS=[('Cell area',{'automated':'cell_area','hil':'cell_area_um2','codex':'cell_area'},'3.5×10$^{-52}$'),
       ('Nucleus area',{'automated':'nucleus_area','hil':'nucleus_area_um2','codex':'nucleus_area'},'0.2'),
       ('Nucleus/cell',{'automated':'nuclear_fraction_of_cell','hil':'nucleus_cell_area_ratio','codex':'nucleus_cell_area_ratio'},'2.4×10$^{-11}$')]
qa=pd.read_csv('/auto_stats/cd177_feature_statistics.csv').set_index('feature').p_adjusted
qh=pd.read_csv('/hil_stats/Paired_Donor_Results.csv').drop_duplicates('metric').set_index('metric').p_fdr_bh
qc=pd.read_csv('/codex_stats/cd177_feature_statistics.csv').set_index('feature').q_value_bh
letter(106,86,'e'); rng=np.random.default_rng(0)
for fi,(name,fm,pub) in enumerate(FEATS):
    e=ax_mm(115+fi*23,90,20,18)
    for mi,m in enumerate(M):
        if m=='automated': v=auto[auto.feature==fm[m]].d.values; q=qa[fm[m]]
        elif m=='hil': v=hp[hp.metric==fm[m]].d.values; q=qh[fm[m]]
        else:
            g2=cx.groupby(['donor','slide','condition'])[fm[m]].median().unstack('condition')
            pos=[k for k in g2.columns if 'pos' in str(k).lower()][0]; neg=[k for k in g2.columns if 'neg' in str(k).lower()][0]
            v=pct(g2[pos],g2[neg]).dropna().values; q=qc[fm[m]]
        v=np.clip(v,-60,150)
        e.scatter(mi+rng.uniform(-.17,.17,len(v)),v,s=3,color=COL[m],alpha=.85,lw=0)
        e.plot([mi-.28,mi+.28],[np.median(v)]*2,color='k',lw=0.7)
        e.text(mi,1.0,f'{q:.2f}',transform=e.get_xaxis_transform(),ha='center',va='bottom',fontsize=4.3)
    e.axhline(0,color='#999999',lw=0.4,ls='--'); e.set_xlim(-.6,2.6); e.set_ylim(-60,150); e.set_yticks([-50,0,50,100,150])
    e.set_xticks(range(3)); e.set_xticklabels([SHORT[m] for m in M],fontsize=4.3,linespacing=1.0)
    e.set_title(name,fontsize=5.5,pad=7)
    e.text(0.5,-0.3,f'published: p = {pub}',transform=e.transAxes,ha='center',va='top',fontsize=4.3,color=COL['published'])
    if fi==0:
        e.set_ylabel('CD177$^+$ vs CD177$^-$ (%)\none dot per donor',fontsize=5)
        e.text(-0.62,1.0,'q:',transform=e.transAxes,ha='left',va='bottom',fontsize=4.3)
    else: e.set_yticklabels([])
for ext in ['pdf','svg','png']: fig.savefig(f'/work/Fig_CD177_halfpage.{ext}',dpi=450 if ext=='png' else None)
