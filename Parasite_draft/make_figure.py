"""Parasite (Leishmania major / TUNEL) validation figure, tasks 1 and 2 (task 3 did not finish for the
automatic runs). Compares the published reference, Agentic-J + human (2 repeats), Agentic-J advanced,
Agentic-J fast and Codex CLI. Same style as Fig_CD177 v7: Nature type sizes (min 5 pt, labels 6 pt,
panel letters 8 pt bold). Footnotes live in the caption.

Sources
  parasite_figures/.../data            Agentic-J run outputs + article source data (JCI Insight 169020, 4C)
  parasite_figures/.../Codex runs      Codex CLI: the session that ran tasks 1 and 2 together (full_auto_task12);
                                       the single-task sessions (auto_1, auto_2) are not used
  inputs/task2_lab                     the lab's task-2 workbooks + ROI sets (EvaluationCalculation1-4.xlsx,
                                       Calculation1-4.zip) copied from the Parasite task folder
Auto-mode detections are matched to the lab's 15 task-2 parasites by optimal 1:1 centroid assignment
within MATCH_PX. Each lab parasite's index comes from its workbook row, tied to its coordinates through the
agent+human repeat-2 table (which records ROI geometry next to the value). Agentic-J run times are the
recorded session times in scripts/figure_workflows.py; Codex times are from each run_manifest.txt."""
import os, sys, csv, statistics as st
import numpy as np, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
import matplotlib.transforms as mtr
from matplotlib.patches import Rectangle, Patch, Ellipse
from matplotlib.lines import Line2D
import tifffile as tf
from skimage.measure import find_contours
from scipy.stats import mannwhitneyu

P = '/p'; sys.path.insert(0, f'{P}/scripts'); os.chdir(f'{P}/scripts')
import ij_roi, parasite_gt as gt
import figure_tasks23_dotplots as dp, score_auto_run as sa, build_workbook as bw
from figure_task1_overlays import load_run, stretch

FS,FL=5,6
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],'font.size':FS,
    'axes.linewidth':0.5,'xtick.major.width':0.5,'ytick.major.width':0.5,'xtick.major.size':2,'ytick.major.size':2,
    'xtick.labelsize':FS,'ytick.labelsize':FS,'axes.labelsize':FS,'axes.titlesize':FS,'axes.spines.top':False,'axes.spines.right':False,
    'svg.fonttype':'none','pdf.fonttype':42,'xtick.major.pad':1.2,'ytick.major.pad':1.2,'axes.labelpad':1.5})
MM=1/25.4; W,H=183,150
fig=plt.figure(figsize=(W*MM,H*MM))
def ax_mm(x,y,w,h,**k): return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H],**k)
def letter(x,y,s): fig.text(x/W,1-y/H,s,fontsize=8,fontweight='bold',va='top')
def head(x,y,s): fig.text(x/W,1-y/H,s,fontsize=FL,va='top')
COL={'published':'#7b3294','hil':'#0571b0','automated':'#1a9641','fast':'#a6611a','codex':'#636363'}
FILL={'manual':'#f7c873','hf':'#8ecdf2','auto':'#a6dba0','other':'#d9d9d9'}
SH={'published':'Published','hil':'Conventional Agentic-J + human','automated':'Conventional Agentic-J','fast':'Exploratory Agentic-J','codex':'Codex'}
TP,TN,UN='#e8575a','#7b68ee','#b0b0b0'          # TUNEL+ / TUNEL- / no host (dots)
TPi,TNi='#ff5c5c','#4fd2ff'                      # the same calls on the dark images
BAD='#c0392b'
MATCH_PX=15                                      # ~1.5 parasite diameters (4.8 um)
TOL=0.1                                          # index agreement tolerance
ROOT=f'{P}/data/parasite_task_results'; FAST=f'{ROOT}/fast-mode/session-20260925-182436'
CX=f'{P}/Codex runs'
CXC=f'{CX}/parasite_parasite_full_auto_task12_20260928_161707/workspace/output'      # both tasks, one session
PUB='published (Fig. 4C)'

# ================= task 1 data =================
nuclei,_=gt.load_task1(f'{ROOT}/oneframe_results/ExampleData_TunelMeasurement.xlsx')
rois,_,_=gt.load_task1_coords(f'{P}/data/one_frame',nuclei)
lab={n['roi']:n for n in nuclei}; rois=sorted(rois,key=lambda r:r['roi'])
hil1={}
with open(f'{ROOT}/hil_results/repeat_1/data/c2_tunel_final/C2_TUNEL_Classified_Nuclei.csv') as fh:   # both repeats: identical values
    for r in csv.DictReader(fh): hil1[int(r['roi_name'].split()[-1])]=(float(r['corrected_tunel_signal']),r['group'])
RUN1={'automated':f'{ROOT}/auto-mode/task_1','fast':f'{FAST}/task_1','codex':CXC}
R1={m:load_run(f'{d}/tunel_nuclei.csv') for m,d in RUN1.items()}
M1={m:tf.imread(f'{d}/tunel_nuclei_mask.tif') for m,d in RUN1.items()}
def nearest(runs,r):
    k,(x,y,v,g)=min(runs.items(),key=lambda kv:(kv[1][0]-r['x'])**2+(kv[1][1]-r['y'])**2); return v,g
T1={'published':[(lab[r['roi']]['corrected'],lab[r['roi']]['group']) for r in rois],'hil':[hil1[r['roi']] for r in rois]}
T1.update({m:[nearest(R1[m],r) for r in rois] for m in R1})
agree={m:sum(a[1]==b[1] for a,b in zip(T1[m],T1['published'])) for m in T1}
r1={m:np.corrcoef([v for v,_ in T1['published']],[v for v,_ in T1[m]])[0,1] for m in T1}
nobj1={m:len(R1[m]) for m in R1}
ov=ij_roi.read(f'{P}/data/one_frame/TUNEL_Overview.roi'); ox,oy,ow,oh=ov['left'],ov['top'],ov['w'],ov['h']
rawd=f'{ROOT}/auto-mode/task_1/leishmania_tunel_position40_task1/raw_images'
dapi=stretch(tf.imread(f'{rawd}/position40_DAPI_raw16bit_ch3.tif')[oy:oy+oh,ox:ox+ow],1.0,99.99)**0.6   # gamma keeps nuclei unsaturated
tun=stretch(tf.imread(f'{rawd}/position40_TUNEL_raw16bit_ch2.tif')[oy:oy+oh,ox:ox+ow],50.0,99.97)

# ================= task 2 data =================
LAB2='/work/inputs/task2_lab'
g2=gt.load_task2_coords(LAB2); w2,_=gt.load_task2(LAB2); lov=ij_roi.read(f'{LAB2}/LastFrame_Overview.roi')
def hil_t2(path,setc,grpc,valc):
    """Pair each agent+human row with the lab workbook row of the same group (and subset where known), nearest index."""
    with open(path) as fh: rows=list(csv.DictReader(fh))
    used,out=set(),[]
    for r in rows:
        s_,g,v=r[setc],r[grpc].strip(),float(r[valc])
        cand=[i for i,w in enumerate(w2) if w['group']==g and i not in used and any(x in s_ for x in w['subsets'])] or \
             [i for i,w in enumerate(w2) if w['group']==g and i not in used]
        i=min(cand,key=lambda i:abs(w2[i]['index']-v)); used.add(i); out.append((w2[i]['index'],v,g,r))
    return out
HP=f'{ROOT}/hil_results'
H2={'run 1':hil_t2(f'{HP}/repeat_1/data/proliferation_index/Proliferation_Index_Unique_Parasites.csv','source_sets','display_group','proliferation_index_value'),
    'run 2':hil_t2(f'{HP}/repeat_2/data/proliferation/position40_proliferation_index_unique_parasites.csv','set_id','group','proliferation_index')}
# coordinates -> lab index: repeat 2 records every ROI's full-field bbox beside the value it reproduced
for g in g2:
    best=min(H2['run 2'],key=lambda t:(float(t[3]['parasite_fullfield_x'])+float(t[3]['parasite_fullfield_w'])/2-g['x'])**2+
                                       (float(t[3]['parasite_fullfield_y'])+float(t[3]['parasite_fullfield_h'])/2-g['y'])**2)
    g['lab_index']=best[0]
RUN2={'automated':f'{ROOT}/auto-mode/task_2','fast':f'{FAST}/task_2','codex':CXC}
A2={}
for m,d in RUN2.items():
    objs=sa.objects_of(bw.auto_table(d,2)); kept,miss,extra=sa.match_by_distance(g2,objs,MATCH_PX)
    A2[m]=dict(objs=objs,kept=kept,extra=extra,
               pairs=[(g['lab_index'],r['value'],g['group'],r['group']) for g,r,_ in kept if r['value'] is not None])
t2,_=dp.collect(ROOT,f'{P}/data/published/jci.insight.169020.sdval.xlsx')
t2={k:t2[k] for k in (PUB,'Agent + human (run 1)','Advanced','Fast')}
for m,k in (('codex','Codex'),):
    t2[k]=dp.group_values([dict(index=o['value'],group=o['group']) for o in A2[m]['objs']])
t2={k:t2[k] for k in (PUB,'Codex','Fast','Advanced','Agent + human (run 1)')}       # figure order
K2=list(t2); MK=['published','codex','fast','automated','hil']
def mwu(g):
    a,b=g.get('TUNEL+'),g.get('TUNEL-')
    return mannwhitneyu(a,b,alternative='two-sided').pvalue if a and b and len(a)>1 and len(b)>1 else None
p2=[mwu(t2[k]) for k in K2]
# per-configuration task-2 scores
S2={}
for k,v in H2.items(): S2['hil '+k]=dict(found=15,extra=0,grp=15,err=[abs(a_-b_) for a_,b_,_,_ in v])
for m in A2: S2[m]=dict(found=len(A2[m]['kept']),extra=len(A2[m]['extra']),grp=sum(gg==rg for _,_,gg,rg in A2[m]['pairs']),
                        err=[abs(a_-b_) for a_,b_,_,_ in A2[m]['pairs']])

# ================= a  pipelines / verdicts =================
letter(0,0,'a'); a=ax_mm(0,0,183,46); a.axis('off'); a.set_xlim(0,183); a.set_ylim(46,0)
steps=['Registration\n(live → fixed)','Host / nucleus\nsegmentation','TUNEL call','Parasite\nlocalisation','Measurement\n& statistics']
rows=[('Published','published',[('by hand','manual'),('hand-drawn\nROIs','manual'),('by inspection','manual'),('15 ROIs + 3\nbackgrounds\neach, by hand','manual'),('ROI Manager →\nExcel;\nMann–Whitney','manual')]),
      ('Codex CLI','codex',[('SIFT + RANSAC,\nECC refined\n(1.6 px)','other'),('T1 DAPI\nthreshold; T2\nCellpose cyto3','other'),('largest-gap\ncutoff','other'),('red-channel\npeaks + water-\nshed: 10 kept','other'),('annulus\nbackground;\nno test','other')]),
      ('Exploratory\nAgentic-J','fast',[('TurboReg rigid\non brightfield','auto'),('T1 StarDist;\nT2 Cellpose\ncyto3 (7 hosts)','auto'),('largest-gap\ncutoff','auto'),('C2 + C3 spot\nevidence:\n22 kept','auto'),('annulus\nbackground;\nno test','auto')]),
      ('Conventional\nAgentic-J','automated',[('TurboReg /\nStackReg rigid','auto'),('T1 DoG + water-\nshed; T2 Cellpose\ncyto3 (6 hosts)','auto'),('2-component\nGMM','auto'),('micro-sam\nvit_b_lm:\n105 → 5 kept','auto'),('annulus\nbackground;\nno test','auto')]),
      ('Conventional\nAgentic-J\n+ human','hil',[('lab ROIs given,\norigin\ntranslated','hf'),('lab ROIs\n(given)','hf'),('largest-gap\ncutoff','hf'),('lab ROIs paired\n1:1 with\nbackgrounds','hf'),('native 16-bit;\nMann–Whitney','hf')])]
TIME={'published':'not recorded','hil':'32 min / n.r.','automated':'55 / 85 min','fast':'22 / 69 min','codex':'18 min (both)'}
def pfmt(p): return f'p = {p:.2g}'
w2n={k:sum(e<=TOL for e in S2['hil '+k]['err']) for k in H2}
sc=S2['codex']; pc=p2[K2.index('Codex')]
VER={'published':[('ref','2/6 nuclei TUNEL+'),('ref',f'8 TUNEL+ vs 7 TUNEL−, {pfmt(p2[K2.index(PUB)])}')],
     'hil':[('ok',f'{agree["hil"]}/6 calls agree, r = {r1["hil"]:.2f}'),('ok',f'15/15; {w2n["run 1"]}/15 within {TOL}; {pfmt(p2[K2.index("Agent + human (run 1)")])}')],
     'automated':[('ok',f'6/6 found, {agree["automated"]}/6 calls, r = {r1["automated"]:.2f}'),('x',f'{S2["automated"]["found"]}/15 found (+{S2["automated"]["extra"]}), no TUNEL+')],
     'fast':[('ok',f'6/6 (+1 split), {agree["fast"]}/6 calls, r = {r1["fast"]:.2f}'),('x',f'{S2["fast"]["found"]}/15 found (+{S2["fast"]["extra"]}), no host')],
     'codex':[('ok',f'6/6 found, {agree["codex"]}/6 calls, r = {r1["codex"]:.2f}'),
              ('ok' if sc['grp']==sc['found'] else 'x',f'{sc["found"]}/15 found, {sc["grp"]}/{sc["found"]} hosts right'+(f'; {pfmt(pc)}' if pc else ''))]}
x0,cw,gap,rh,y0=14.5,15.3,0.5,6.4,6.4
for i,s_ in enumerate(steps): a.text(x0+i*(cw+gap)+cw/2,3.2,s_,ha='center',va='center',fontsize=FS,fontweight='bold',linespacing=1.0)
xt=x0+5*(cw+gap)+0.8; vx=[xt+15.5,xt+15.5+36]
a.text(xt+7,3.2,'Time\nT1 / T2',ha='center',va='center',fontsize=FS,fontweight='bold',linespacing=1.0)
for j,t in enumerate(['Task 1  TUNEL classification\n(6 nuclei, fixed frame)','Task 2  proliferation, Fig. 4C\n(15 parasites, last live frame)']):
    a.text(vx[j]+17,3.2,t,ha='center',va='center',fontsize=FS,fontweight='bold',linespacing=1.0)
for r,(name,key,cells) in enumerate(rows):
    y=y0+r*(rh+0.5)
    a.text(0.3,y+rh/2,name,va='center',fontsize=FS,fontweight='bold',color=COL[key],linespacing=1.0)
    for i,(tt,kind) in enumerate(cells):
        xx=x0+i*(cw+gap)
        a.add_patch(Rectangle((xx,y),cw,rh,fc=FILL[kind],ec='k',lw=0.35))
        a.text(xx+cw/2,y+rh/2,tt,ha='center',va='center',fontsize=FS,linespacing=1.0)
    a.text(xt+7,y+rh/2,TIME[key],ha='center',va='center',fontsize=FS,color='#555555' if key=='published' else 'k')
    for j,(ok,txt) in enumerate(VER[key]):
        if ok!='ref':
            a.text(vx[j]+1,y+rh/2,'✓' if ok=='ok' else '✗',ha='center',va='center',fontsize=7,fontweight='bold',
                   color=COL[key],fontfamily='DejaVu Sans')
        a.text(vx[j]+3,y+rh/2,txt,va='center',fontsize=FS,color=COL['published'] if ok=='ref' else 'k')
xl=x0; rend=fig.canvas.get_renderer()
for k,t in [('manual','Manual'),('hf','Human feedback'),('auto','Automatic (Agentic-J)'),('other','Automatic (other agent)')]:
    a.add_patch(Rectangle((xl,42.4),2.4,2.2,fc=FILL[k],ec='k',lw=0.3)); tx=a.text(xl+3.1,43.5,t,va='center',fontsize=FS)
    xl+=3.1+tx.get_window_extent(rend).width/fig.dpi*25.4+2.2
LEG_END=xl-2.2; STEP_END=x0+5*(cw+gap)-gap

# ================= b  task-1 gallery =================
yb=48; letter(0,yb,'b'); head(7,yb+0.6,'Task 1: the same six nuclei, segmented and called by each pipeline')
def draw_calls(ax,kind):
    if kind in ('published','hil'):
        for r,(v,g) in zip(rois,T1[kind]):
            ax.add_patch(Ellipse((r['x']-ox,r['y']-oy),r['w'],r['h'],fill=False,lw=0.6,ec=TPi if g=='TUNEL+' else TNi))
    else:
        for k,(x,y,v,g) in R1[kind].items():
            for c in find_contours((M1[kind]==int(k)).astype(float),0.5):
                ax.plot(c[:,1]-ox,c[:,0]-oy,lw=0.6,color=TPi if g=='TUNEL+' else TNi)
cols=['raw','published','codex','fast','automated','hil']
clab=['raw TUNEL','Published','Codex','Exploratory\nAgentic-J','Conventional\nAgentic-J','Conventional\nAgentic-J\n+ human']
ccol=['k']+[COL[c] for c in cols[1:]]
c3=[r for r in rois if r['roi']==3][0]; c6=[r for r in rois if r['roi']==6][0]
views=[('Field',(0,ow,0,oh)),('Cell 3',(c3['x']-ox-45,c3['x']-ox+45,c3['y']-oy-45,c3['y']-oy+45)),
       ('Cell 6',(c6['x']-ox-45,c6['x']-ox+45,c6['y']-oy-45,c6['y']-oy+45))]
tw,tg,rp=14.2,0.45,16.0
for ri,(rl,(xa,xb,ya,yb_)) in enumerate(views):
    for ci,kind in enumerate(cols):
        ax=ax_mm(9+ci*(tw+tg),yb+8.5+ri*rp,tw,tw)
        if kind=='raw': ax.imshow(tun,cmap='magma',interpolation='nearest')
        else: ax.imshow(dapi,cmap='gray',interpolation='nearest'); draw_calls(ax,kind)
        ax.set_xlim(xa,xb); ax.set_ylim(yb_,ya); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(0.3)
        if ri==0:
            ax.set_title(clab[ci],fontsize=FS,color=ccol[ci],pad=1.5,linespacing=0.95)
            if kind in ('published','hil'): pts=[(r['x'],r['y'],r['roi']) for r in rois]
            elif kind!='raw':               # the run's own objects, labelled with the nearest lab nucleus
                near={}
                for x,y,_v,_g in R1[kind].values(): near.setdefault(min(rois,key=lambda r:(r['x']-x)**2+(r['y']-y)**2)['roi'],[]).append((x,y))
                pts=[]
                for n,v in near.items():             # a nucleus split into pieces: every piece labelled, 3a / 3b, left to right
                    v=sorted(v)
                    pts+=[(x,y,f'{n}{"ab"[j]}' if len(v)>1 else str(n),'left' if (len(v)>1 and j==0) else 'right') for j,(x,y) in enumerate(v)]
            else: pts=[]
            if kind in ('published','hil'): pts=[(x,y,str(n),'right') for x,y,n in pts]
            for x,y,n,side in pts:
                if side=='left': ax.text(x-ox+2,y-oy+14,n,color='w',fontsize=4.2,ha='center',va='top')
                else: ax.text(x-ox+20,y-oy+(30 if n=='1' else -6),n,color='w',fontsize=4.2)
        if ci==0: ax.set_ylabel(rl,fontsize=FS,rotation=0,ha='right',va='center',labelpad=1.5)
        if ci==0 and ri<2:   # scale bars: 20 um on the field, 5 um on the zooms (0.3223 um/px)
            L,lab_=(20,'20 µm') if ri==0 else (5,'5 µm'); px=L/gt.PIXEL_UM; span=xb-xa
            ax.plot([xa+0.06*span,xa+0.06*span+px],[yb_-0.08*span]*2,color='w',lw=0.9)
            ax.text(xa+0.06*span,yb_-0.12*span,lab_,color='w',fontsize=4.2,va='bottom')
lx=ax_mm(9,yb+56.3,88,3); lx.axis('off'); lx.set_xlim(0,88); lx.set_ylim(0,1)
for x,c,t in [(2,TPi,'labeled TUNEL+'),(22,TNi,'labeled TUNEL−')]:
    lx.plot([x,x+4],[0.5,0.5],color=c,lw=1); lx.text(x+5,0.5,t,va='center',fontsize=FS)

# corrected TUNEL signal per nucleus
s=ax_mm(12,yb+64.5,85,31.5)
MARK={'published':'o','hil':'D','automated':'s','fast':'^','codex':'v'}
xs=np.arange(6); ML=['published','codex','fast','automated','hil']
for m,off in zip(ML,np.linspace(-0.2,0.2,5)):
    v=[a_[0] for a_ in T1[m]]; g=[a_[1] for a_ in T1[m]]
    s.scatter(xs+off,v,marker=MARK[m],s=9,color=COL[m],lw=0,zorder=3,label=SH[m])
    hit=[i for i in range(6) if g[i]=='TUNEL+']
    s.scatter(xs[hit]+off,np.array(v)[hit],marker=MARK[m],s=34,facecolors='none',edgecolors=TP,lw=0.7,zorder=4)
s.set_xticks(xs); s.set_xticklabels([f'Cell {r["roi"]}' for r in rois]); s.set_xlim(-0.5,5.5)
s.set_ylabel('Corrected TUNEL signal (a.u.)'); s.set_ylim(0,31); s.set_yticks([0,10,20,30])
s.grid(axis='y',lw=0.25,color='#dddddd'); s.set_axisbelow(True)
h_,l_=s.get_legend_handles_labels()
h_.append(Line2D([],[],marker='o',ms=4,lw=0,mfc='none',mec=TP,mew=0.7)); l_.append('labeled TUNEL+')
s.legend(h_,l_,frameon=False,fontsize=FS,loc='upper left',ncol=1,handletextpad=0.3,borderaxespad=0.2,labelspacing=0.3)
s.set_title('Corrected TUNEL signal per nucleus',fontsize=FS,pad=2)

# ================= c  agreement with the lab: scorecard =================
xc=101; letter(xc,yb,'c'); head(xc+7,yb+0.6,'Agreement with the published analysis')
SCM=[('codex','codex','Codex'),('fast','fast','Exploratory\nAgentic-J'),('automated','automated','Conventional\nAgentic-J'),('hil run 1','hil','Conventional\nAgentic-J\n+ human')]
def t1val(m): return agree['hil' if m.startswith('hil') else m]
def within(m): return sum(e<=TOL for e in S2[m]['err'])
SROWS=[('Task 1','TUNEL calls agree',lambda m:(t1val(m),6)),
       ('Task 2','Lab parasites found',lambda m:(S2[m]['found'],15)),
       ('','Host group correct',lambda m:(S2[m]['grp'],S2[m]['found'])),
       ('',f'Index within {TOL} of lab',lambda m:(within(m),len(S2[m]['err']))),
       ('','False detections',lambda m:S2[m]['extra'])]
sc=ax_mm(xc,yb+4.5,82,28); sc.axis('off'); sc.set_xlim(0,82); sc.set_ylim(28,0)
LX,CX0,CW,CG,RY0,RH,RG=2,32,11.9,0.6,6.6,3.5,0.55
for j,(m,cm,lab_) in enumerate(SCM):
    sc.text(CX0+j*(CW+CG)+CW/2,3.0,lab_,ha='center',va='center',fontsize=FS,color=COL[cm],linespacing=0.95,fontweight='bold')
y=RY0
for i,(grp,lab_,fn) in enumerate(SROWS):
    if grp=='Task 2': y+=1.2                                     # small break between the tasks
    if grp: sc.text(LX,y+RH/2,grp,va='center',fontsize=FS,fontweight='bold')
    sc.text(LX+9,y+RH/2,lab_,va='center',fontsize=FS,color='#555555' if lab_=='False detections' else 'k')
    for j,(m,cm,_) in enumerate(SCM):
        x=CX0+j*(CW+CG); v=fn(m)
        if lab_=='False detections':                              # lower is better: plain count, no colour scale
            sc.add_patch(Rectangle((x,y),CW,RH,fc='#f2f2f2',ec='none'))
            sc.text(x+CW/2,y+RH/2,str(v),ha='center',va='center',fontsize=FS,color='#333333')
            continue
        k,n=v
        if n==0:
            sc.add_patch(Rectangle((x,y),CW,RH,fc='w',ec='#bbbbbb',lw=0.3,ls='--'))
            sc.text(x+CW/2,y+RH/2,'–',ha='center',va='center',fontsize=FS,color='#888888'); continue
        fr=k/n; sc.add_patch(Rectangle((x,y),CW,RH,fc=COL[cm],alpha=0.12+0.88*fr,ec='none'))
        sc.text(x+CW/2,y+RH/2,f'{k} / {n}',ha='center',va='center',fontsize=FS,color='w' if fr>0.55 else 'k',fontweight='bold' if fr>0.55 else 'normal')
    y+=RH+RG

# ================= d  task 2 = Fig. 4C =================
ye=yb+35; letter(xc,ye,'d'); head(xc+7,ye+0.6,'Task 2: proliferation index distributions')
SRC=['Published','Codex','Exploratory','Conventional','Conventional\n+ human']
def stars(p): return '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
rng=np.random.default_rng(0)
EH=40; d=ax_mm(xc+9,ye+6,57,EH); lim=(0,10.6)
for i,k in enumerate(K2):
    g=t2[k]
    if i%2==0: d.axvspan(i-0.5,i+0.5,color='#f1f1f1',lw=0,zorder=0)      # alternating bands separate the runs
    for grp,dx,c in [('TUNEL+',-0.19,TP),('TUNEL-',0.19,TN),('Unassigned',0,UN)]:
        v=g.get(grp)
        if not v: continue
        n_out=sum(x_<lim[0] for x_ in v)
        d.scatter(i+dx+rng.uniform(-0.05,0.05,len(v)),np.clip(v,*lim),s=6,color=c,lw=0,alpha=0.85,zorder=2)
        d.plot([i+dx-0.14,i+dx+0.14],[st.median(v)]*2,color='k',lw=0.8,zorder=4,solid_capstyle='butt')
        if n_out: d.text(i+dx,0.15,f'↓{n_out}',ha='center',va='bottom',fontsize=4.2,color=c)
    if p2[i] is not None: d.text(i,lim[1],stars(p2[i]),ha='center',va='top',fontsize=FS)
d.set_xticks(range(len(K2))); d.set_xticklabels(SRC,fontsize=FS,rotation=45,ha='right',rotation_mode='anchor',linespacing=0.95)
for tl,m in zip(d.get_xticklabels(),MK): tl.set_color(COL[m])
d.tick_params(axis='x',length=0); d.set_ylim(*lim); d.set_yticks([0,5,10]); d.set_xlim(-0.5,len(K2)-0.5)
d.grid(axis='y',lw=0.25,color='#dddddd'); d.set_axisbelow(True); d.set_ylabel('Proliferation index\n(10 − R/G)',linespacing=0.95)
d.legend(handles=[Line2D([],[],marker='o',lw=0,ms=2.5,color=TP,label='TUNEL+ host'),Line2D([],[],marker='o',lw=0,ms=2.5,color=TN,label='TUNEL− host'),
                  Line2D([],[],marker='o',lw=0,ms=2.5,color=UN,label='no host'),Line2D([],[],color='k',lw=0.8,label='median')],
         frameon=False,fontsize=FS,loc='upper left',bbox_to_anchor=(1.02,1.0),handlelength=0.9,handletextpad=0.3,borderaxespad=0,labelspacing=0.35)
def dmed(g,B=4000):
    a_,b_=np.asarray(g.get('TUNEL+',[])),np.asarray(g.get('TUNEL-',[]))
    if len(a_)<2 or len(b_)<2: return None
    r=np.random.default_rng(1); bs=[np.median(r.choice(a_,len(a_)))-np.median(r.choice(b_,len(b_))) for _ in range(B)]
    return np.median(a_)-np.median(b_),*np.percentile(bs,[2.5,97.5])
DM=[dmed(t2[k]) for k in K2]
ty=mtr.blended_transform_factory(d.transData,fig.transFigure); yrow=1-(ye+6+EH+13.2)/H; yci=1-(ye+6+EH+15.0)/H
d.text(-0.5,yrow,'Δ median  ',ha='right',va='center',fontsize=FS,fontweight='bold',transform=ty)
d.text(-0.5,yci,'95% CI  ',ha='right',va='center',fontsize=4.3,color='#555555',transform=ty)
for i,r in enumerate(DM):
    npos=len(t2[K2[i]].get('TUNEL+',[]))
    d.text(i,yrow,(f'{r[0]:.2f}' if r is not None else ('no\nTUNEL+' if npos==0 else f'only {npos}\nTUNEL+')),ha='center',va='center',
           fontsize=FS if r is not None else 4.2,color=COL[MK[i]],style='normal' if r is not None else 'italic',linespacing=0.9,transform=ty)
    if r is not None: d.text(i,yci,f'[{r[1]:.1f}, {r[2]:.1f}]',ha='center',va='center',fontsize=4.0,color='#555555',transform=ty)

for ext in ['pdf','svg','png']: fig.savefig(f'/work/Fig_Parasite.{ext}',dpi=450 if ext=='png' else None)
# ---- print what was plotted, so a rebuild can be checked against the figure
print('legend end',round(LEG_END,1),'steps end',round(STEP_END,1)); print('task1 objects',nobj1,' calls agree',agree,' r',{m:round(v,3) for m,v in r1.items()})
for m in T1: print('  ',m,[f'{v:.2f}{"+" if g=="TUNEL+" else ""}' for v,g in T1[m]])
for m,v in S2.items(): print('task2',f'{m:10}',{k:(v[k] if k!='err' else [round(e,3) for e in v[k]]) for k in v})
for m in A2: print('   ',m,'matched',[(g['name'],g['subset'],g['group'],r['group'],round(dd,1)) for g,r,dd in A2[m]['kept']])
for k,p,r in zip(K2,p2,DM): print(f'task2 {k:24}',{g:(len(v),round(st.median(v),3)) for g,v in t2[k].items()},'p=',None if p is None else f'{p:.3g}','dmed/CI=',None if r is None else np.round(r,3))
