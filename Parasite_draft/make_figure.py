"""Parasite (Leishmania major / TUNEL) validation figure, tasks 1 and 2 only (task 3 did not finish for
the automatic runs). Same style as Fig_CD177 v7: Nature type sizes (min 5 pt, labels 6 pt, letters 8 pt bold).

Sources: run outputs + article source data (JCI Insight 169020, sheet 4C) in parasite_figures/, loaded with
its own scripts; the lab's task-2 workbooks and ROI sets (EvaluationCalculation1-4.xlsx, Calculation1-4.zip,
copied from the Parasite task folder) in inputs/task2_lab/. Auto-mode detections are matched to the lab's
15 parasites by optimal 1:1 centroid assignment within MATCH_PX. Run times are the recorded session times
listed in parasite_figures/.../scripts/figure_workflows.py."""
import os, sys, csv, statistics as st
import numpy as np, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch, Ellipse
from matplotlib.lines import Line2D
import tifffile as tf
from skimage.measure import find_contours
from scipy.stats import mannwhitneyu

P = '/p'; sys.path.insert(0, f'{P}/scripts'); os.chdir(f'{P}/scripts')
import ij_roi, parasite_gt as gt, xlsx_min
import figure_tasks23_dotplots as dp, score_auto_run as sa, build_workbook as bw
from figure_task1_overlays import load_run, stretch

FS,FL=5,6
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],'font.size':FS,
    'axes.linewidth':0.5,'xtick.major.width':0.5,'ytick.major.width':0.5,'xtick.major.size':2,'ytick.major.size':2,
    'xtick.labelsize':FS,'ytick.labelsize':FS,'axes.labelsize':FS,'axes.titlesize':FS,'axes.spines.top':False,'axes.spines.right':False,
    'svg.fonttype':'none','pdf.fonttype':42,'xtick.major.pad':1.2,'ytick.major.pad':1.2,'axes.labelpad':1.5})
MM=1/25.4; W,H=183,145
fig=plt.figure(figsize=(W*MM,H*MM))
def ax_mm(x,y,w,h,**k): return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H],**k)
def letter(x,y,s): fig.text(x/W,1-y/H,s,fontsize=8,fontweight='bold',va='top')
def head(x,y,s): fig.text(x/W,1-y/H,s,fontsize=FL,va='top')
COL={'published':'#7b3294','hil':'#0571b0','automated':'#1a9641','fast':'#a6611a'}
FILL={'manual':'#f7c873','hf':'#8ecdf2','auto':'#a6dba0','none':'#ffffff'}
TP,TN,UN='#e8575a','#7b68ee','#b0b0b0'          # TUNEL+ / TUNEL- / unassigned (dots)
TPi,TNi='#ff5c5c','#4fd2ff'                      # same calls on the dark images
OK,BAD='#333333','#c0392b'
ROOT=f'{P}/data/parasite_task_results'; FAST=f'{ROOT}/fast-mode/session-20260925-182436'

# ================= data =================

MATCH_PX=15                                   # ~1.5 parasite diameters (4.8 um)
# ================= data =================
t2,_=dp.collect(ROOT,f'{P}/data/published/jci.insight.169020.sdval.xlsx')
K2=['published (Fig. 4C)','Agent + human (run 1)','Agent + human (run 2)','Advanced','Fast']
MK=['published','hil','hil','automated','fast']
def mwu(g):
    a,b=g.get('TUNEL+'),g.get('TUNEL-')
    return mannwhitneyu(a,b,alternative='two-sided').pvalue if a and b and len(a)>1 and len(b)>1 else None
p2=[mwu(t2[k]) for k in K2]

# task 1: lab ROIs + calls, HIL (both repeats gave identical values; repeat 1 file), Advanced, Fast
nuclei,_=gt.load_task1(f'{ROOT}/oneframe_results/ExampleData_TunelMeasurement.xlsx')
rois,_,_=gt.load_task1_coords(f'{P}/data/one_frame',nuclei)
lab={n['roi']:n for n in nuclei}; rois=sorted(rois,key=lambda r:r['roi'])
hil1={}
with open(f'{ROOT}/hil_results/repeat_1/data/c2_tunel_final/C2_TUNEL_Classified_Nuclei.csv') as fh:
    for r in csv.DictReader(fh): hil1[int(r['roi_name'].split()[-1])]=(float(r['corrected_tunel_signal']),r['group'])
auto=load_run(f'{ROOT}/auto-mode/task_1/tunel_nuclei.csv'); fast=load_run(f'{FAST}/task_1/tunel_nuclei.csv')
amask=tf.imread(f'{ROOT}/auto-mode/task_1/tunel_nuclei_mask.tif'); fmask=tf.imread(f'{FAST}/task_1/tunel_nuclei_mask.tif')
def nearest(runs,r):
    k,(x,y,v,g)=min(runs.items(),key=lambda kv:(kv[1][0]-r['x'])**2+(kv[1][1]-r['y'])**2); return v,g
T1={'published':[(lab[r['roi']]['corrected'],lab[r['roi']]['group']) for r in rois],
    'hil':[hil1[r['roi']] for r in rois],'automated':[nearest(auto,r) for r in rois],'fast':[nearest(fast,r) for r in rois]}
agree={m:sum(a[1]==b[1] for a,b in zip(T1[m],T1['published'])) for m in T1}
ov=ij_roi.read(f'{P}/data/one_frame/TUNEL_Overview.roi'); ox,oy,ow,oh=ov['left'],ov['top'],ov['w'],ov['h']
rawd=f'{ROOT}/auto-mode/task_1/leishmania_tunel_position40_task1/raw_images'
dapi=stretch(tf.imread(f'{rawd}/position40_DAPI_raw16bit_ch3.tif')[oy:oy+oh,ox:ox+ow],1.0,99.99)**0.6   # gamma keeps nuclei unsaturated in the zooms
tun=stretch(tf.imread(f'{rawd}/position40_TUNEL_raw16bit_ch2.tif')[oy:oy+oh,ox:ox+ow],50.0,99.97)

# task 2: the lab's 15 parasites (coordinates from the ROI sets, index from the workbooks)
LAB2='/work/inputs/task2_lab'
g2=gt.load_task2_coords(LAB2); w2,_=gt.load_task2(LAB2)
lov=ij_roi.read(f'{LAB2}/LastFrame_Overview.roi')
def hil_t2(path,setc,grpc,valc):
    with open(path) as fh: rows=[(r[setc],r[grpc].strip(),float(r[valc])) for r in csv.DictReader(fh)]
    used,pairs=set(),[]
    for s_,g,v in rows:
        cand=[i for i,w in enumerate(w2) if w['group']==g and i not in used and any(x in s_ for x in w['subsets'])] or \
             [i for i,w in enumerate(w2) if w['group']==g and i not in used]
        i=min(cand,key=lambda i:abs(w2[i]['index']-v)); used.add(i); pairs.append((w2[i]['index'],v,g))
    return pairs
HP=f'{ROOT}/hil_results'
H2={'run 1':hil_t2(f'{HP}/repeat_1/data/proliferation_index/Proliferation_Index_Unique_Parasites.csv','source_sets','display_group','proliferation_index_value'),
    'run 2':hil_t2(f'{HP}/repeat_2/data/proliferation/position40_proliferation_index_unique_parasites.csv','set_id','group','proliferation_index')}
A2={}
for m,rel in (('automated','auto-mode/task_2'),('fast','fast-mode/session-20260925-182436/task_2')):
    objs=sa.objects_of(bw.auto_table(f'{ROOT}/{rel}',2)); kept,miss,extra=sa.match_by_distance(g2,objs,MATCH_PX)
    A2[m]=dict(objs=objs,kept=kept,extra=extra)
ctrl={w['label']:w['index'] for w in w2 if w['group']=='TUNEL-'}
def err_auto(m):   # index error only where the lab label is unambiguous (controls C<n> = Control<n>)
    return [(ctrl[g['gt_label']],r['value']) for g,r,d in A2[m]['kept'] if g['gt_label'] in ctrl and r['value'] is not None]


# ================= a  pipelines / verdicts =================
letter(0,0,'a'); a=ax_mm(0,0,183,40); a.axis('off'); a.set_xlim(0,183); a.set_ylim(40,0)
steps=['Registration\n(live → fixed)','Host / nucleus\nsegmentation','TUNEL call','Parasite\nlocalisation','Measurement\n& statistics']
rows=[('Published','published',[('by hand','manual'),('hand-drawn\nROIs','manual'),('by inspection','manual'),('15 ROIs + 3\nbackgrounds\neach, by hand','manual'),('ROI Manager →\nExcel;\nMann–Whitney','manual')]),
      ('Agentic-J\n+ human','hil',[('lab ROIs given,\norigin\ntranslated','hf'),('lab ROIs\n(given)','hf'),('largest-gap\ncutoff','hf'),('lab ROIs paired\n1:1 with\nbackgrounds','hf'),('native 16-bit;\nMann–Whitney','hf')]),
      ('Agentic-J\n(advanced)','automated',[('TurboReg /\nStackReg rigid','auto'),('T1 DoG + water-\nshed; T2 Cellpose\ncyto3 (6 hosts)','auto'),('2-component\nGMM','auto'),('micro-sam\nvit_b_lm:\n105 → 5 kept','auto'),('annulus\nbackground;\nno test','auto')]),
      ('Agentic-J\n(fast)','fast',[('TurboReg rigid\non brightfield','auto'),('T1 StarDist;\nT2 Cellpose\ncyto3 (7 hosts)','auto'),('largest-gap\ncutoff','auto'),('C2 + C3 spot\nevidence:\n22 kept','auto'),('annulus\nbackground;\nno test','auto')])]
TIME={'published':'not recorded','hil':'32 min / n.r.','automated':'55 / 85 min','fast':'22 / 69 min'}
def pfmt(p): return f'p = {p:.2g}'
r1={m:np.corrcoef([v for v,_ in T1['published']],[v for v,_ in T1[m]])[0,1] for m in T1}
TOL=0.1
e2={k:max(abs(a_-b_) for a_,b_,_ in v) for k,v in H2.items()}
w2n={k:sum(abs(a_-b_)<=TOL for a_,b_,_ in v) for k,v in H2.items()}
n2={m:len(A2[m]['kept']) for m in A2}; x2={m:len(A2[m]['extra']) for m in A2}
VER={'published':[('ref','2/6 nuclei TUNEL+'),('ref',f'8 TUNEL+ vs 7 TUNEL−, {pfmt(p2[0])}')],
     'hil':[('ok',f'{agree["hil"]}/6 calls agree, r = {r1["hil"]:.2f}'),('ok',f'15/15; {w2n["run 1"]} | {w2n["run 2"]} within {TOL}; {pfmt(p2[1])}')],
     'automated':[('ok',f'6/6 found, {agree["automated"]}/6 calls, r = {r1["automated"]:.2f}'),('x',f'{n2["automated"]}/15 found (+{x2["automated"]}), no TUNEL+')],
     'fast':[('ok',f'6/6 found (+1 split), {agree["fast"]}/6 calls'),('x',f'{n2["fast"]}/15 found (+{x2["fast"]}), no host')]}
x0,cw,gap,rh,y0=13,15.6,0.5,6.4,6.4
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
            a.text(vx[j]+1,y+rh/2,'\u2713' if ok=='ok' else '\u2717',ha='center',va='center',fontsize=7,fontweight='bold',
                   color=COL[key] if ok=='ok' else BAD,fontfamily='DejaVu Sans')
        a.text(vx[j]+3,y+rh/2,txt,va='center',fontsize=FS,color=COL['published'] if ok=='ref' else 'k')
xl=x0
for k,t in [('manual','Manual (analyst)'),('hf','Human feedback (user)'),('auto','Automatic (Agentic-J)')]:
    a.add_patch(Rectangle((xl,35.9),2.4,2.2,fc=FILL[k],ec='k',lw=0.3)); a.text(xl+3.2,37.0,t,va='center',fontsize=FS); xl+=3.2+len(t)*0.98+3.5

# ================= b  task-1 gallery =================
yb=41.5; letter(0,yb,'b'); head(7,yb+0.6,'Task 1: the same six nuclei, segmented and called by each pipeline')
def draw_calls(ax,kind):
    if kind in ('published','hil'):
        for r,(v,g) in zip(rois,T1[kind]):
            ax.add_patch(Ellipse((r['x']-ox,r['y']-oy),r['w'],r['h'],fill=False,lw=0.6,ec=TPi if g=='TUNEL+' else TNi))
    else:
        mask,runs=(amask,auto) if kind=='automated' else (fmask,fast)
        for k,(x,y,v,g) in runs.items():
            for c in find_contours((mask==int(k)).astype(float),0.5):
                ax.plot(c[:,1]-ox,c[:,0]-oy,lw=0.6,color=TPi if g=='TUNEL+' else TNi)
cols=['raw','published','hil','automated','fast']
clab=['raw TUNEL','Published','Agentic-J\n+ human','Agentic-J\n(advanced)','Agentic-J\n(fast)']
ccol=['k',COL['published'],COL['hil'],COL['automated'],COL['fast']]
c3=[r for r in rois if r['roi']==3][0]; c6=[r for r in rois if r['roi']==6][0]
views=[('Field',(0,ow,0,oh)),('Cell 3',(c3['x']-ox-45,c3['x']-ox+45,c3['y']-oy-45,c3['y']-oy+45)),
       ('Cell 6',(c6['x']-ox-45,c6['x']-ox+45,c6['y']-oy-45,c6['y']-oy+45))]
notes={(1,4):'split in two',(2,2):'called TUNEL−',(2,3):'called TUNEL−',(2,4):'called TUNEL−',(2,1):'TUNEL+ by eye'}
tw,tg,rp=15.6,0.45,17.6
for ri,(rl,(xa,xb,ya,yb_)) in enumerate(views):
    for ci,kind in enumerate(cols):
        ax=ax_mm(9+ci*(tw+tg),yb+7+ri*rp,tw,tw)
        if kind=='raw': ax.imshow(tun,cmap='magma',interpolation='nearest')
        else: ax.imshow(dapi,cmap='gray',interpolation='nearest'); draw_calls(ax,kind)
        ax.set_xlim(xa,xb); ax.set_ylim(yb_,ya); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(0.3)
        if ri==0:
            ax.set_title(clab[ci],fontsize=FS,color=ccol[ci],pad=1.5,linespacing=0.95)
            if kind in ('published','hil'):
                for r in rois: ax.text(r['x']-ox+20,r['y']-oy+(30 if r['roi']==1 else -6),str(r['roi']),color='w',fontsize=4.2)
        if ci==0: ax.set_ylabel(rl,fontsize=FS,rotation=0,ha='right',va='center',labelpad=1.5)
        if (ri,ci) in notes: ax.text(0.5,-0.03,notes[(ri,ci)],transform=ax.transAxes,ha='center',va='top',fontsize=FS,
                                     color=COL['published'] if ci==1 else BAD)
        if ci==0 and ri<2:   # scale bars: 20 um on the field, 5 um on the zooms (0.3223 um/px)
            L,lab_=(20,'20 µm') if ri==0 else (5,'5 µm'); px=L/gt.PIXEL_UM; span=xb-xa
            ax.plot([xa+0.06*span,xa+0.06*span+px],[yb_-0.08*span]*2,color='w',lw=0.9)
            ax.text(xa+0.06*span,yb_-0.12*span,lab_,color='w',fontsize=4.2,va='bottom')
lx=ax_mm(9,yb+61.2,88,3); lx.axis('off'); lx.set_xlim(0,88); lx.set_ylim(0,1)
for x,c,t in [(2,TPi,'called TUNEL+'),(22,TNi,'called TUNEL−')]:
    lx.plot([x,x+4],[0.5,0.5],color=c,lw=1); lx.text(x+5,0.5,t,va='center',fontsize=FS)
lx.text(42,0.5,'(outlines on DAPI; IDs are the lab\'s)',va='center',fontsize=FS,color='#555555')

# corrected TUNEL signal per nucleus
s=ax_mm(12,yb+70,83,27.5)
MARK={'published':'o','hil':'D','automated':'s','fast':'^'}; SHN={'published':'Published','hil':'Agentic-J + human','automated':'Agentic-J (advanced)','fast':'Agentic-J (fast)'}
xs=np.arange(6)
for m,off in zip(['published','hil','automated','fast'],[-0.15,-0.05,0.05,0.15]):
    v=[a_[0] for a_ in T1[m]]; g=[a_[1] for a_ in T1[m]]
    s.plot(xs+off,v,'-',color=COL[m],lw=0.6,alpha=0.6,zorder=2)
    s.scatter(xs+off,v,marker=MARK[m],s=9,color=COL[m],lw=0,zorder=3,label=SHN[m])
    hit=[i for i in range(6) if g[i]=='TUNEL+']
    s.scatter(xs[hit]+off,np.array(v)[hit],marker=MARK[m],s=34,facecolors='none',edgecolors=TP,lw=0.7,zorder=4)
s.set_xticks(xs); s.set_xticklabels([f'Cell {r["roi"]}' for r in rois]); s.set_xlim(-0.5,5.5)
s.set_ylabel('Corrected TUNEL signal (a.u.)'); s.set_ylim(0,31); s.set_yticks([0,10,20,30])
s.grid(axis='y',lw=0.25,color='#dddddd'); s.set_axisbelow(True)
h_,l_=s.get_legend_handles_labels()
h_.append(Line2D([],[],marker='o',ms=4,lw=0,mfc='none',mec=TP,mew=0.7)); l_.append('called TUNEL+')
s.legend(h_,l_,frameon=False,fontsize=FS,loc='upper left',ncol=1,handletextpad=0.3,borderaxespad=0.2,labelspacing=0.3)

# ================= c  agreement with the lab =================
xc=101; letter(xc,yb,'c'); head(xc+7,yb+0.6,'Agreement with the published analysis')
M4=['hil','hil','automated','fast']
e2a={m:[abs(a_-b_) for a_,b_ in err_auto(m)] for m in A2}
panels=[('Task 1: TUNEL\ncalls (of 6)',[agree['hil']/6]*2+[agree['automated']/6,agree['fast']/6],[f'{agree["hil"]}']*2+[f'{agree["automated"]}',f'{agree["fast"]}'],(0,1.25),[0,0.5,1],['','','','']),
        ('Task 1: signal vs\nlab (Pearson r)',[r1['hil']]*2+[r1['automated'],r1['fast']],[f'{r1["hil"]:.2f}']*2+[f'{r1["automated"]:.2f}',f'{r1["fast"]:.2f}'],(0,1.25),[0,0.5,1],['','','','']),
        ('Task 2: lab parasites\nfound (of 15)',[1,1,n2['automated']/15,n2['fast']/15],['15','15',str(n2['automated']),str(n2['fast'])],(0,1.25),[0,0.5,1],['','',f'+{x2["automated"]}',f'+{x2["fast"]}']),
        (f'Task 2: index within\n{TOL} of lab',[w2n['run 1']/15,w2n['run 2']/15,(sum(x<=TOL for x in e2a['automated'])/len(e2a['automated'])) if e2a['automated'] else 0,0],
         [f'{w2n["run 1"]}',f'{w2n["run 2"]}',f'{sum(x<=TOL for x in e2a["automated"])}/{len(e2a["automated"])}','–'],(0,1.25),[0,0.5,1],['','','',''])]
for j,(t,vals,labs,yl,yt,ann) in enumerate(panels):
    d=ax_mm(xc+7.5+j*19.2,yb+9.5,15.5,13)
    for i,(m,v) in enumerate(zip(M4,vals)):
        d.bar(i,v,0.66,color=COL[m],lw=0,alpha=1 if i!=1 else 0.6)
        d.text(i,v+0.03*(yl[1]),labs[i],ha='center',va='bottom',fontsize=4.3)
        if ann[i]: d.text(i,v+0.12*yl[1],ann[i],ha='center',va='bottom',fontsize=4.3,color=BAD if j==2 else '#555555')
    d.set_ylim(*yl); d.set_yticks(yt); d.set_xlim(-0.6,3.6); d.set_xticks([]); d.spines['left'].set_bounds(yt[0],yt[-1])
    d.set_title(t,fontsize=FS,pad=2,linespacing=0.95); d.grid(axis='y',lw=0.25,color='#dddddd'); d.set_axisbelow(True)
    if j: d.set_yticklabels([])
fig.legend(handles=[Patch(fc=COL['hil'],label='Agentic-J + human (run 1 | run 2)'),Patch(fc=COL['automated'],label='Agentic-J (advanced)'),Patch(fc=COL['fast'],label='Agentic-J (fast)')],
           frameon=False,fontsize=FS,ncol=3,loc='upper left',bbox_to_anchor=((xc+7)/W,1-(yb+24.3)/H),handlelength=1.0,handletextpad=0.4,columnspacing=1.0,borderaxespad=0)

# ================= d  task 2: where parasites were found, and the value per parasite =================
yd=yb+29; letter(xc,yd,'d'); head(xc+7,yd+0.6,'Task 2: parasite detection and per-parasite index')
UM=gt.PIXEL_UM; lx0,ly0=lov['left'],lov['top']
mp=ax_mm(xc+8,yd+6,30,30)
for g in g2:
    mp.add_patch(plt.Circle(((g['x']-lx0)*UM,(g['y']-ly0)*UM),MATCH_PX*UM,fc='none',ec='#cccccc',lw=0.3,ls=':'))
    mp.scatter((g['x']-lx0)*UM,(g['y']-ly0)*UM,s=9,color=TP if g['group']=='TUNEL+' else TN,lw=0,zorder=3)
nout={}
for m,mk in (('automated','s'),('fast','^')):
    hit={id(r) for _g,r,_d in A2[m]['kept']}; o=A2[m]['objs']
    xs_=np.array([(r['x']-lx0)*UM for r in o]); ys_=np.array([(r['y']-ly0)*UM for r in o])
    inside=(xs_>=0)&(xs_<=lov['w']*UM)&(ys_>=0)&(ys_<=lov['h']*UM); nout[m]=int((~inside).sum())
    for r,x_,y_ in zip(o,xs_,ys_):
        mp.scatter(x_,y_,marker=mk,s=11,fc=COL[m] if id(r) in hit else 'none',ec=COL[m],lw=0.5,zorder=4)
mp.set_xlim(0,lov['w']*UM); mp.set_ylim(lov['h']*UM,0); mp.set_aspect('equal')
mp.set_xlabel('x (µm)'); mp.set_ylabel('y (µm)'); mp.set_xticks([0,100,200]); mp.set_yticks([0,100,200])
for s_ in ('top','right'): mp.spines[s_].set_visible(True)
mp.legend(handles=[Line2D([],[],marker='o',lw=0,ms=2.8,color=TP,label='lab TUNEL+'),Line2D([],[],marker='o',lw=0,ms=2.8,color=TN,label='lab TUNEL−'),
                   Line2D([],[],marker='s',lw=0,ms=2.8,mfc=COL['automated'],mec=COL['automated'],label='advanced'),Line2D([],[],marker='^',lw=0,ms=3,mfc=COL['fast'],mec=COL['fast'],label='fast'),
                   Line2D([],[],marker='s',lw=0,ms=2.8,mfc='none',mec='#555555',mew=0.5,label=f'open: no lab\nparasite within\n{MATCH_PX} px (4.8 µm)')],
          frameon=False,fontsize=4.5,loc='upper left',bbox_to_anchor=(1.02,1.02),handletextpad=0.2,borderaxespad=0,labelspacing=0.35)
idp=ax_mm(xc+61,yd+8,21,21)
idp.plot([0,10.5],[0,10.5],color='#999999',lw=0.4,ls=':',zorder=1)
for (k,pairs),mk,al in zip(H2.items(),('D','o'),(1,0.6)):
    idp.scatter([a_ for a_,_,_ in pairs],[b_ for _,b_,_ in pairs],marker=mk,s=7 if mk=='D' else 12,fc=COL['hil'] if mk=='D' else 'none',ec=COL['hil'],lw=0.5,alpha=al,zorder=3)
ea=err_auto('automated'); idp.scatter([a_ for a_,_ in ea],[b_ for _,b_ in ea],marker='s',s=9,color=COL['automated'],lw=0,zorder=4)
idp.set_xlim(2,10.5); idp.set_ylim(2,10.5); idp.set_xticks([2,4,6,8,10]); idp.set_yticks([2,4,6,8,10]); idp.set_aspect('equal')
idp.set_xlabel('Lab index (workbook)'); idp.set_ylabel('Run index')
idp.set_title('Per-parasite index',fontsize=FS,pad=2)
idp.legend(handles=[Line2D([],[],marker='D',lw=0,ms=2.5,color=COL['hil'],label='run 1'),Line2D([],[],marker='o',lw=0,ms=3,mfc='none',mec=COL['hil'],label='run 2'),
                    Line2D([],[],marker='s',lw=0,ms=2.8,color=COL['automated'],label='advanced')],frameon=False,fontsize=4.5,loc='upper left',handletextpad=0.2,borderaxespad=0.1,labelspacing=0.3)

# ================= e  task 2 = Fig. 4C =================
ye=yd+42.5; letter(xc,ye,'e'); head(xc+7,ye+0.6,'Task 2 result vs the published Fig. 4C')
SRC=['Published','+ human\nrun 1','+ human\nrun 2','advanced','fast']
def stars(p): return '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
rng=np.random.default_rng(0)
EH=15; d=ax_mm(xc+9,ye+6,52,EH); lim=(0,10.6)
for i,k in enumerate(K2):
    g=t2[k]
    for grp,dx,c in [('TUNEL+',-0.19,TP),('TUNEL-',0.19,TN),('Unassigned',0,UN)]:
        v=g.get(grp)
        if not v: continue
        d.scatter(i+dx+rng.uniform(-0.06,0.06,len(v)),np.clip(v,*lim),s=3,color=c,lw=0,alpha=0.85,zorder=2)
        d.plot([i+dx-0.14,i+dx+0.14],[st.median(v)]*2,color='k',lw=0.8,zorder=4,solid_capstyle='butt')
    if p2[i] is not None: d.text(i,lim[1],stars(p2[i]),ha='center',va='top',fontsize=FS)
    elif not g.get('TUNEL+'): d.text(i,0.6,'no TUNEL+\nparasite',ha='center',va='bottom',fontsize=4.5,color='#8a8a8a',style='italic',linespacing=0.95)
d.set_xticks(range(5)); d.set_xticklabels(SRC,fontsize=FS,linespacing=0.95)
for tl,m in zip(d.get_xticklabels(),MK): tl.set_color(COL[m])
d.tick_params(axis='x',length=0); d.set_ylim(*lim); d.set_yticks([0,5,10]); d.set_xlim(-0.55,4.55)
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
import matplotlib.transforms as mtr
ty=mtr.blended_transform_factory(d.transData,fig.transFigure); yrow=1-(ye+6+EH+6.2)/H; yci=1-(ye+6+EH+8.0)/H
d.text(-0.55,yrow,'Δ median  ',ha='right',va='center',fontsize=FS,fontweight='bold',transform=ty)
d.text(-0.55,yci,'95% CI  ',ha='right',va='center',fontsize=4.3,color='#555555',transform=ty)
for i,r in enumerate(DM):
    d.text(i,yrow,'–' if r is None else f'{r[0]:.2f}',ha='center',va='center',fontsize=FS,color=COL[MK[i]],transform=ty)
    if r is not None: d.text(i,yci,f'[{r[1]:.1f}, {r[2]:.1f}]',ha='center',va='center',fontsize=4.3,color='#555555',transform=ty)

for ext in ['pdf','svg','png']: fig.savefig(f'/work/Fig_Parasite.{ext}',dpi=450 if ext=='png' else None)
# ---- print what was plotted, so a rebuild can be checked against the figure
print('within',TOL,w2n); print('task1 calls agree:',agree,' r vs lab:',{m:round(v,3) for m,v in r1.items()})
for m in T1: print('  ',m,[f'{v:.2f}{"+" if g=="TUNEL+" else ""}' for v,g in T1[m]])
for k,v in H2.items(): print('task2 hil',k,'max|d|=%.4f'%e2[k],[(round(a_,3),round(b_,3),g) for a_,b_,g in v])
for m in A2: print('task2',m,'objects',len(A2[m]['objs']),'matched',[(g['name'],g['subset'],g['group'],r['group'],round(dd,1)) for g,r,dd in A2[m]['kept']],
                   'extra',len(A2[m]['extra']),'outside overview',nout[m],'control errors',err_auto(m) if m=='automated' else '')
for k,p,r in zip(K2,p2,DM): print(f'task2 {k:24}',{g:(len(v),round(st.median(v),3)) for g,v in t2[k].items()},'p=',p,'dmed/CI=',None if r is None else np.round(r,3))

