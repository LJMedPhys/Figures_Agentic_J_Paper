"""Score nucleus masks of every method against the expert's nucleus annotations (11 clean tiles,
lobes regrouped at GAP_PX=15, rerun_lobes_merged/labels) with the same metric as the expert's CV
(elf mean_segmentation_accuracy). HIL = held-out 4-fold CV predictions (model never saw the tile)."""
import glob, json, os, numpy as np, pandas as pd, tifffile, warnings, logging
from scipy import ndimage as ndi
from elf.evaluation import mean_segmentation_accuracy
warnings.filterwarnings('ignore'); logging.getLogger('tifffile').setLevel(logging.ERROR)
FT='/ft'; RR=f'{FT}/rerun_lobes_merged'; GAP_PX,MIN_FRAG,MIN_OBJ=15,10,25
def merge_lobes(gt):  # copied from rerun_cv.py
    frags=[]
    for i in np.unique(gt):
        if i==0: continue
        lab,n=ndi.label(gt==i,structure=np.ones((3,3)))
        frags+=[lab==k for k in range(1,n+1) if (lab==k).sum()>=MIN_FRAG]
    parent=list(range(len(frags)))
    def find(a):
        while parent[a]!=a: a=parent[a]
        return a
    dts=[ndi.distance_transform_edt(~m) for m in frags]
    for a in range(len(frags)):
        for b in range(a+1,len(frags)):
            if dts[a][frags[b]].min()<=GAP_PX: parent[find(a)]=find(b)
    out=np.zeros(gt.shape,np.uint32); nxt=1
    for root in sorted({find(a) for a in range(len(frags))}):
        m=np.any([frags[a] for a in range(len(frags)) if find(a)==root],axis=0)
        if m.sum()>=MIN_OBJ: out[m]=nxt; nxt+=1
    return out
def per_cell_nuclei(cells,nucmask):
    out=np.zeros(cells.shape,np.uint32)
    out[nucmask&(cells>0)]=cells[nucmask&(cells>0)]
    lab=np.zeros_like(out); k=1          # drop tiny slivers, relabel
    for v in np.unique(out)[1:]:
        m=out==v
        if m.sum()>=MIN_OBJ: lab[m]=k; k+=1
    return lab
tiles={t['name']:t for t in json.load(open(f'{FT}/manifest.json'))['tiles']}
hilpred={os.path.basename(p)[:-9]:p for p in glob.glob(f'{RR}/jobs/vit_b_lm_fold*/preds/*_best.tif')}
stockpred={os.path.basename(p)[:-10]:p for p in glob.glob(f'{RR}/jobs/vit_b_lm_fold*/preds/*_stock.tif')}
rows=[]
for gp in sorted(glob.glob(f'{RR}/labels/*.tif')):
    t=os.path.basename(gp)[:-4]; gt=tifffile.imread(gp); m=tiles[t]; s=os.path.splitext(m['source'])[0]
    sl=(slice(m['y0'],m['y0']+m['height']),slice(m['x0'],m['x0']+m['width']))
    inst=tifffile.imread(f'/ref/full/instance_segmentation/{s}.TIF').astype(np.int64)
    dil=np.zeros_like(inst)
    for lab in np.unique(inst)[1:]:
        dil[ndi.binary_dilation(inst==lab,np.ones((3,3)),iterations=3)]=lab
    probf=[p for p in [f'/ref/full/segmentation/{s}_Probabilities.tiff',f'/ref/full/segmentation/New Folder With Items/{s}_Probabilities.tiff'] if os.path.exists(p)][0]
    prob=tifffile.imread(probf)
    preds={'published':per_cell_nuclei(dil[sl],prob[...,0][sl]>0.5),
           'automated':per_cell_nuclei(tifffile.imread(f'/auto/{s}_cell_labels.tif')[sl].astype(np.int64),tifffile.imread(f'/auto/{s}_nucleus_labels.tif')[sl]>0),
           'fast':(lambda a: per_cell_nuclei(a[0][sl].astype(np.int64),a[1][sl]>0))(tifffile.imread(f'/fast/{s}_labels.tif')),
           'codex':(lambda a: per_cell_nuclei(a[0][sl].astype(np.int64),a[1][sl]>0))(tifffile.imread(f'/codex/{s}_labels.tif')),
           'hil_stock':merge_lobes(tifffile.imread(stockpred[t])),
           'hil':merge_lobes(tifffile.imread(hilpred[t]))}
    for k,p in preds.items():
        msa,acc=mean_segmentation_accuracy(p,gt,return_accuracies=True)
        rows.append(dict(tile=t,group=m['group'],method=k,msa=float(msa),sa50=float(acc[0]),n_pred=int(p.max()),n_true=int(gt.max()),**{f'acc_{th:.2f}':float(a) for th,a in zip(np.arange(0.5,1.0,0.05),acc)}))
df=pd.DataFrame(rows); df.to_csv('/work/nuclei_vs_expert.csv',index=False)
print(df.groupby('method')[['msa','sa50']].mean().round(3))
print(df.groupby('method')[[c for c in df if c.startswith('acc_')]].mean().round(3).T)
