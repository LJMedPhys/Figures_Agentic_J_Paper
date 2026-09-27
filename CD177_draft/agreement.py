"""Score existing masks of each run against the published (Cellpose + ilastik) masks.
Reference cell = paper Cellpose instance (DILATE=1: dilated 3x(3x3) as in s2_process_all_cells.py; DILATE=0: undilated), kept if 50<=bbox w,h<=120;
reference nucleus = ilastik nucleus prob > 0.5 inside that cell (s3_extract_feats.py)."""
import glob, os, numpy as np, pandas as pd, tifffile, cv2, warnings
from multiprocessing import Pool
from scipy import ndimage as ndi
warnings.filterwarnings('ignore'); import logging; logging.getLogger('tifffile').setLevel(logging.ERROR)
R='/ref/full'
DILATE=os.environ.get('DILATE','1')=='1'   # 1 = paper's measured (dilated) cell, 0 = raw Cellpose instance
OUT='/work/inputs/mask_agreement_per_cell.csv' if DILATE else '/work/inputs/mask_agreement_per_cell_undilated.csv'
def prob_path(s):
    for d in ['segmentation','segmentation/New Folder With Items']:
        p=f'{R}/{d}/{s}_Probabilities.tiff'
        if os.path.exists(p): return p
auto_tab=pd.read_csv('/autotab/cd177_cell_features.csv')
auto_meas=set(zip(auto_tab.image.map(lambda x: os.path.splitext(os.path.basename(x))[0]), auto_tab.image_cell_label)) if 'size_analysis_eligible' not in auto_tab else \
          set(zip(auto_tab[auto_tab.size_analysis_eligible.astype(str)=='True'].image.map(lambda x: os.path.splitext(os.path.basename(x))[0]), auto_tab[auto_tab.size_analysis_eligible.astype(str)=='True'].image_cell_label))
fast_tab=pd.read_csv('/fasttab/cell_features.csv')
fast_meas=set(zip(fast_tab.image.map(lambda x: os.path.splitext(os.path.basename(x))[0]), fast_tab.cell_label))
hil_tab=pd.read_csv('/hiltab/per_cell_morphology_eligible.csv')
hil_meas=set(zip(hil_tab.basename.map(lambda x: os.path.splitext(x)[0]), hil_tab.original_cell_label))
def load_method(m,s):
    if m=='automated':
        c=tifffile.imread(f'/auto/{s}_cell_labels.tif'); n=tifffile.imread(f'/auto/{s}_nucleus_labels.tif')>0
    elif m=='fast':
        a=tifffile.imread(f'/fast/{s}_labels.tif'); c,n=a[0],a[1]>0
    elif m=='codex':
        a=tifffile.imread(f'/codex/{s}_labels.tif'); c,n=a[0],a[1]>0
    else:
        c=tifffile.imread(f'/hil/cells/{s}_cells.tif'); n=tifffile.imread(f'/hil/lobes/{s}_lobes.tif')>0
    return c.astype(np.int64),n
def one(s):
    rows=[]
    inst=tifffile.imread(f'{R}/instance_segmentation/{s}.TIF').astype(np.int64)
    prob=tifffile.imread(prob_path(s)); refnuc=prob[...,0]>0.5
    refs=[]
    for lab,sl in enumerate(ndi.find_objects(inst),1):
        if sl is None: continue
        h=sl[0].stop-sl[0].start; w=sl[1].stop-sl[1].start
        if w<50 or h<50 or w>120 or h>120: continue
        b=(inst==lab).astype(np.uint8)
        refs.append((lab,cv2.dilate(b,np.ones((3,3),np.uint8),iterations=3)>0 if DILATE else b>0))
    for m in ['automated','fast','codex','hil']:
        try: pc,pn=load_method(m,s)
        except FileNotFoundError: continue
        meas = auto_meas if m=='automated' else hil_meas if m=='hil' else fast_meas if m=='fast' else None
        npred=len(np.unique(pc))-1
        for lab,rc in refs:
            ov=pc[rc]; ov=ov[ov>0]
            row=dict(image=s,method=m,ref_label=lab,ref_cell_px=int(rc.sum()),ref_nuc_px=int((refnuc&rc).sum()),n_pred_cells_image=npred)
            if ov.size:
                best=np.bincount(ov).argmax(); pcell=pc==best
                inter=(pcell&rc).sum(); iou=inter/(pcell|rc).sum()
                a=pn&pcell; b=refnuc&rc; d=2*(a&b).sum()/max(a.sum()+b.sum(),1)
                row.update(pred_label=int(best),cell_iou=iou,nuc_dice=d,pred_cell_px=int(pcell.sum()),pred_nuc_px=int(a.sum()),
                           measured=(True if meas is None else (s,int(best)) in meas))
            else: row.update(cell_iou=0.0)
            rows.append(row)
    return rows
if __name__=='__main__':
    stems=[os.path.splitext(os.path.basename(f))[0] for f in glob.glob(f'{R}/raw/*.TIF') if '60x' not in f]
    with Pool(16) as p: out=p.map(one,sorted(stems))
    df=pd.DataFrame([r for rs in out for r in rs]); df.to_csv(OUT,index=False)
    df['matched']=df.cell_iou>=0.5; df['batch']=np.where(df.image.str.contains('Top'),'Top','Bottom')
    print(df.groupby(['method','batch']).agg(n_ref=('matched','size'),recall=('matched','mean'),
          iou_med=('cell_iou',lambda x:x[x>=0.5].median()),nuc_dice_med=('nuc_dice',lambda x: x[df.loc[x.index,'matched']].median())).round(3))
