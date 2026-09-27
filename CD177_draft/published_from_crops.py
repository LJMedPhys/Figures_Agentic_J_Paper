"""Per-cell published measurements from the paper's own crops (full/extracted_cells/*_semantic.tiff),
with the paper's s3_extract_feats.py definitions: cell = sum(prob) > 0.1, nucleus = prob[...,0] > 0.5,
n/c = nucleus px / cell px (G_nuc_to_cyto)."""
import glob, os, numpy as np, pandas as pd, tifffile, warnings, logging
from multiprocessing import Pool
warnings.filterwarnings('ignore'); logging.getLogger('tifffile').setLevel(logging.ERROR)
def one(p):
    seg=tifffile.imread(p); cell=np.sum(seg,axis=-1)>0.1; nuc=seg[...,0]>0.5
    b=os.path.basename(p)[:-len('_semantic.tiff')]
    return dict(crop=b,cond=os.path.basename(os.path.dirname(p)),is60x='60x' in b,cell_px=int(cell.sum()),nuc_px=int(nuc.sum()),nc=nuc.sum()/(cell.sum()+1e-8))
if __name__=='__main__':
    fs=sorted(glob.glob('/ref/full/extracted_cells/*/*_semantic.tiff'))
    with Pool(16) as pl: df=pd.DataFrame(pl.map(one,fs))
    df.to_csv('/work/published_cells_from_crops.csv',index=False)
    print(len(df),'crops;',df.is60x.sum(),'from 60x'); d=df[~df.is60x]
    print(d.groupby('cond')[['cell_px','nuc_px','nc']].median().round(3))
    from scipy.stats import mannwhitneyu
    for k in ['cell_px','nc','nuc_px']:
        a,b=d[d.cond=='pos'][k],d[d.cond=='neg'][k]; U,p=mannwhitneyu(a,b); print(k,'MWU p=%.2g'%p,'r=%.3f'%abs(1-2*U/(len(a)*len(b))))
