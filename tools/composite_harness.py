import sys; sys.path.insert(0,'tools')
import numpy as np
from make_channel_restart import natural_line, subdivide

# --- real channel metric at the y110 interface ---
ny=64; ly=2.0; y0=natural_line(ny,ly,16.0,0.5)      # level-0 nodes
y1=subdivide(y0)                                     # level-1 (fine) nodes
# bottom band = 3 block rows (nb=8) -> 24 level-0 cells refined; interface at level-0 node 24
i_if=24
hC=y0[i_if+1]-y0[i_if]            # coarse cell height just ABOVE the interface
hF=y1[2*i_if]-y1[2*i_if-1]        # fine cell height just BELOW the interface (level-1)
print(f"real metric @ interface (y110): coarse height hC={hC:.5e}, fine height hF={hF:.5e}, ratio={hC/hF:.3f}")
hx=hF                              # tangential fine width (uniform x): use hF as representative; coarse=2*hx

# --- composite 2D (x-y) patch: periodic x, Neumann y. Coarse rows above, fine rows below. ---
Ncx=3                              # coarse cells in x  (fine = 2*Ncx)
Nfx=2*Ncx
NcyA=3                             # coarse rows above interface
NfyB=3                             # fine rows below interface
# real cell heights: take a few coarse heights above and fine heights below the interface
hc=[y0[i_if+r+1]-y0[i_if+r] for r in range(NcyA)]                 # coarse rows upward
hf=[y1[2*i_if-1-r]-y1[2*i_if-2-r] for r in range(NfyB)]           # fine rows downward (top fine row first)

# cell ids
cells={}; cid=0
for r in range(NfyB):        # fine rows below (r=0 nearest interface)
    for c in range(Nfx): cells[('F',r,c)]=cid; cid+=1
for r in range(NcyA):        # coarse rows above (r=0 nearest interface)
    for c in range(Ncx): cells[('C',r,c)]=cid; cid+=1
ncell=cid
Vc=np.zeros(ncell)
for (t,r,c),idx in cells.items():
    Vc[idx]=(hx*hf[r]) if t=='F' else (2*hx*hc[r])

# face DOFs: build D as dict rows. faces: ('u',...) periodic x ; ('v',...) y-normal.
faces={}; fid=0
def face(key):
    global fid
    if key not in faces: faces[key]=fid; fid+=1
    return faces[key]
rows=[]  # (cell, face, coeff)  ; coeff = +-area/Vcell
# fine region u-faces (periodic in x), heights hf[r], width hx
for r in range(NfyB):
    for c in range(Nfx):
        cl=cells[('F',r,c)]; A=hf[r]
        fL=face(('uF',r,c)); fR=face(('uF',r,(c+1)%Nfx))
        rows.append((cl,fL,-A/Vc[cl])); rows.append((cl,fR,+A/Vc[cl]))
# coarse region u-faces
for r in range(NcyA):
    for c in range(Ncx):
        cl=cells[('C',r,c)]; A=hc[r]
        fL=face(('uC',r,c)); fR=face(('uC',r,(c+1)%Ncx))
        rows.append((cl,fL,-A/Vc[cl])); rows.append((cl,fR,+A/Vc[cl]))
# fine-fine internal v-faces (between fine row r and r+1), area hx ; top of patch fine handled by interface
for r in range(NfyB-1):
    for c in range(Nfx):
        cu=cells[('F',r,c)]; cd=cells[('F',r+1,c)]; A=hx
        f=face(('vF',r,c))   # face between fine row r (above) and r+1 (below)
        rows.append((cu,f,-A/Vc[cu]))   # bottom face of upper cell
        rows.append((cd,f,+A/Vc[cd]))   # top face of lower cell
# coarse-coarse internal v-faces
for r in range(NcyA-1):
    for c in range(Ncx):
        cu=cells[('C',r+1,c)]; cd=cells[('C',r,c)]; A=2*hx
        f=face(('vC',r,c))
        rows.append((cu,f,-A/Vc[cu])); rows.append((cd,f,+A/Vc[cd]))
# INTERFACE v-faces: the Nfx fine faces (authoritative). Fine top row r=0; coarse bottom row r=0.
for c in range(Nfx):
    fF=cells[('F',0,c)]; cc=c//2; fC=cells[('C',0,cc)]; A=hx
    f=face(('vIF',c))
    rows.append((fF,f,+A/Vc[fF]))      # fine cell top face (outward +)
    rows.append((fC,f,-A/Vc[fC]))      # coarse cell bottom face: SUM of the 2 fine faces (each area hx)
# (far top of coarse, far bottom of fine = Neumann v=0: not DOFs)
nface=fid
D=np.zeros((ncell,nface))
for cl,f,co in rows: D[cl,f]+=co
# face control volumes Wf = area * normal spacing(center-to-center)
Wf=np.zeros(nface)
for key,f in faces.items():
    if key[0] in ('uF','uC'):
        Wf[f]= (hf[key[1]] if key[0]=='uF' else hc[key[1]]) * (hx if key[0]=='uF' else 2*hx)
    elif key[0]=='vF':
        r=key[1]; Wf[f]=hx*(hf[r]/2+hf[r+1]/2)
    elif key[0]=='vC':
        r=key[1]; Wf[f]=2*hx*(hc[r]/2+hc[r+1]/2)
    elif key[0]=='vIF':
        Wf[f]=hx*(hf[0]/2+hc[0]/2)     # fine-below to coarse-above center gap
Wc=Vc.copy()
# G = -Dᵀ in weighted inner product:  G = -Wf^{-1} Dᵀ Wc ;  L = D G
G=-(np.diag(1/Wf) @ D.T @ np.diag(Wc))
L=D@G
# checks
WcL=np.diag(Wc)@L
asym=np.abs(WcL-WcL.T).max()
ev=np.linalg.eigvalsh(0.5*(WcL+WcL.T))   # generalized symmetric part
# adjointness residual D + Gᵀ = 0 in weighted ip:  Wc D + (Wf G)ᵀ = 0
adj=np.abs(np.diag(Wc)@D + (np.diag(Wf)@G).T).max()
print(f"\nDOFs: {ncell} cells, {nface} faces")
print(f"adjointness  max|Wc D + (Wf G)^T| = {adj:.3e}")
print(f"L symmetry   max|WcL - (WcL)^T|   = {asym:.3e}")
print(f"L eigenvalues (Wc-sym): max={ev.max():.3e}  min={ev.min():.3e}")
nnull=int(np.sum(np.abs(ev)<1e-9*max(1,abs(ev.min()))))
print(f"near-zero eigenvalues (null space dim): {nnull}   (expect 1 = constant pressure)")
# interface gradient consistency vs (p_C - p_F)/(hf/2+hc/2)
fIF=faces[('vIF',0)]; fF=cells[('F',0,0)]; fC=cells[('C',0,0)]
print(f"\ninterface gradient G[v_if, p_C]={G[fIF,fC]:+.5e}  G[v_if, p_F]={G[fIF,fF]:+.5e}")
print(f"expected 1/(hf/2+hc/2) = {1/(hf[0]/2+hc[0]/2):+.5e}  (= (2/3)/h at uniform)")
