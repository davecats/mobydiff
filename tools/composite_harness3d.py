# P-b.0: composite D / G=-D^T / L=-D D^T adjointness + SPD harness on 3D layouts.
# General builder: a coarse base grid; a refinement mask marks coarse cells that
# split into 2x2x2 fine cells. Builds the conservative control-volume divergence
# D (coarse cell facing finer neighbour sums the 4 covering fine faces), the
# weighted transpose G=-D^T, L=DG, and checks D+G^T=0, L symmetric NSD, null=1.
import numpy as np, itertools

def build(refined, NC, h=1.0, periodic=(True,True,True)):
    # coarse cell size h_c=2h (so fine=h). refined: set of (I,J,K). NC=(Nx,Ny,Nz) coarse.
    Nx,Ny,Nz=NC; hc=2*h; hf=h
    # ---- enumerate pressure cells ----
    cell={}; cid=0
    def addcell(key):
        nonlocal cid
        cell[key]=cid; cid+=1
    for I in range(Nx):
        for J in range(Ny):
            for K in range(Nz):
                if (I,J,K) in refined:
                    for a,b,c in itertools.product((0,1),repeat=3): addcell(('F',I,J,K,a,b,c))
                else: addcell(('C',I,J,K))
    ncell=cid
    Vc=np.zeros(ncell)
    for key,idx in cell.items():
        Vc[idx]=hf**3 if key[0]=='F' else hc**3
    # ---- faces: iterate every cell, its +d faces; register a DOF; coarse facing
    #      finer neighbour references the 4 fine faces (no own DOF). ----
    faces={}; fid=0
    def face(key):
        nonlocal fid
        if key not in faces: faces[key]=fid; fid+=1
        return faces[key]
    rows=[]   # (cell, face, signed area/Vcell)
    Wf={}
    def wrap(I,d,N): return (I% N) if periodic[d] else I
    def is_ref(I,J,K):
        if not(0<=I<Nx and 0<=J<Ny and 0<=K<Nz):
            if not periodic: return None
            I,J,K=wrap(I,0,Nx),wrap(J,1,Ny),wrap(K,2,Nz)
        return (I,J,K) in refined
    # helper: fine subface key on the boundary between coarse cell (I,J,K) and its
    # +d refined neighbour, fine sub-index (s,t) in the tangential plane.
    def add_interface(dirn, cI,cJ,cK, nI,nJ,nK):
        # coarse cell (cI..) has its +dirn face covered by 4 fine faces of refined (nI..)
        d=dirn
        ccell=cell[('C',cI,cJ,cK)]
        for s,t in itertools.product((0,1),repeat=2):
            # fine cell in neighbour on the touching layer
            sub=[0,0,0]
            if d==0: sub=[0,s,t]                 # +x face: fine a=0 layer
            elif d==1: sub=[s,0,t]
            else: sub=[s,t,0]
            fk=('F',nI%Nx,nJ%Ny,nK%Nz,*sub)
            fcell=cell[fk]
            fkey=('IF',d,cI,cJ,cK,s,t)
            f=face(fkey); A=hf*hf
            rows.append((ccell,f,+A/Vc[ccell]))   # coarse +d outward
            rows.append((fcell,f,-A/Vc[fcell]))    # fine -d (its low face) inward
            # control volume: area * center gap (coarse center to fine center)
            Wf[f]=A*(hc/2+hf/2)
    # main face loop: for each cell, the +x,+y,+z face
    for key,idx in list(cell.items()):
        if key[0]=='C':
            I,J,K=key[1:]
            for d,(dI,dJ,dK) in enumerate([(1,0,0),(0,1,0),(0,0,1)]):
                nb=is_ref(I+dI,J+dJ,K+dK)
                if nb is None: continue            # non-periodic boundary: Neumann (no DOF)
                if nb:   # neighbour refined -> 4 fine subfaces
                    add_interface(d,I,J,K,I+dI,J+dJ,K+dK)
                else:    # coarse-coarse face
                    nc=cell[('C',(I+dI)%Nx,(J+dJ)%Ny,(K+dK)%Nz)]
                    f=face(('CC',d,I,J,K)); A=hc*hc
                    rows.append((idx,f,+A/Vc[idx])); rows.append((nc,f,-A/Vc[nc]))
                    Wf[f]=A*hc
        else: # fine cell
            I,J,K,a,b,c=key[1:]; sub=(a,b,c)
            for d,(dI,dJ,dK) in enumerate([(1,0,0),(0,1,0),(0,0,1)]):
                s=sub[d]
                if s==0:  # internal-to-block +d face goes to sub=1 ; handled by the s==1 cell's -d? 
                    nsub=list(sub); nsub[d]=1
                    nc=cell[('F',I,J,K,*nsub)]
                    f=face(('FFi',d,I,J,K,sub[(d+1)%3],sub[(d+2)%3])); A=hf*hf
                    rows.append((idx,f,+A/Vc[idx])); rows.append((nc,f,-A/Vc[nc])); Wf[f]=A*hf
                else:     # s==1: +d face leaves the block to neighbour coarse cell at I+dI
                    nb=is_ref(I+dI,J+dJ,K+dK)
                    if nb is None: continue
                    if nb:  # neighbour also refined -> fine-fine across blocks
                        nsub=list(sub); nsub[d]=0
                        nc=cell[('F',(I+dI)%Nx,(J+dJ)%Ny,(K+dK)%Nz,*nsub)]
                        f=face(('FFx',d,I,J,K,sub[(d+1)%3],sub[(d+2)%3])); A=hf*hf
                        rows.append((idx,f,+A/Vc[idx])); rows.append((nc,f,-A/Vc[nc])); Wf[f]=A*hf
                    # else neighbour coarse: this fine face is part of the coarse cell's
                    # interface, already added by the COARSE side (add_interface). skip.
    nface=fid
    D=np.zeros((ncell,nface))
    for cl,f,co in rows: D[cl,f]+=co
    Wfv=np.array([Wf[i] for i in range(nface)])
    return D,Vc,Wfv

def check(name,D,Vc,Wf):
    G=-(np.diag(1/Wf)@D.T@np.diag(Vc)); L=D@G; WcL=np.diag(Vc)@L
    adj=np.abs(np.diag(Vc)@D+(np.diag(Wf)@G).T).max()
    asym=np.abs(WcL-WcL.T).max(); ev=np.linalg.eigvalsh(0.5*(WcL+WcL.T))
    nnull=int(np.sum(np.abs(ev)<1e-9*max(1,abs(ev.min()))))
    # conservation: sum_c Vc*(D u) for random u = boundary flux; for all-periodic = 0
    print(f" {name:28s} cells={len(Vc):4d} faces={D.shape[1]:4d}  adj={adj:.2e}  Lsym={asym:.2e}  ev[max,min]=[{ev.max():+.2e},{ev.min():+.2e}]  null={nnull}")

# (i) 3D PLANAR y-interface: fine slab (one coarse y-layer refined), periodic x,z
NC=(2,3,2); ref={(I,1,K) for I in range(2) for K in range(2)}
D,Vc,Wf=build(ref,NC,periodic=(True,True,True)); check("3D planar slab (4 subfaces)",D,Vc,Wf)
# (ii) embedded fine CUBE: single refined coarse cell -> faces+12 edges+8 corners
NC=(3,3,3); ref={(1,1,1)}
D,Vc,Wf=build(ref,NC,periodic=(True,True,True)); check("3D embedded cube (edge+corner)",D,Vc,Wf)
# (iii) refined cell touching the periodic boundary
NC=(3,3,3); ref={(0,1,1)}
D,Vc,Wf=build(ref,NC,periodic=(True,True,True)); check("refined @ periodic boundary",D,Vc,Wf)
# (iv) two adjacent refined cells (shared fine-fine across-block + edges)
NC=(3,3,3); ref={(1,1,1),(1,1,2)}
D,Vc,Wf=build(ref,NC,periodic=(True,True,True)); check("two adjacent refined cells",D,Vc,Wf)

# Explicit conservation: for an all-periodic patch, sum_c Vc*(D u) = 0 for ANY u
# (the control-volume divergence telescopes; the coarse interface flux is exactly
# the area-sum of the fine faces, so no double counting). Instrument it.
print("\nConservation  sum_c Vc*(D u) for random u (all-periodic -> 0 exactly):")
for name,NC,ref in [("planar slab",(2,3,2),{(I,1,K) for I in range(2) for K in range(2)}),
                    ("embedded cube",(3,3,3),{(1,1,1)}),
                    ("two adjacent",(3,3,3),{(1,1,1),(1,1,2)})]:
    D,Vc,Wf=build(ref,NC,periodic=(True,True,True))
    rng=np.arange(1,D.shape[1]+1,dtype=float); u=np.sin(0.7*rng)+0.3*np.cos(1.9*rng)
    print(f"  {name:16s}  |sum Vc*(D u)| = {abs(Vc@(D@u)):.3e}")
