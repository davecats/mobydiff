import sys, glob; sys.path.insert(0,'../assets/postpro')
import rans_bl as r, numpy as np
d=r.dns(); X=(15,20,30,40,50,60,70,80,100,110)
print('x      ', ' '.join(f'{x:7.0f}' for x in X))
j=[np.argmin(abs(d['x']-x)) for x in X]
print('DNS cf ', ' '.join(f'{d["cf"][k]:7.5f}' for k in j)); print('DNS Rth', ' '.join(f'{d["Re_theta"][k]:7.0f}' for k in j))
for c in sys.argv[1:]:
  fs=sorted(glob.glob(c+'_*.h5'),key=lambda s:int(s.split('_')[-1][:-3])); f=r.load(fs[-1]); R=r.integrals(f['x'],f['y'],f['U'],450)
  i=[np.argmin(abs(f['x']-x)) for x in X]
  print(c,fs[-1]); print(' cf    ', ' '.join(f'{R["cf"][k]:7.5f}' for k in i)); print(' Rth   ', ' '.join(f'{R["Re_theta"][k]:7.0f}' for k in i))
