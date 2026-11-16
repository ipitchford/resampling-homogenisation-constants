from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))  # sibling verify.py (needed under python -I)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from verify import occupancy, poly_modulus
BASE=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#322c38','text.color':'#322c38','axes.edgecolor':'#857e89','savefig.facecolor':'white'})
purple='#68416D';ochre='#B57A22';teal='#2A7B78';coral='#BB6559'
fig,axs=plt.subplots(3,1,figsize=(6.15,8.0),layout='constrained')
ns=np.arange(2,31);cs=[float(sum(w*max(k,l) for k,l,w in occupancy(int(n)))) for n in ns]
axs[0].plot(ns,cs,color=purple,lw=1.7);axs[0].scatter([2,3,4],[cs[0],cs[1],cs[2]],color=[teal,ochre,purple],s=24,zorder=4)
axs[0].axhline(14/9,color=ochre,ls=':',lw=1);axs[0].annotate('14/9 at n = 3',(3,14/9),xytext=(9,1.552),color=ochre,arrowprops={'arrowstyle':'-','color':ochre})
axs[0].set(xlabel='Number of components n',ylabel='Sharp amplification factor',title='A  Uniform sampling with replacement',xlim=(1.5,30.5),ylim=(1.496,1.561))
ds=np.linspace(.0001,.9999,400);vals=np.array([poly_modulus(float(d)) for d in ds]);f3=14*ds/9-2*ds**2/3+2*ds**3/27
axs[1].plot(ds,vals[:,0],color=purple,lw=2,label='All asymmetric erasure pairs')
axs[1].plot(ds,f3,color=coral,lw=1.5,ls='--',label='Only two differing pairs (n = 3)')
axs[1].plot(ds,ds,color='#999298',ls=':',lw=1,label='H = D')
d0=.3248691294333441;d1=(5-np.sqrt(13))/2;dstar=(9-np.sqrt(33))/8
for d,lab in [(d0,'D₀'),(d1,'D₁')]:axs[1].axvline(d,color=ochre,lw=1,ls=':');axs[1].text(d,.03,lab,color=ochre,ha='center')
axs[1].axvline(dstar,color=coral,lw=1,ls=':');axs[1].text(dstar,.93,'D*',color=coral,ha='center')
axs[1].set(xlabel='Original total variation D',ylabel='Largest homogenised TV, H',title='B  Exact optimum for n = 3',xlim=(0,1),ylim=(0,1));axs[1].legend(frameon=False,fontsize=8,loc='lower right')
axs[2].plot(ds,vals[:,1],color=purple,lw=1.8,label='x*');axs[2].plot(ds,vals[:,2],color=teal,lw=1.8,label='y*')
d0=.3248691294333441;d1=(5-np.sqrt(13))/2
for d in [d0,d1]:axs[2].axvline(d,color=ochre,lw=1,ls=':')
axs[2].text(d0,.49,'D₀',color=ochre,ha='center');axs[2].text(d1,.49,'D₁',color=ochre,ha='center')
axs[2].set(xlabel='Original total variation D',ylabel='Optimal logarithmic allocation',title='C  Two changes in the extremiser',xlim=(0,1),ylim=(-.02,1.02));axs[2].legend(frameon=False,loc='center left')
for ax in axs:ax.grid(axis='y',alpha=.14)
for ext in ['pdf','png','svg']:fig.savefig(BASE/'figures'/('sharp_homogenisation.'+ext),dpi=220)
print('Saved scientific figure in PDF, PNG and SVG.')
