import csv, sys, matplotlib
from matplotlib.lines import Line2D
matplotlib.use('Agg'); import matplotlib.pyplot as plt
# Fig. 4: dot-line plots of Pearson r (95% bootstrap CI) between IAKA and IHKA/PILO log2FC vectors.
# Usage: python Figure4_dotline.py <07.CrossDataset dir> <output prefix>
# Inputs: 03.DEGCorrelation/Pearson_correlation_summary.csv,
#         03.DEGCorrelation/Pearson_correlation_summary_EpiGenes.csv,
#         04.IHKAByDataset/Pearson_correlation_by_dataset.csv
E=sys.argv[1].rstrip('/')+'/'; OUT=sys.argv[2]
def load(f,key):
    d={}
    for r in csv.DictReader(open(f)):
        d[key(r)]=(float(r['r']),float(r['ci_lo']),float(r['ci_hi']))
    return d
W=load(E+'03.DEGCorrelation/Pearson_correlation_summary.csv',lambda r:(r['model'],r['iaka'],r['comparison']))
G=load(E+'03.DEGCorrelation/Pearson_correlation_summary_EpiGenes.csv',lambda r:(r['model'],r['iaka'],r['comparison']))
D={}
for r in csv.DictReader(open(E+'04.IHKAByDataset/Pearson_correlation_by_dataset.csv')):
    gs='W' if r['gene_set'].startswith('Whole') else 'G'
    D[(gs,r['dataset'].split()[0],r['iaka'],r['comparison'])]=(float(r['r']),float(r['ci_lo']),float(r['ci_hi']))
P=['AC','IM','CR']; C=['HA','AC','IM','CR']; CL=['Hyperacute','Acute','Intermediate','Chronic']
col={'AC':'#CCCCCC','IM':'#888888','CR':'#404040'}
RED='#A30234'
mk={'AC':'o','IM':'s','CR':'^'}
ms={'AC':4.3,'IM':4.0,'CR':5.4}
lab={'AC':'IAKA acute (day 3)','IM':'IAKA intermediate (day 7)','CR':'IAKA chronic (day 14)'}
plt.rcParams.update({'font.family':['Arial','Liberation Sans','DejaVu Sans'],'font.size':7,'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6})
fig=plt.figure(figsize=(7.2,5.4))
gsp=fig.add_gridspec(2,5,width_ratios=[1,1,0.28,1,1],wspace=0.1,hspace=0.62,left=0.075,right=0.99,top=0.9,bottom=0.15)
def panel(ax,get,title,cols,showy):
    for k,p in enumerate(P):
        xs=[]; v=[]
        for j,c in enumerate(cols):
            t=get(p,c)
            if t is None: continue
            xs.append(j+(k-1)*0.13); v.append((t,c))
        ax.plot(xs,[a[0][0] for a in v],'-',color=col[p],lw=1.3,zorder=1)
        ax.errorbar(xs,[a[0][0] for a in v],yerr=[[a[0][0]-a[0][1] for a in v],[a[0][2]-a[0][0] for a in v]],fmt=mk[p],color=RED,mfc=RED,mec='white',mew=0.4,ms=ms[p],capsize=0,elinewidth=0.8,ecolor=RED,zorder=2,label=lab[p])
        for x,(t,c) in zip(xs,v):
            if c==p:
                others=[get(q,c) for q in P if q!=p and get(q,c) is not None]
                below=any(0<o[0]-t[0]<0.09 for o in others)
                off=(0,-11) if below else (0,7)
                ax.scatter(x,t[0],s=58,facecolors='none',edgecolors='black',lw=0.7,zorder=3)
                ax.annotate(f'{t[0]:.2f}',(x,t[0]),xytext=off,textcoords='offset points',ha='center',va='center',fontsize=5.8)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([CL[C.index(c)] for c in cols],rotation=35,ha='right')
    ax.set_xlim(-0.5,len(cols)-0.5); ax.set_ylim(0,1.0); ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
    ax.grid(axis='y',color='#e6e6e6',lw=0.5); ax.spines[['top','right']].set_visible(False)
    ax.set_title(title,fontweight='bold',fontsize=7.5,pad=4)
    if showy: ax.set_ylabel('Pearson r (95% CI)',labelpad=2)
    else: ax.tick_params(labelleft=False)
axs={}
spec=[('a',0,0,'Whole transcriptome'),('b',0,3,'Epileptogenesis-related genes'),('c',1,0,'Whole transcriptome'),('d',1,3,'Epileptogenesis-related genes')]
for L,row,c0,sub in spec:
    a1=fig.add_subplot(gsp[row,c0]); a2=fig.add_subplot(gsp[row,c0+1],sharey=a1)
    if row==0:
        S=W if L=='a' else G
        panel(a1,lambda p,c,S=S:S.get(('IHKA',p,c)),'IHKA',C,True)
        panel(a2,lambda p,c,S=S:S.get(('PILO',p,c)),'PILO',C,False)
    else:
        gs='W' if L=='c' else 'G'
        panel(a1,lambda p,c,gs=gs:D.get((gs,'CA1',p,c)),'CA1 (GSE99577)',C,True)
        panel(a2,lambda p,c,gs=gs:D.get((gs,'DG',p,c)),'DG (GSE213393)',C[1:],False)
    pos=a1.get_position()
    fig.text(pos.x0-0.055,pos.y1+0.045,L,fontsize=11,fontweight='bold')
    fig.text((pos.x0+a2.get_position().x1)/2,pos.y1+0.045,sub,ha='center',fontsize=7.5)
    axs[L]=a1
h=[Line2D([],[],color=col[p],lw=1.3,marker=mk[p],mfc=RED,mec='white',ms=ms[p]+1) for p in P]; l=[lab[p] for p in P]
h.append(Line2D([],[],marker='o',ls='',mfc='none',mec='black',ms=6)); l.append('Phase-matched comparison')
fig.legend(h,l,loc='lower center',ncol=4,frameon=False,fontsize=7,bbox_to_anchor=(0.53,-0.01))
fig.savefig(OUT+'.png',dpi=600); fig.savefig(OUT+'.svg'); fig.savefig(OUT+'.pdf')
