"""Read existing 510-minute snapshots and plot death/growth and radial profiles."""
from pathlib import Path
import sys, os, csv, json
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
from compare_biological_outputs import read, scalar
os.environ.setdefault('MPLCONFIGDIR', '/tmp/physicell-comparison-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
rows = []
for col, (label, folder) in enumerate([('PhysiCell-X', ROOT/'output'), ('PhysiCelldFBA', ROOT.parent/'PhysiCelldFBA/output')]):
    s = read(folder/'output00000017.xml')
    pos = s['labels']['position']
    dead = scalar(s, 'dead').astype(bool)
    growth = scalar(s, 'growth_rate')
    radius = np.linalg.norm(pos[:2], axis=0)
    visible = abs(pos[2]) <= scalar(s, 'radius')
    colors = plt.cm.Greens(np.clip(growth/.832260843339191, 0, 1)*.8+.1)
    colors[dead] = [56/255,38/255,0,1]
    axes[0,col].scatter(*pos[:2,visible], c=colors[visible], s=1, rasterized=True)
    axes[0,col].set(title=f'{label}: z=0 section, 510 min', xlabel='x (µm)', ylabel='y (µm)', aspect='equal', xlim=(-120,120), ylim=(-120,120))
    bins = np.arange(0,151,5)
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (radius>=lo)&(radius<hi)
        if mask.any():
            rows.append(dict(project=label,radius_um=(lo+hi)/2,cells=int(mask.sum()),dead_fraction=float(dead[mask].mean()),mean_stored_growth=float(growth[mask].mean())))
    selected = [r for r in rows if r['project']==label]
    axes[1,0].plot([r['radius_um'] for r in selected], [r['dead_fraction'] for r in selected], label=label)
    axes[1,1].plot([r['radius_um'] for r in selected], [r['mean_stored_growth'] for r in selected], label=label)
axes[1,0].set(xlabel='Distance from colony origin (µm)',ylabel='Dead fraction (all z layers)')
axes[1,1].set(xlabel='Distance from colony origin (µm)',ylabel='Mean stored growth rate (1/h)')
for ax in axes[1]: ax.legend(); ax.grid(alpha=.2)
fig.savefig(HERE/'death_patterns.png', dpi=170)
with (HERE/'radial_death.csv').open('w') as f:
    w=csv.DictWriter(f, fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps(sorted([r for r in rows if r['project']=='PhysiCelldFBA'],key=lambda r:r['dead_fraction'],reverse=True)[:6],indent=2))
