import sys
from pylab import *
from parse_csv import parse_csv
import common
import matplotlib.pyplot as plt

npz = np.load('data.R.npz')
Ic = npz['Ic']
ps = npz['seeg_log_power']
W_ = npz['SC']
nodes = npz['nodes']
gain = npz['gain']

data = parse_csv(sys.argv[1])
opt = data['x0'].shape[0]==1
n = data['x0'].shape[1]



# connectivity
#im = lambda i: imshow(i, interpolation='none', aspect='auto', cmap='binary')
#figure(); subplot(121); im(W_); subplot(122); im(data['FC'].mean(axis=0));
#savefig('figure-x2.png', dpi=300)


# figure X
figure(figsize=(8, 5))
showpicks = r_[0]
amp, offset = [data[_].reshape((-1, 1, 1)) for _ in 'amp offset'.split()]
ps_est =  amp*(log(exp(data['x']).dot(gain.T)) + offset)
pelo, pehi = percentile(ps_est, [5, 95], axis=0)
t = r_[:pelo.shape[0]]

np.save('pehi.npy', pehi)
np.save('pelo.npy', pelo)
np.save('t.npy', t)

subplot(211)
for i, (lo, hi) in enumerate(zip(pelo.T, pehi.T)):
    if i in showpicks:
        fill_between(t, lo, hi, alpha=0.5, facecolor='b', edgecolor='none')
plot(ps[:, showpicks], 'k', alpha=0.5)
leg = []
#[9, 46, 60, 130, 138, 139]
for i in [15]:
    cn = common.bipnames[i]
    leg.append(cn)
legend(leg)
xlabel('Time (s)'); ylabel('sEEG(t)'); 
subplot(212)
violinplot(data['x0'], widths=0.7);
plt.axhline(y=-2.05, linewidth=.5, color = 'r', linestyle='--')
text(0.2, -2, "$x_{0c}$", size=10, color = 'r')
plt.axhline(y=-1.65, linewidth=.5, color = 'g', linestyle='--')
#text(0.2, -1.6, "$\Delta x$", size=10, color = 'g')

import pandas as pd


node_name=[]
node_region=[]
node_abbrev = []
for i in nodes:
    rn = common.reg_names[i]
    node_region.append(rn)
    ra = ''.join([p[0] for p in rn.split('-')[1:] if p]).upper()
    node_abbrev.append(ra)
    print (ra, ' '.join(rn.split('-')))
    node_name.append(ra)
    

    #df = pd.DataFrame(rn, ra)
    #df.to_csv('x.csv')
    
#import pandas as pd
#df = pd.DataFrame( node_region, node_abbrev)
#df.to_csv('nodes.csv')
      
xticks(r_[1:len(nodes)+1], node_abbrev, rotation=30)
ylabel("$P(x_0)$")
ylim([-4.5, -0.5])
tight_layout()
text(-2, 5, "A", size=14)
text(-2, 0, "B", size=14)
savefig('figure-x.png', dpi=300)

with open('report.txt', 'w') as fd:
    fd.write("index,x0mean,x0median,x0std,abbrev,name\n")
    fmt = "%d,%f,%f,%f,%s,%s\n"
    for i, idx in enumerate(nodes):
        i_x0 = data['x0'][:, i]
        args = idx, np.mean(i_x0), np.median(i_x0), np.std(i_x0), node_abbrev[i], common.reg_names[idx]
        fd.write(fmt % args)


with open('reportx0.txt', 'w') as fd:
    fd.write("#index,x0mean,x0std\n")
    fmt = "%d,%f,%f\n"
    for i, idx in enumerate(nodes):
        i_x0 = data['x0'][:, i]
        args = idx, np.mean(i_x0), np.std(i_x0)
        fd.write(fmt % args)

with open('reportnodes.txt', 'w') as fd:
    fd.write("#index,abbrev,name\n")
    fmt = "%d,%s,%s\n"
    for i, idx in enumerate(nodes):
        i_x0 = data['x0'][:, i]
        args = idx, node_abbrev[i], common.reg_names[idx]
        fd.write(fmt % args)

show()
"""
[common.reg_names[i] for i in nodes]
get_ipython().magic('pinfo violinplot')
xticks(r_[1:len(nodes)+1], [common.reg_names[i] for i in nodes], rotation=45)
cla()
violinplot(data['x0'], vert=True)
violinplot(data['x0'], vert=False)
cla()
violinplot(data['x0'], vert=False)
yticks(r_[1:len(nodes)+1], [common.reg_names[i] for i in nodes])
tight_layout()
tight_layout()
subplot(211);
tight_layout()
yticks(r_[1:len(nodes)+1], [[p[0] for p in common.reg_names[i].split()[1:]] for i in nodes])
yticks(r_[1:len(nodes)+1], [[p[0] for p in common.reg_names[i].split('-')[1:]] for i in nodes])
yticks(r_[1:len(nodes)+1], [''.join([p[0] for p in common.reg_names[i].split('-')[1:]]) for i in nodes])
"""
