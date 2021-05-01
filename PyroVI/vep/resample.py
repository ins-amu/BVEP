from pylab import *
import bz2, glob
import math
import os
import torch
import torch.distributions.constraints as constraints
import pyro
from pyro.optim import Adam
from pyro.infer import SVI, Trace_ELBO
import pyro.distributions as dist
from pyro.contrib.autoguide import AutoDiagonalNormal
import torch, torch.optim, torch.autograd
from torch.autograd import Variable

from .data import rload

dtype = torch.FloatTensor


# small override of ADVI guide to ensure the arrays it creates
# match dtype & device of our data, e.g. float32 on cuda
class AdviDev(AutoDiagonalNormal):
    def sample_latent(self, *args, **kwargs):
        zeros, ones = args[0].new_zeros, args[0].new_ones
        loc = pyro.param("{}_loc".format(self.prefix), lambda: zeros(self.latent_dim))
        scale = pyro.param(
            "{}_scale".format(self.prefix),
            lambda: ones(self.latent_dim),
            constraint=constraints.positive,
        )
        return pyro.sample(
            "_{}_latent".format(self.prefix),
            dist.Normal(loc, scale).independent(1),
            infer={"is_auxiliary": True},
        )


nc, nw, ns = 256, 64, 32  # 1 MB / patient
npz_fname = f"resamp-{nc}-{nw}-{ns}.npz"
try:
    npz = np.load(npz_fname)
    print("loaded cached data")
    S, G, W = [npz[key] for key in "SGW"]
except IOError:
    print("reading data")
    b, *_ = bzs = glob.glob("data/*/*/*.bz2")
    for bn in bzs:
        with bz2.BZ2File(bn) as b:
            with open(bn.split(".bz2")[0], "wb") as fd:
                fd.write(b.read())
    Rs = []
    for rn in glob.glob("data/*/*/data.R"):
        Rs.append(rload(rn))

    print("resampling")

    def wut(G, Wut):
        _, nn = G.shape
        W = torch.zeros(nn, nn)
        p = 0
        for i in range(nn):
            for j in range(nn):
                if i > j:
                    W[i, j] = W[j, i] = Wut[p]
                    p += 1
        W /= W.max()
        return W

    S, G, W = [], [], []
    for R in Rs:
        ic = np.random.permutation(R["ns"])[np.r_[:nc] % R["ns"]]
        it = np.random.permutation(R["nt"])[np.r_[:nw] % R["nt"]]
        it = (np.c_[it] + np.r_[:ns]) % R["nt"]
        s = R["seeg"]
        s = s[ic]
        s = s[:, it]
        S.append(s)
        G.append(R["gain"][ic])
        W.append(wut(R["gain"], R["counts_triu"]).numpy())
    S, G, W = np.array(S), np.array(G), np.array(W)
    np.savez(npz_fname, S=S, G=G, W=W)
S = S.transpose((3, 0, 1, 2)).copy()
print(S.shape, G.shape, W.shape)
print("dataset is %d KB" % (sum([_.nbytes for _ in (S, G, W)]) >> 10,))

S[~np.isfinite(S)] = 0.0
assert all([np.isfinite(_).all() for _ in (S, G, W)])

Gt = torch.tensor(G).type(dtype)
St = torch.tensor(S).type(dtype)
Wt = torch.tensor(W).type(dtype)
dev = torch.device(os.environ.get("dev", "cpu"))
pyro.clear_param_store()
Gt, St, Wt = [_.to(dev) for _ in (Gt, St, Wt)]
print(Gt.size(), St.size())  # , Wt.size()

print("starting model")

# ideally, this can become another cohort model
# TODO rewrite resampling as part of data loading
# TODO write new cohort class which applies resampling
# TODO adapt other classes to new shapes or write new ones
# TODO add flags to make use of cohort
# TODO TODO TODO put in place eavluation against hypotheses
def model(Gt, St, Wt):
    T, Z, O, N = Gt.new_tensor, Gt.new_zeros, Gt.new_ones, dist.Normal
    (nt, np, nc, nw), (_, _, nn) = St.size(), Gt.size()
    # src bounded -2 to 14
    Src = pyro.sample(
        "src", dist.Uniform(Z(nt, np, nn, nw) - T(2.0), T(14.0)).independent(4)
    )
    # src[1:] ~ normal(src[:-1]*0.9 + w.dot(src[:-1])*0.5), x0)
    mu = 0.9 * Src + 0.5 * Wt.matmul(Src)
    x0 = pyro.sample("x0", N(Z(np, nn, 1), T(0.2)))
    pyro.sample("step", N(mu[:-1], x0.exp()), obs=Src[1:])
    Mu = Gt.matmul(Src)
    pyro.sample("S", N(Mu, T(0.1)).independent(4), obs=St)


# question is how to interpret the result: if we assume a linear system
# in source space, the stability of the each node as predicted by an LSA
# should be reflected in the amplitude of the source space.
# Two difficulties (1) the prior from dataset & (2) regularize
# W/o regularization, we should be able to match results of LSA alone

# TODO generate new new data w/ non-temporal averaging
# TODO adapt to new datasets w/ hypothesis, train to predict hypotheses
# TODO this is the ICEI project, if it works, it's a big deal

import time

guide = AdviDev(model)
svi = SVI(model, guide, Adam({"lr": 0.001}), loss=Trace_ELBO())
trace = []
print("first step")
last_loss = svi.step(Gt, St, Wt)
print("starting iters")
tic = time.time()
n_iter = 100
for step in range(n_iter):
    loss = svi.step(Gt, St, Wt)
    print(step, loss, abs(loss - last_loss))
    last_loss = loss
print(guide.median(Gt, St, Wt)["x0"])
toc = time.time()
et = toc - tic
per_iter = et / n_iter
print("elapsed time {0} ({1} s / iter)".format(et, per_iter))
