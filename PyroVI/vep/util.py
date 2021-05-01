import os
import glob
import numpy as np


def _load_one(npz):
    npz = np.load(npz)
    # seeg time series, gain matrix, weights, EZ hypothesis, onset delay
    return [npz[_] for _ in "seeg_t G w ezh od ezs".split(" ")]


def load_data(sid=None):
    npzs = {}
    for npz in glob.glob("data/*.npz"):
        npzs[os.path.basename(npz).split(".npz")[0]] = npz
    if sid is not None:
        return _load_one(npzs[sid])
    return {sid: _load_one(npz) for sid, npz in npzs.items()}


def roi_names():
    with open("data/destrieux.txt", "r") as fd:
        lines = fd.readlines()
    return np.array([_.strip() for _ in lines])


def prec_recl(ezh, vep, vt=0.1):
    """Compute precision and recall, for "truth" `ezh`, score `vep` and threshold `vt`"""
    tp = ezh.astype(bool) & (vep > vt)
    fp = ~ezh.astype(bool) & (vep > vt)
    fn = ezh.astype(bool) & ~(vep > vt)
    prec = tp.sum() / (tp.sum() + fp.sum())
    recl = tp.sum() / (tp.sum() + fn.sum())
    return prec, recl


def prec_recl_auto_pct(ezh, vep):
    pr = []
    # TODO biset O(log n) instead of sweep O(n)
    for pct in np.r_[90:100:30j]:
        pr.append(prec_recl(ezh, vep, np.percentile(vep, pct)) + (pct,))
    return max(pr, key=lambda t: t[0] + (1 if t[1] > 0.5 else 0))


def prec_recl_auto(ezh, vep):
    p, r, pct = prec_recl_auto_pct(ezh, vep)
    return p, r


def best_pr_over_sz(score, ezh):
    """Compute best P, R over seizure scores `score` and "true" ez hypothesis `ezh`"""
    pr = []
    for score_ in score:
        fpr, tpr, thresholds = sklearn.metrics.roc_curve(ezh, score_)
        # XXX thresholds[-2] is not magical, should rewrite
        pr.append(prec_recl(ezh, score_, thresholds[-2]))
    return pr[np.argmax([p for p, r in pr])]


def make_prs_per_sz(scores, ezhs_):
    """For scores and "true" EZ of all patients `scores` and `ezhs_`, compute PRs"""
    p, r = np.array([best_pr_over_sz(*_) for _ in zip(scores, ezhs_)]).T
    return p[np.isfinite(p)].mean(), p, r


def create_region_bias():
    data = load_data()
    # construct per region bias as mean of all ez hypotheses
    es = np.array([e for _, _, _, e, _ in data.values()]).mean(axis=0)
    es += es[es > 0].min()
    np.save("data/region_bias.npy", es)


def region_bias():
    return np.load("data/region_bias.npy")


if __name__ == "__main__":
    create_region_bias()
    print(region_bias())
