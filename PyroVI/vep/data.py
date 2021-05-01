# trivial data loading for VEP data

import os
import bz2
import pickle
import zipfile
import glob
import logging
import numpy as np
import multiprocessing
import torch
import pyro

logger = logging.getLogger(__name__)
here = os.path.abspath(os.path.dirname(__file__))

_param_store_key = "__pyro_param_store_state__"


def save_zip(zfname, locals_):
    pyro.get_param_store().save(zfname + ".ps")
    with zipfile.ZipFile(zfname, "w") as zf:
        # logger.debug('saving pyro param store')
        # with zf.open(_param_store_key, 'w') as fd:
        #     torch.save(pyro.get_param_store().get_state(), fd)
        # logger.debug('saving vep objects')
        for key, val in locals_.items():
            try:
                with zf.open(key, "w") as fd:
                    pickle.dump(val, fd)
                logger.debug("dumped %s to %s:%s", key, zfname, key)
            except (AttributeError, TypeError) as exc:
                logger.debug("failed to dump %s %s: %s", key, val, exc)


def load_zip(zfname):
    objs = {}
    pyro.get_param_store().load(zfname + ".ps")
    with zipfile.ZipFile(zfname, "r") as zf:
        # logger.debug('restoring state of pyro param store')
        # with zf.open(_param_store_key, 'r') as fd:
        #     pyro.get_param_store().set_state(pickle.load(fd))
        # logger.debug('loading vep objects')
        for zi in zf.filelist:
            zi: zipfile.ZipInfo
            try:
                with zf.open(zi.filename, "r") as fd:
                    objs[zi.filename] = pickle.load(fd)
            except EOFError as exc:
                logger.debug(
                    "load_zip %s saw %s for entry %s", zfname, exc, zi.filename
                )
    return objs


def rload(fname):
    """Load a dict of data from an R dump format file.
    """
    if isinstance(fname, str):
        if fname.endswith(".bz2"):
            fd = bz2.open(fname, "rt")
        elif fname.endswith(".R"):
            fd = open(fname, "r")
        else:
            raise ValueError("unknown file type for {}".format(fname))
    else:
        fd = fname
    with fd:
        lines = fd.readlines()
    data = {}
    for line in lines:
        lhs, rhs = [_.strip() for _ in line.split("<-")]
        if rhs.startswith("structure"):
            *_, vals, dim = rhs.replace("(", " ").replace(")", " ").split("c")
            vals = [float(v) for v in vals.split(",")[:-1]]
            dim = [int(v) for v in dim.split(",")]
            val = np.array(vals).reshape(dim[::-1]).T
        elif rhs.startswith("c"):
            val = np.array([float(_) for _ in rhs[2:-1].split(",")])
        else:
            try:
                val = int(rhs)
            except:
                try:
                    val = float(rhs)
                except:
                    raise ValueError(rhs)
        data[lhs] = val
    return data


def load_data_dir_id(id, **kwargs):
    parts = here, "..", "data", "id*%03d*" % id, "vep", "data.R.bz2"
    pattern = os.path.join(*parts)
    logger.debug("load_data_dir_id pattern %r", pattern)
    matches = glob.glob(pattern)
    if len(matches) != 1:
        msg = "{} datasets found for patient {}".format(len(matches), id)
        raise IOError(msg)
    logger.debug("load_data_dir_id match %r", matches)
    return Patient(matches[0], **kwargs)


def _load_many_data_dir_id_helper(args):
    id, kwargs = args
    try:
        result = load_data_dir_id(id, **kwargs)
    except IOError as e:
        logger.warn(e)
        result = None
    return id, result


def load_many_data_dir_id(ids, **kwargs):
    dev = kwargs.pop("dev", None)
    argss = [(id, kwargs) for id in ids]
    with multiprocessing.Pool() as p:
        results = p.map(_load_many_data_dir_id_helper, argss)
    results = {id: result.to(dev) for id, result in results if result}
    return results


def load_gsw(data):
    "Load gain matrix, seeg data and weights matrix."
    if isinstance(data, str):
        data = rload(data)

    Sn = data["seeg"]
    Sn[~np.isfinite(Sn)] = 0.0

    dtype = torch.FloatTensor
    G = torch.tensor(data["gain"]).type(dtype)
    S = torch.tensor(Sn).type(dtype)  # [:, np.r_[0:43, 107:153]]
    Wut = data["counts_triu"]
    Wut /= Wut.max()

    ns, nn = G.shape
    ns, nt = S.shape
    assert G.shape[0] == S.shape[0]
    W = torch.zeros(nn, nn)
    p = 0
    for i in range(nn):
        for j in range(nn):
            if i > j:
                W[i, j] = W[j, i] = Wut[p]
                p += 1
    return G, S, W


# TODO this is getting big enough that we should have a test harness in place
# but there is no time for that right now


class PatientData:
    def __init__(self, gain, sensors, weights, dev=None):
        self.gain = gain
        self.sensors = sensors
        self.weights = weights
        self.dev = dev
        self.to_dev()

    def to_dev(self):
        if self.dev is not None:
            for key in "gain sensors weights".split(" "):
                setattr(self, key, getattr(self, key).to(self.dev))
        return self


class Patient(PatientData):
    def __init__(self, rfile, start=0, end=-1, dev=None):
        self.rfile = rfile
        gain, sensors, weights = load_gsw(rfile)
        sensors = sensors[:, start:end]
        super().__init__(gain, sensors, weights, dev=dev)


class ResampledPatient:
    def __init__(self, patient, rs_sensors, rs_windows, rs_winsize):
        self.patient = patient
        self.rs_sensors = rs_sensors
        self.rs_windows = rs_windows
        self.rs_winsize = rs_winsize

    def __iter__(self):
        return self

    def __next__(self):
        ix_sensors = np.random.randint(0, self.patient.n_sensor, self.rs_sensors)
        ix_win_t0 = np.random.randint(0, self.patient.sensors.shape[1], self.rs_windows)
        ix_win_ts = np.c_[ix_win_t0] + np.r_[: self.rs_winsize]
        ix_win_ts %= self.patient.sensors.shape[1]
        gain = self.patient.gain[ix_sensors]
        sensors = self.patient.sensors[ix_sensors][:, ix_win_ts]
        weights = self.patient.weights
        return PatientData(gain, sensors, weights)
