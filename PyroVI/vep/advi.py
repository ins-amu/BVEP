# provides ADVI from Pyro pieces

import time
import logging
import numpy as np
import torch
import torch.distributions.constraints as constraints
import pyro
import pyro.distributions as dist
from pyro.contrib.autoguide import biject_to, AutoDiagonalNormal, AutoDelta
from pyro.optim import Adam
from pyro.infer import SVI, Trace_ELBO, JitTrace_ELBO

from .psis import psisloo

logger = logging.getLogger(__name__)


class NonFiniteLoss(ValueError):
    pass


class GuideExtrasMixin:
    """Mixes in with a guide to provide a method to
    """

    def sample(self):
        """Sample from the guide's latent distribution, unpacking the
        parameters similarly to the :meth:`median` method.
        """
        sample = self.sample_latent().detach()
        result = {}
        return {
            site["name"]: biject_to(site["fn"].support)(unconstrained_value)
            for site, unconstrained_value in self._unpack_latent(sample)
        }


class MeanFieldGuide(AutoDiagonalNormal, GuideExtrasMixin):
    """MeanFieldGuide implements a mean-field variational distribution based
    on `pyro.contrib.autoguide.AutoDiagonalNormal`, with some extras
    """

    def __init__(self, *args, data_like=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.data_like = data_like

    def sample_latent(self, *args, **kwargs):
        if self.data_like is None:
            # there has to be a better way..
            for arg in args:
                if hasattr(arg, "new_zeros"):
                    self.data_like = arg
                elif hasattr(arg, "data_like") and arg.data_like is not None:
                    self.data_like = arg.data_like
                else:
                    self.data_like = torch.tensor(0.0)
        zeros, ones = self.data_like.new_zeros, self.data_like.new_ones
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


class MapGuide(AutoDelta, GuideExtrasMixin):
    pass


class Advi:
    def __init__(
        self,
        model,
        guide=None,
        elbo=None,
        opt=None,
        lr=0.002,
        num_particles=1,
        n_steps=10,
        data_like=None,
        jit=False,
    ):
        self.model = model
        self.guide = guide or MeanFieldGuide(model, data_like=data_like)
        Elbo = JitTrace_ELBO if jit else Trace_ELBO
        self.elbo = elbo or Elbo(num_particles=num_particles)
        self.opt = opt or Adam({"lr": lr})
        self.svi = SVI(self.model, self.guide, self.opt, loss=self.elbo)
        self.steps_taken = 0
        self.n_steps = n_steps
        self.loss = []
        self.ll = []  # log likelihood

    def dloss(self):
        if len(self.loss) > 2:
            return np.abs(self.loss[-1] - self.loss[-2])
        return np.inf

    def dll(self):
        if len(self.ll) > 2:
            return np.abs(np.diff(np.array(self.ll[-self.n_steps :]))).mean()
        return np.inf

    def _check_loss(self, loss, *last_loss):
        if np.isfinite(loss):
            return loss
        msg = "Non-finite loss (%f), last value %s"
        msg %= loss, last_loss
        raise NonFiniteLoss(msg)

    def step(self, *args, n_steps=0, log=True):
        for i in range(n_steps or self.n_steps):
            loss = self.svi.step(*args)
            self._check_loss(loss, self.loss[-1:])
            self.steps_taken += 1
            self.loss.append(loss)
            if hasattr(self.model, "log_lik"):
                self.ll.append(self.model.log_lik(self.guide))
            else:
                self.ll.append(loss)
        if log and self.steps_taken > 2:
            logger.info(
                "step %5d elbo %12.2f delbo %6.2f mll %6.3f",
                self.steps_taken,
                self.loss[-1],
                self.dloss(),
                self.ll[-1],
            )

    def step_until(
        self, *args, max_steps=0, delbo_tol=0.0, mll_tol=None, max_time=0, psis=True
    ):
        if not (
            max_steps > 0 or delbo_tol > 0.0 or mll_tol is not None or max_time > 0
        ):
            logger.warn(
                "no termination condition provided to step_until,"
                " please use Ctrl-C to stop."
            )
        tic = time.time()
        stop = False
        # TODO add condition for norm of elbo gradient instead of delbo
        while True:
            try:
                if max_steps > 0 and self.steps_taken > max_steps:
                    cond = "max steps exceeded"
                    stop = True
                if delbo_tol > 0.0 and self.loss and self.dloss() < delbo_tol:
                    cond = "delbo tolerance exceeded"
                    stop = True
                if mll_tol is not None and self.ll and self.ll[-1] > mll_tol:
                    cond = "mll tolerance achieved"
                    stop = True
                if max_time > 0 and (time.time() - tic) > max_time:
                    cond = "max time exceeded"
                    stop = True
                if stop:
                    logger.info("step_until stopping on condition %s", cond)
                    break
                self.step(*args, log=True)
            except KeyboardInterrupt:
                logger.info("step_until stopping on condition keyboard " "interrupt")
                break
        toc = time.time()
        et = toc - tic
        logger.info(
            "elapsed time %.3f s; %.1f ms/step", et, 1000 * et / self.steps_taken
        )
        if psis:
            self._compute_psis()

    def _compute_psis(self):
        if isinstance(self.guide, MapGuide):
            msg = "skipping requested PSIS for MAP estimate"
            logger.warning(msg)
            return
        logger.debug("computing PSIS LOO")
        self.psis = self.compute_psis()
        logger.info(
            "psis-loo %.3f max-k %.3f 90%%-k %.3f",
            self.psis["loo"],
            self.psis["ks"].max(),
            np.percentile(self.psis["ks"], 90),
        )
        return self.psis

    def checkpoint(self):
        # save stuff to disk in case we run out of time!
        # http://pyro.ai/examples/dmm.html#Checkpointing
        raise NotImplementedError

    def compute_psis(self, n_sample=100):
        lls = []
        for i in range(n_sample):
            ll = self.model.log_lik(self.guide, mean=False, mode="sample")
            ll = ll.detach().cpu().numpy().flat[:]
            lls.append(ll)
        lls = np.array(lls)
        return psisloo(lls)
