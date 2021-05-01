import logging
import numpy as np
import torch
import pyro
from pyro.distributions import Uniform, Normal


logger = logging.getLogger(__name__)


class Component:
    def __init__(self, prefix=None, data_like=None):
        self.data_like = data_like
        self.prefix = prefix or self.__class__.__name__.lower()


class Coupling(Component):

    def __init__(self, weights, lag=1, **kwargs):
        super().__init__(**kwargs)
        self.weights = weights
        self.lag = lag

    def __call__(self, sources):
        coupling = self.weights.mm(sources[:, :-self.lag])
        return coupling


class DiffCoupling(Coupling):

    def __call__(self, sources):
        n_node, n_time = sources.size()
        from_state = sources[:, :-self.lag].view(n_node,      1, n_time - 1)
        to_state   = sources[:, :-self.lag].view(     1, n_node, n_time - 1)
        diff_state = from_state - to_state
        weighted_diff = self.weights.view(n_node, n_node, 1) * diff_state
        coupling = weighted_diff.sum(0)
        return coupling


class Gain(Component):

    def __init__(self, gain, amp=1.0, offset=0.0, **kwargs):
        kwargs['data_like'] = gain
        super().__init__(**kwargs)
        self.gain = gain
        self.data_like = gain
        self.amp = amp
        self.offset = offset

    def propagate(self, sources):
        return self.gain.mm(sources)

    def sensor_mean(self, sources):
        return self.amp * self.propagate(sources) + self.offset

    def dist(self, sources, eps):
        sensor_mean = self.sensor_mean(sources)
        sensor_std = eps.exp() * eps.new_tensor(0.7)
        return Normal(sensor_mean, sensor_std)

    def __call__(self, sources, sensors, eps=0.0):
        sample = pyro.sample(
            "{}_obs".format(self.prefix),
            self.dist(sources, eps=eps).independent(2),
            obs=sensors
        )
        return sample


class ExpGain(Gain):

    def propagate(self, sources):
        return self.gain.mm(sources.exp()).log()


class ARStep(Component):

    def __init__(self, n_node, n_time, coupling=None, **kwargs):
        super().__init__(**kwargs)
        self.n_node = n_node
        self.n_time = n_time
        self.coupling = coupling

    @property
    def sources_key(self):
        return '{}_sources'.format(self.prefix)

    def sources_dist(self):
        lower = self.data_like.new_zeros(self.n_node, self.n_time) - 2.0
        upper = lower + 14.0
        return Uniform(lower, upper).independent(2)

    def sources_sample(self):
        sample = pyro.sample(self.sources_key, self.sources_dist())
        return sample

    def step_dist(self, sources, coupling, ar, k, sig):
        k = k.exp() * k.new_tensor(0.1)
        mean = sources[:, :-1] * ar + k * coupling(sources)
        std = self.data_like.new_tensor(sig.exp() * 1.4)
        return Normal(mean, std).independent(2)

    def step_sample(self, sources, coupling, ar, k, sig, **kwargs):
        return pyro.sample(
            "{}_step".format(self.prefix),
            self.step_dist(sources, coupling, ar, k, sig),
            **kwargs
        )

    def __call__(self, coupling=None, ar=0.9, k=0.0, sig=0.0):
        sources = self.sources_sample()
        coupling = coupling or self.coupling
        self.step_sample(sources, coupling, ar, k, sig, obs=sources[:, 1:])
        return sources


class Model:

    def __init__(self, source_est, observer, sensors):
        self.source_est = source_est
        self.sensors = sensors
        self.observer = observer

    def log_lik(self, guide, mean=True, mode='median'):
        sources_mean = getattr(guide, mode)()[self.source_est.sources_key]
        # TODO not good
        eps_mean = getattr(guide, mode)()['cohort_eps']
        lp = self.observer.dist(sources_mean, eps_mean).log_prob(self.sensors)
        if mean:
            return lp.detach().mean().item()
        else:
            return lp

    def __call__(self, ar=0.9, k=0.0, eps=0.0, sig=0.0):
        sources = self.source_est(ar=ar, k=k, sig=sig)
        sample = self.observer(sources=sources,
                               sensors=self.sensors,
                               eps=eps)
        return sample

    def qc_plot(self, guide, fig=(15, 5)):
        import pylab as pl
        if fig:
            pl.figure(figsize=fig)
        to_np = lambda a: a.detach().cpu().numpy()
        eps_mean = guide.median()['cohort_eps']
        sources = guide.median()[self.source_est.sources_key]
        sensors_mu = self.observer.dist(sources, eps_mean).mean
        sources, sensors_mu = [to_np(_) for _ in (sources, sensors_mu)]
        ax1 = pl.subplot(141)
        pl.imshow(self.sensors, aspect='auto', interpolation='none')
        pl.title('Sensors')
        ax4 = pl.subplot(144, sharex=ax1)
        pl.imshow(sources, aspect='auto', interpolation='none')
        pl.title('Sources')
        ax3 = pl.subplot(143, sharex=ax1)
        pl.plot(sources.T, 'w', alpha=0.1)
        pl.title('Sources')
        pl.grid(1)
        ax2 = pl.subplot(142, sharex=ax1)
        pl.imshow(sensors_mu, aspect='auto', interpolation='none')
        pl.title('Sensor prediction')


class Cohort(Component):

    def __init__(self, patient_models):
        super().__init__(
            prefix="cohort", data_like=patient_models[0].observer.gain)
        self.patient_models = patient_models

    def cohort_params(self):
        t = self.data_like.new_tensor
        N = lambda mu, sd: Normal(t(mu), t(sd))
        # TODO would be useful to provide guide with prior
        ar = pyro.sample("{}_ar".format(self.prefix), N(0.9, 0.2))
        k = pyro.sample("{}_k".format(self.prefix), N(0.0, 0.2))
        eps = pyro.sample("{}_eps".format(self.prefix), N(0.0, 0.2))
        sig = pyro.sample("{}_sig".format(self.prefix), N(0.0, 0.2))
        return ar, k, eps, sig

    def __call__(self, *args, **kwargs):
        # TODO propagating scalars like this is a pain. modules?
        ar, k, eps, sig = self.cohort_params()
        for patient in self.patient_models:
            patient(ar=ar, k=k, eps=eps, sig=sig)

    def log_lik(self, guide, mean=True, mode='median'):
        lls = []
        for patient in self.patient_models:
            lls.append(patient.log_lik(guide, mean, mode))
            if not mean:
                lls[-1] = lls[-1].view(-1)
        return np.mean(lls) if mean else torch.cat(lls)

    def report_params(self, guide):
        qf = getattr(guide, 'quantiles', None)
        if qf:
            self.report_quantiles(guide)
        else:
            self.report_map(guide)

    def report_quantiles(self, guide, qs=frozenset([0.05, 0.5, 0.95])):
        logger.info('reporting quantiles %s for cohort', qs)
        for key, val in guide.quantiles(sorted(qs)).items():
            msg = '%s' + '\t%+06.3f' * len(val)
            if key.startswith(self.prefix):
                logger.info(msg, key, *val)

    def report_map(self, guide):
        for key, val in guide.median().items():
            msg = '%s MAP estimate ' + '\t%+06.3f'
            if key.startswith(self.prefix):
                logger.info(msg, key, val)