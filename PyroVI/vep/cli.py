import argparse
import logging
import zipfile
import pickle
import torch
import pyro
import pyro.optim

from .data import load_many_data_dir_id, save_zip
from .advi import Advi, MapGuide, NonFiniteLoss
from .model import ExpGain, Coupling, ARStep, Model, Cohort


def build_parser():
    parser = argparse.ArgumentParser(description="Run VEP model")

    # debug options
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Print debug information"
    )
    parser.add_argument(
        "-V", "--validate", action="store_true", help="Enable Pyro validation"
    )

    # algo options
    parser.add_argument(
        "-d", "--device", type=str, default="cpu", help="Device to use (cpu or cuda)"
    )
    parser.add_argument(
        "-C",
        "--cpus",
        type=int,
        default=1,
        help="Number of OpenMP threads to use " "(default: all available)",
    )
    parser.add_argument("-l", "--mll_tol", default=-2.0, type=float)
    parser.add_argument("-m", "--max_steps", default=0, type=int)
    parser.add_argument("-n", "--n_steps", default=100, type=int)
    parser.add_argument(
        "-o",
        "--map_est",
        default=False,
        action="store_true",
        help="Compute a MAP estimate (instead of SVI)",
    )
    parser.add_argument(
        "-T",
        "--max_time",
        default="0.0",
        type=str,
        help="Max walltime in seconds, minutes or hours" "(23m, 5.5h)",
    )
    parser.add_argument(
        "-r", "--rate", default=0.1, type=float, help="Learning rate for optimizer"
    )
    parser.add_argument(
        "--optim", default="Adam", type=str, help="Name of optimizer to use"
    )
    parser.add_argument(
        "-N",
        "--num_particles",
        default=1,
        type=int,
        help="Number of gradient evaluations in MC integration" "of the ELBO",
    )
    parser.add_argument(
        "--restarts",
        nargs=1,
        default=10,
        type=int,
        help="Number of times to restart due to instability.",
    )

    # input data options
    parser.add_argument(
        "-P", "--new-cohort", action="store_true", help="Run new cohort data & model"
    )
    parser.add_argument(
        "-p",
        "--patient_id",
        nargs="+",
        default=1,
        type=int,
        help="Patient id from data/ directory",
    )
    parser.add_argument("-s", "--start", default=0, type=int)
    parser.add_argument("-e", "--end", default=-1, type=int)

    # output options
    parser.add_argument(
        "-f", "--figure", type=str, help="Filename for saving QC figure"
    )
    parser.add_argument(
        "-z", "--zip", type=str, help="Create ZIP archive of pickled objects"
    )
    parser.add_argument(
        "--psis", action="store_true", help="Compute & print PSIS LOO values"
    )

    return parser


def main(args):
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)
    logger = logging.getLogger(__name__)
    logger.debug("args %s", args)
    pyro.clear_param_store()

    if args.cpus and args.cpus > 0:
        torch.set_num_threads(args.cpus)

    if args.validate:
        logger.info("enabling Pyro validation " "(performance penalties may ensue)")
        pyro.enable_validation()

    dev = torch.device(args.device)
    logger.debug("have device %r", dev)

    logger.info("loading data for pids %s", args.patient_id)
    patients = []
    all_data = load_many_data_dir_id(
        args.patient_id, dev=dev, start=args.start, end=args.end
    )
    for pid, data in all_data.items():
        logger.debug("pid %d sensors size %s", pid, data.sensors.size())

        observer = ExpGain(data.gain, prefix="p{}_forward".format(pid))
        coupling = Coupling(data.weights, prefix="p{}_coupling".format(pid))
        arstep = ARStep(
            data.n_node,
            data.n_time,
            coupling=coupling,
            data_like=data.gain,
            prefix="p{}_sources".format(pid),
        )
        model = Model(arstep, observer, sensors=data.sensors)
        patients.append(model)

    cohort = Cohort(patients)

    # parse time
    if args.max_time.lower().endswith("m"):
        max_time_f = float(args.max_time[:-1]) * 60.0
    elif args.max_time.lower().endswith("h"):
        max_time_f = float(args.max_time[:-1]) * 3600.0
    else:
        max_time_f = float(args.max_time)
    logger.debug("max_time %r parsed to %.3f seconds", args.max_time, max_time_f)

    # ideally a restarter would be a clever wrapper, but it's easier to do it here
    # checkpoints would be the way to go
    restarts_left = args.restarts
    while restarts_left > 0:
        try:
            Optim = getattr(pyro.optim, args.optim)
            optim = Optim({"lr": args.rate})
            svi = Advi(
                cohort,
                opt=optim,
                num_particles=args.num_particles,
                guide=MapGuide(cohort) if args.map_est else None,
                data_like=patients[-1].observer.gain,
                n_steps=args.n_steps,
            )
            svi.step_until(
                patients,
                max_steps=args.max_steps,
                mll_tol=args.mll_tol,
                psis=args.psis,
                max_time=max_time_f,
            )
            break
        except NonFiniteLoss as oops:
            msg = "%s encountered, %d restarts left"
            msg %= oops, restarts_left
            logger.warning(msg)
            restarts_left -= 1
            pyro.clear_param_store()
    else:
        msg = "%d restarts failed to remain stable! quitting early"
        msg %= (args.restarts,)
        logger.error(msg)
        return {}

    cohort.report_params(svi.guide)

    if args.figure:
        logger.info("plotting results, saving to %s", args.figure)
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.pyplot import style

        style.use("dark_background")
        with PdfPages(args.figure) as pdf:
            for pid, patient in zip(args.patient_id, patients):
                patient.qc_plot(svi.guide)
                pdf.savefig()

    if args.zip:
        # TODO need to persist params as well
        logger.info("dumping objects to %s", args.zip)
        save_zip(args.zip, locals())

    return locals()


def do(line):
    return main(build_parser().parse_args(line.split(" ")))
