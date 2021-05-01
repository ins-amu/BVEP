import os
import argparse
from subprocess import check_output, STDOUT

pulled_images = {
    'pytorch': 'nvcr.io/nvidia/pytorch:18.06-py3'
}


def build_slurm_opts(
        partition='debug', constraint='gpu', account='ich001',
        time='00:15:00'):
    return [f'--{key}={value}'
            for key, value in locals().items()]


def build_slurm_cmd(*cmd, slurm_cmd='srun', **kwargs):
    return [slurm_cmd] + build_slurm_opts(**kwargs) + list(cmd)


def build_shifter_cmd(*cmd, image_key='pytorch'):
    shifter = ['shifter', 'run']
    home = os.environ['HOME']
    shifter.append(f'--mount=source={home},destination={home},type=bind')
    shifter.append(pulled_images[image_key])
    return shifter + list(cmd)


def build_parser():
    from argparse import ArgumentParser, REMAINDER
    parser = ArgumentParser()
    parser.add_argument('-p', '--partition', default='debug')
    parser.add_argument('-t', '--time', default='00:15:00')
    parser.add_argument('-s', '--sweep', nargs='+', action='append')
    parser.add_argument('-f', '--format', nargs=2, action='append', default=[])
    parser.add_argument('-l', '--local', nargs='?', const=-1, type=int)
    parser.add_argument('rest', nargs=REMAINDER)
    return parser


def build_one_sweep_sequence(key_value):
    if len(key_value) == 1:
        return []
    key, *values = key_value
    dashes = '-' if len(key) == 1 else '--'
    for val in values:
        yield f'{dashes}{key} {val}'.split(' ')


def build_sweep_sequence(sweep):
    sequences = map(build_one_sweep_sequence, sweep)
    from itertools import product
    for args in product(*sequences):
        zero = []
        for arg in args:
            zero += arg
        yield zero


def main():
    parser = build_parser()
    args = parser.parse_args()
    print(args)
    rest = args.rest[:]
    # use -- as sentinel value to split jobs args from vep args
    if rest and rest[0] == '--':
        rest = rest[1:]
    sweep = args.sweep or []
    cmds = []
    for i, sarg in enumerate(build_sweep_sequence(sweep)):
        fmt_args = []
        sarg_slug = ' '.join(sarg).replace('-', '').replace(' ', '-')
        for key, fmt in args.format:
            fmt_args.append(f'{"-" if len(key) == 1 else "--"}{key}')
            fmt_args.append(fmt.format(i=sarg_slug))
        vep_cmd = 'python -m vep'.split(' ') + rest + list(sarg) + fmt_args

        log_fname = f'log-{sarg_slug}.txt'
        if args.local:
            if args.local > 0:
                if not ('-C' in vep_cmd or '--cpus' in vep_cmd):
                    vep_cmd.append('--cpus')
                    vep_cmd.append('1')
            print(' '.join(vep_cmd))
            cmds.append((vep_cmd, log_fname))
        else:
            shifter_cmd = build_shifter_cmd(*vep_cmd)
            slurm_cmd = build_slurm_cmd(
                *shifter_cmd,
                partition=args.partition,
                time=args.time,
            )
            print(' '.join(slurm_cmd))
            cmds.append((slurm_cmd, log_fname))

    if args.local:
        if args.local < 0:
            [run_and_log(vep_cmd, log_fname) for vep_cmd, log_fname in cmds]
        else:

            import multiprocessing
            with multiprocessing.Pool(args.local) as p:
                p.map(_run_and_log_mp_f, cmds)

    else:
        threads = []
        for slurm_cmd, log_fname in cmds:
            threads.append(async_check_call(slurm_cmd, log_fname))
        [_.join() for _ in threads]


def _run_and_log_mp_f(args):
    return run_and_log(*args)

def run_and_log(cmd, out_fname):
    with open(out_fname, 'wb') as fd:
        fd.write(check_output(cmd, stderr=STDOUT))

def async_check_call(cmd, log_fname):
    from threading import Thread
    thread = Thread(target=run_and_log, args=(cmd, log_fname))
    thread.start()
    return thread


if __name__ == '__main__':
    main()
