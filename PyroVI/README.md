## PyroVI

*This is a historical, experimental implementation of VEP models
in Pyro*

While the data used here can't be included, it should be straightforward
to read through the `vep.cli` module to see how the model components defined
in `vep.model` are built and used for variational inference.

### Usage

To make use in batch systems easier, a command line interface is provided, e.g. run patient 23, until time point 40, converge with change in log lik < 10.0, save figure to qc.png, objects to obj.pickle and use a GPU,
```bash
python -m vep -p 23 -e 40 -l 10.0 -f qc.png -P obj.pickle -d cuda
```
but this is also accessible from Python as `vep.cli.do` returning all objects in a dictionary which can be used as follows:


```python
import vep.cli
globals().update(vep.cli.do('-p 23 -e 40 -l 1.0'))
```

    INFO:vep.advi:step   100 elbo     25251.26 delbo 667.00 ll    -24317.36 dll  23.27
    INFO:vep.advi:step   200 elbo     24667.64 delbo 424.49 ll    -22025.49 dll  29.43
    INFO:vep.advi:step   300 elbo     25082.92 delbo 638.43 ll    -19707.26 dll  31.05
    INFO:vep.advi:step   400 elbo     23468.36 delbo 461.01 ll    -17256.87 dll  19.35
    INFO:vep.advi:step   500 elbo     23422.73 delbo 188.64 ll    -15313.34 dll  23.92
    INFO:vep.advi:step   600 elbo     23230.78 delbo 215.40 ll    -13668.80 dll   7.63
    INFO:vep.advi:step   700 elbo     22811.75 delbo 103.30 ll    -12573.07 dll   3.98
    INFO:vep.advi:step   800 elbo     22236.07 delbo   0.81 ll    -11902.70 dll   1.50
    INFO:vep.advi:step   900 elbo     22111.61 delbo 166.43 ll    -11492.06 dll   4.37
    INFO:vep.advi:step  1000 elbo     21845.88 delbo  16.06 ll    -11170.81 dll   3.67
    INFO:vep.advi:step  1100 elbo     21801.57 delbo 152.22 ll    -10997.58 dll   2.71
    INFO:vep.advi:step  1200 elbo     21768.13 delbo  91.87 ll    -10894.36 dll   2.36
    INFO:vep.advi:step  1300 elbo     21390.81 delbo 118.17 ll    -10780.66 dll   0.09
    INFO:vep.advi:elapsed time 17.279 s; 13.3 ms/step
    INFO:vep.advi:psis-loo -11192.246 max-k 1.106 90%-k 0.247



```python
%pylab inline
model.qc_plot(svi.guide)
```

    Populating the interactive namespace from numpy and matplotlib



![png](README_files/README_2_1.png)

*can't include figure, data is not anonymized*


- plain AR model clearly strong, use uniform + centered to get same effect
- letting it run for longer provides good improvement
- interesting that we can minibatch sensors, nodes or time points separately (in principle)
- could use posterior predictive mean + 6D Epileptor to show what model understood
- per node variance has lower LOO
- pow scaling W doesn't help

Cohort level modeling would potentially help fill in info on under/over-estimated connections per DWI, in terms of improving explanation of pathology.

Full help for CLI (ignore part about ipykernel_launcher.py)


```python
vep.cli.do('--help')
```
    
    (pyro) [duke@dukebox VEP]$ python -m vep --help
    usage: __main__.py [-h] [-p PATIENT_ID [PATIENT_ID ...]] [-v] [-V] [-d DEVICE]
                       [-s START] [-e END] [-l MLL_TOL] [-m MAX_STEPS]
                       [-n N_STEPS] [-f FIGURE] [-P PICKLE] [--psis] [-T MAX_TIME]
    
    Run VEP model
    
    optional arguments:
      -h, --help            show this help message and exit
      -p PATIENT_ID [PATIENT_ID ...], --patient_id PATIENT_ID [PATIENT_ID ...]
                            Patient id from data/ directory
      -v, --verbose         Print debug information
      -V, --validate        Enable Pyro validation
      -d DEVICE, --device DEVICE
                            Device to use (cpu or cuda)
      -s START, --start START
      -e END, --end END
      -l MLL_TOL, --mll_tol MLL_TOL
      -m MAX_STEPS, --max_steps MAX_STEPS
      -n N_STEPS, --n_steps N_STEPS
      -f FIGURE, --figure FIGURE
                            Filename for saving QC figure
      -P PICKLE, --pickle PICKLE
                            Filename for saving pickle of objects
      --psis                Compute & print PSIS LOO values
      -T MAX_TIME, --max_time MAX_TIME
                            Max walltime in seconds, minutes or hours(23m, 5.5h)



