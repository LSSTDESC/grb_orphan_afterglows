# Orphans

Study Gamma-ray Bursts orphan afterglows with Rubin LSST and FINK

## Authors
- Marina Masson <marina.masson@lpsc.in2p3.fr>
- Johan Bregeon <bregeon@in2p3.fr>

## Installing from source
This is the recommended installation procedure in a conda environment.
It shall take care of the `afterglowpy` dependency cleanly.

```
> git clone git@github.com:LSSTDESC/grb_orphan_afterglows.git
or
> https://github.com/LSSTDESC/grb_orphan_afterglows.git
then
> cd orphans
> conda env create -f environment.yml
> conda activate orphans
> pip install -r requirements.txt
> conda install markupsafe==2.0.1
> pip install -e .
```

then if needed, one shall also install the `rubin_sim` package as follows:
```
# from the orphans directory
> cd ..
> git clone https://github.com/lsst/rubin_sim.git
> cd rubin_sim
> conda install -c conda-forge --file=requirements.txt
> pip install -e .
> cd ../orphans
```
in which case, one shall follow instructions to download the scheduler data, and then set:
```
> export RUBIN_SIM_DATA_DIR=$MYWORKDIR/rubin_sim/rubin_sim_data"
```
