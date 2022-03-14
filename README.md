# Orphans

Study Gamma-ray Bursts orphan afterglows with Rubin LSST and FINK


## Installing from source
This is the recommended installation procedure in a conda environment.
It shall take care of the `afterglowpy` dependency cleanly.

```
> git clone git@gitlab.in2p3.fr:johan-bregeon/orphans.git
or
> https://gitlab.in2p3.fr/johan-bregeon/orphans.git
then
> cd orphans
> conda env create -f environment.yml
> conda activate orphans
> pip install -r requirements.txt
> pip install -e .
```