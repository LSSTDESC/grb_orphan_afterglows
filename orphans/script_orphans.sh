#!/bin/sh


# SLURM options

#SBATCH --job-name=orphans		# Nom du job
#SBATCH --output=orphans_simu_%j.log    # Standard output et error log

#SBATCH --partition=htc                 # Choix de partition (htc par défaut)

#SBATCH --ntasks=1                      # Exécuter une seule tâche
#SBATCH --mem=2000                      # Mémoire en MB par défaut
#SBATCH --time=0-00:15:00               # Délai max = 7 jours


# Commandes à soumettre :

ccenv anaconda
conda activate /pbs/home/m/mmasson/lsst/orphans_env
	
export RUBIN_SIM_DATA=/pbs/home/m/mmasson/lsst/rubin_sim_data
export RUBIN_SIM=/pbs/home/m/mmasson/lsst/rubin_sim
export DUSTMAPS=/pbs/home/m/mmasson/lsst/orphans/data/schlafly_dust_factor.csv
export SIMU=/pbs/home/m/mmasson/lsst/orphans/data/simulations
export OBS=/pbs/home/m/mmasson/lsst/orphans/data/pseudo_obs
	
#module load python
python generate_grb_pop.py