import subprocess

for i in range(1000):
	sub.process('sbatch script_orphans.sh', shell=True)
