#!/bin/bash
#SBATCH -J bmntd_R       		# job name
#SBATCH -o log_slurm.o%j    # output and error file name (%j expands to jobID)
#SBATCH -n 48 			    # total number of tasks requested
#SBATCH -N 1 			    # number of nodes you want to run on
#SBATCH -p bsudfq			# queue (partition)
#SBATCH -t 150:00:00 		# run time (hh:mm:ss)
#SBATCH --mail-user=jessicabernardin@boisestate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Load the module
module load borah-base/default
module load borah-misc r/4.2.2

#Replace myscript.R with the name of the script
Rscript bmntd.R
