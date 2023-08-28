#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=00:59:00
#SBATCH --output=convergence_plot.sh.%j.out
#SBATCH --error=convergence_plot.sh.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=K.Tempest@physik.uni-muenchen.de 

# load modules
module load python

# run the script
python make_convergence_eescomparison_plot.py --time=11 # 16 hours
