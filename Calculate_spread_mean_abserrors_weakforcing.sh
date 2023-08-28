#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --time=11:59:00
#SBATCH --output=stat_calc.sh.%j.out
#SBATCH --error=stat_calc.sh.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=K.Tempest@physik.uni-muenchen.de 

# load modules
module load python

# run the model
python Calculate_spread_mean_abserror.py
