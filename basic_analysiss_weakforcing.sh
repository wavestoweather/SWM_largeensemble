#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=220G
#SBATCH --time=03:00:00
#SBATCH --output=basic_analysis.sh.%j.out
#SBATCH --error=basic_analysis.sh.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=K.Tempest@physik.uni-muenchen.de 

# load modules
module load python

# run the script
python basic_analysis.py 
