#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --time=11:59:00
#SBATCH --output=make_neighbourhood_distribution_data.sh.%j.out
#SBATCH --error=make_neighbourhood_distribution_data.sh.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=K.Tempest@physik.uni-muenchen.de 

# load modules
module load python

# run the model
python make_neighbourhood_distribution_data.py --grid_point=500
