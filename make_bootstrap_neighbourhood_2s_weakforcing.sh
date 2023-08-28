#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=23:59:00
#SBATCH --output=neigh_bootstrap_2.sh.%j.out
#SBATCH --error=neigh_bootstrap_2.sh.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=K.Tempest@physik.uni-muenchen.de 

# load modules
module load python

# run the model
python make_bootstrap_neighbourhood_2.py --grid_point=500 --times=24 --times_diff=15 --times_start=14 --boots=10000 --n_half_size=0
python make_bootstrap_neighbourhood_2.py --grid_point=500 --times=24 --times_diff=15 --times_start=14 --boots=10000 --n_half_size=1
python make_bootstrap_neighbourhood_2.py --grid_point=500 --times=24 --times_diff=15 --times_start=14 --boots=10000 --n_half_size=2
python make_bootstrap_neighbourhood_2.py --grid_point=500 --times=24 --times_diff=15 --times_start=14 --boots=10000 --n_half_size=3
python make_bootstrap_neighbourhood_2.py --grid_point=500 --times=24 --times_diff=15 --times_start=14 --boots=10000 --n_half_size=4
python make_bootstrap_neighbourhood_2.py --grid_point=500 --times=24 --times_diff=15 --times_start=14 --boots=10000 --n_half_size=5
python make_bootstrap_neighbourhood_2.py --grid_point=500 --times=24 --times_diff=15 --times_start=14 --boots=10000 --n_half_size=6
