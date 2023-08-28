#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=90G
#SBATCH --time=23:59:00
#SBATCH --output=run_0.sh.%j.out
#SBATCH --error=run_0.sh.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=K.Tempest@physik.uni-muenchen.de 

# load modules
module load python

# run the model

python run_script.py --DA=1

python run_script.py --start_time=0 --batch_number=0 --members_per_batch=1000 --DA=0
python run_script.py --start_time=0 --batch_number=1 --members_per_batch=1000 --DA=0
python run_script.py --start_time=0 --batch_number=2 --members_per_batch=1000 --DA=0
python run_script.py --start_time=0 --batch_number=3 --members_per_batch=1000 --DA=0
python run_script.py --start_time=0 --batch_number=4 --members_per_batch=1000 --DA=0

python run_script.py --start_time=60 --batch_number=0 --members_per_batch=1000 --DA=0
python run_script.py --start_time=60 --batch_number=1 --members_per_batch=1000 --DA=0
python run_script.py --start_time=60 --batch_number=2 --members_per_batch=1000 --DA=0 
python run_script.py --start_time=60 --batch_number=3 --members_per_batch=1000 --DA=0 
python run_script.py --start_time=60 --batch_number=4 --members_per_batch=1000 --DA=0 

python run_script.py --start_time=120 --batch_number=0 --members_per_batch=1000 --DA=0
python run_script.py --start_time=120 --batch_number=1 --members_per_batch=1000 --DA=0
python run_script.py --start_time=120 --batch_number=2 --members_per_batch=1000 --DA=0
python run_script.py --start_time=120 --batch_number=3 --members_per_batch=1000 --DA=0
python run_script.py --start_time=120 --batch_number=4 --members_per_batch=1000 --DA=0

python run_script.py --start_time=180 --batch_number=0 --members_per_batch=1000 --DA=0
python run_script.py --start_time=180 --batch_number=1 --members_per_batch=1000 --DA=0
python run_script.py --start_time=180 --batch_number=2 --members_per_batch=1000 --DA=0
python run_script.py --start_time=180 --batch_number=3 --members_per_batch=1000 --DA=0
python run_script.py --start_time=180 --batch_number=4 --members_per_batch=1000 --DA=0

python run_script.py --start_time=240 --batch_number=0 --members_per_batch=1000 --DA=0
python run_script.py --start_time=240 --batch_number=1 --members_per_batch=1000 --DA=0
python run_script.py --start_time=240 --batch_number=2 --members_per_batch=1000 --DA=0
python run_script.py --start_time=240 --batch_number=3 --members_per_batch=1000 --DA=0
python run_script.py --start_time=240 --batch_number=4 --members_per_batch=1000 --DA=0

python run_script.py --start_time=300 --batch_number=0 --members_per_batch=1000 --DA=0
python run_script.py --start_time=300 --batch_number=1 --members_per_batch=1000 --DA=0
python run_script.py --start_time=300 --batch_number=2 --members_per_batch=1000 --DA=0
python run_script.py --start_time=300 --batch_number=3 --members_per_batch=1000 --DA=0
python run_script.py --start_time=300 --batch_number=4 --members_per_batch=1000 --DA=0
