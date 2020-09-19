#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -n 1
source ~/.bashrc
source ~/py3/bin/activate
export PYTHONPATH=$PYTHONPATH:/home/samgeen/Programming/pymses_python3
python stellarextinctions.py
