#!/bin/bash
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=p_mc1689_1
#SBATCH --job-name=Hsub-S31
#SBATCH --output=slurm.Hsub-S31.out
#SBATCH --error=slurm.Hsub-S31.err
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END
#SBATCH --mail-user=brian.keane@rutgers.edu
#Run the MATLAB command
cd /home/keanebp/projects/NeuralMech/docs/scripts/
/home/keanebp/projects/NeuralMech/docs/scripts//Step1_hcp_preprocessing_msmall_amarel.sh --server='amarel' --preFS='false' --FS='false' --postFS='false' --fmriVol='true' --fmriSurf='false' --restFix='false' --msmAll='false' --dedriftResample='false' --subj='sub-S31'
