#!/bin/bash
#

# (node-specific) constraint options from /etc/slurm/slurm.conf:
# ib,avx,avx512,lnvo6g

working_dir=${PWD}
job_name=$1

output_dir=${PWD}
input_dir=${PWD}

cat > .sbatch_script.${job_name}.$$ << EOF
#!/bin/bash
#

#SBATCH --job-name=${job_name}
#SBATCH --output=${output_dir}/${job_name}.sc
#SBATCH --error=${output_dir}/${job_name}.sc
#SBATCH --constraint="ib"
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --partition=core

cd ${working_dir}

date 

python3 NextToNewestGeneration.py $1 $1 
echo "I finished the basis optimization."
date

python3 A3_lit_M.py $1 $1 
echo "I calculated the coupling matrices."
date

#\rm sbatch_script.${job_name}.$$

EOF

sbatch .sbatch_script.${job_name}.$$
rm .sbatch_script.${job_name}.$$