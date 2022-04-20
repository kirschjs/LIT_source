#!/bin/bash
#

working_dir=${PWD}
job_name=$1

output_dir=/tmp/compton_output
input_dir=${PWD}

cat > .sbatch_script.${job_name}.$$ << EOF
#!/bin/bash
#

#SBATCH --job-name=${job_name}
#SBATCH --output=${output_dir}/${job_name}.sc
#SBATCH --error=${output_dir}/${job_name}.sc
#SBATCH --partition=core
#SBATCH --constraint="avx"
#SBATCH --ntasks=1

cd ${working_dir}

date 

mpiexec python3 NextToNewestGeneration.py 

echo "done!"
date
#\rm sbatch_script.${job_name}.$$

EOF

sbatch .sbatch_script.${job_name}.$$
#rm .sbatch_script.${job_name}.$$