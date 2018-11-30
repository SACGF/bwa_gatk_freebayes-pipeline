#!/bin/bash

### Job name
#PBS -N console

### Join queuing system output and error files into a single output file
#PBS -j oe

### Send email to user when job ends or aborts
#PBS -m ae

### email address for user
#PBS -M jinghua.feng@adelaide.edu.au

### Request nodes, memory
#PBS -l nodes=1:ppn=1
#PBS -l mem=500mb,vmem=1gb
#PBS -l walltime=500:00:00

# this is the one testing with qualimap; NOTE the ralative path used for the '-o' parameter of qsub
sacgf_py3_env_runner.sh snakemake -s ./test.snake --printshellcmds -j 10 --jobname "{params.job_name}.{jobid}" --cluster "qsub -o ~/../../../scratch/jfeng/test_snake/logs/ -j oe {params.email} -l nodes=1:ppn={params.t} -l mem={params.mem}mb,vmem={params.vmem}mb  -l walltime={params.walltime}"



sacgf_py3_env_runner.sh snakemake -s /home/users/jfeng/bioinformatics/scripts/pipelines/mapping_and_variant_calling/rules/bwa_gatk_workflow.2.snake \
    -j 4 --printshellcmds --jobname "{params.job_name}.{jobid}" \
    --cluster-config /home/users/jfeng/bioinformatics/scripts/pipelines/mapping_and_variant_calling/templates/torque.json \
    --cluster "qsub -j oe -l {cluster.ppn} -l {cluster.mem} -l {cluster.walltime}"

