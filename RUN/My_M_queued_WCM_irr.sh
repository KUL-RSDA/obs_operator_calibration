#!/usr/bin/bash

#SBATCH --job-name=My_M_queued_WCM
#SBATCH --time=59:59:00
#SBATCH --nodes=1 --ntasks-per-node=31
#SBATCH --account=lp_ees_swm_ls_001
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=gabrielle.delannoy@kuleuven.be
#SBATCH --clusters=genius
#SBATCH -o script_log_irr.txt
#SBATCH -e script_out_irr.txt

cd $SLURM_SUBMIT_DIR

NCPU=`wc -l < $SLURM_JOB_NODELIST`
echo ------------------------------------------------------
echo 'This job is allocated on '$NCPU' cpu(s)'
echo 'Job is running on node(s): '
cat $SLURM_JOB_NODELIST
echo ------------------------------------------------------

module load cluster/genius/centos7
module load matlab
#module load matlab

path_script=/data/leuven/314/vsc31402/GIT_REPO/WCM_CALI_main/
script=CalWCM_main_PAR_CPU.m
exp_option=po_ol_hymap_irr
tag=q_irr
N_cpu=31

#==========================================

in_script="$path_script$script"
I_file="./I_$tag$script"
IO_file="./IO_$tag$script"
IOO_file="./tmp_$tag$script"

./replace.pl $in_script $I_file EXP $exp_option

./replace.pl $I_file $IO_file PAR_N_CPU $N_cpu

rm -f $I_file

if [ $N_cpu -gt 0 ]; then
   k=1

   while [ $k -le $N_cpu ]
   do
          echo "preparing processor $k out of $N_cpu..."

          cd /data/leuven/314/vsc31402/GIT_REPO/WCM_CALI_main/RUN/

          ./replace.pl $IO_file $IOO_file PAR_CPU $k

          log_file="/scratch/leuven/314/vsc31402/4DMED/OUT_LOG/out.cali_$tag$k.$$"

          matlab  -nojvm -nosplash -singleCompThread -display null < $IOO_file > $log_file &

          sleep 1

          echo "CALIBRATION: track in $log_file"

          rm -f ./$IOO_file

          k=$(($k+1))
   done
fi

rm -f $IO_file

echo Done!

wait

#==========================================================
# Launch the job with 
# >> sbatch My_M_queued_WCM_irr.sh
#==========================================================
