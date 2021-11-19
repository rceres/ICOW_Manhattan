  #PBS -l nodes=1:ppn=20
	#PBS -l walltime=13:00:00
	#PBS -j oe
	#PBS -o output_test.txt
  #PBS -A cef13_b_g_sc_default
	echo "set workdir"
	cd $PBS_O_WORKDIR
	echo "load openmpi"
  module load gcc/5.3.1 openmpi/1.10.1

  
  echo " "
  echo " "
  echo "Job started on `hostname` at `date` , `time`"
##  echo "Running on nodes: " 
##  cat $PBS_NODEFILE
  
  echo         "../BorgParICOW.exe ~/scratch/CFK_050000/CFK_050000_output.txt ~/scratch/CFK_050000/ 0.935 0.0 0.260 0.0 0.232 0.0 12 1000 .05 .04 .005 0.02 .0002 0.01 .02 .02 .02"
  mpirun -np 20 ../BorgParICOW.exe ~/scratch/CFK_050000/CFK_050000_output.txt ~/scratch/CFK_050000/ 0.935 0.0 0.260 0.0 0.232 0.0 12 1000 .05 .04 .005 0.02 .0002 0.01 .02 .02 .02
  echo " "
  echo "Job Ended at `date` `time`"
  echo " "