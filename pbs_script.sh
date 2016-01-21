# This is just a sample PBS script. You should adapt the
# following example script, so that you can run all your experiments.
# For more information about Torque and PBS, check the following sources:
#    https://support.cc.gatech.edu/facilities/instructional-labs/how-to-run-jobs-on-the-instructional-hpc-clusters
#    http://www.democritos.it/activities/IT-MC/documentation/newinterface/pages/runningcodes.html


# allocate 4 of the (24 total) sixcore nodes for up to 5 minutes

#PBS -q class
#PBS -l nodes=4:sixcore
#PBS -l walltime=00:05:00
#PBS -N cse6220-nqueens

# TODO: change this to your project directory relative to your home directory
#       (= $HOME)
export PBS_O_WORKDIR=$HOME/path/to/your/project


# TODO:
# for your experiments, run your program for different `n`, different master-depth `k`
# and with different number of processors. This will allow you to properly
# analyze the scalability behaviour of each parameter setting.
# Also: you should time each parameter setting multiple times and then average
#       the results to achieve a smoother and more accurate analysis.
MASTER_DEPTH=4
EXE=$PBS_O_WORKDIR/nqueens

# loop over number of processors (just an example, uncomment to use the loop)
# for p in 2 4 8 12 16 24
# do
p=24 # (four sixcore processors => 24 processes total)
N=8 # TODO: loop over N as well!
OMPI_MCA_mpi_yield_when_idle=0 mpirun --hostfile $PBS_NODEFILE -np $p $EXE -t $N $MASTER_DEPTH
# done

