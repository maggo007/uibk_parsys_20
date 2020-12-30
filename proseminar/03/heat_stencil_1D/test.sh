#!/bin/bash

# Execute job in the queue "std.q" unless you have special requirements.
#$ -q std.q

# The batch system should use the current directory as working directory.
#$ -cwd

# Name your job. Unless you use the -o and -e options, output will
# go to a unique file name.ojob_id for each job.
#$ -N kopp_heat_1D

##$ -M markus.kopp@student.uibk.ac.at
##$ -m e

# Join the error stream to the output stream.
#$ -j yes

M=8

#module load openmpi/4.0.3
for X in 2 4 8 16 32 64; do
    if (("$X" < "$M"));
    then
        echo "#$ -pe openmpi-"$X"perhost $X"
        #$ -pe openmpi-$Xperhost $X
    else
        echo "#$ -pe openmpi-8perhost $X"
        #$ -pe openmpi-8perhost $X
    fi
    echo "new test with N=$X"
    time mpiexec -n $X --oversubscribe --report-bindings --display-devel-map ./heat_stencil_1D_mpi 1024
    echo "--------------------"
done

