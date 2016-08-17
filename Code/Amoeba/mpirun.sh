#!/bin/csh -f
if ($#argv != 6) then
  echo "Usage: threadrun minthread maxthread timedir exec inputs machfile"
  exit 1
endif

set minProcs = $1
set maxProcs = $2

set timedir  = $3
set exec     = $4
set inputs   = $5
set machfile = $6

set procs = $minProcs
echo "executable name = $exec"
echo "input file name = $inputs"
echo "time directory  = $timedir"

if (! -e $timedir) then
    mkdir $timedir
endif

while ($procs <= $maxProcs)

 echo "number of procs = $procs"
 echo "setenv OMP_NUM_THREADS 1"
 setenv OMP_NUM_THREADS 1
 echo "setenv CH_OUTPUT_INTERVAL 10000"
 setenv CH_OUTPUT_INTERVAL 10000

 echo "machfile = $machfile"
 echo "exec = $exec"
 echo "mpirun -np $procs -machinefile $machfile $exec $inputs"
 mpirun -np $procs -machinefile $machfile $exec $inputs

 set pprocs = `echo $procs | awk '{printf("%06d",$1);}'`    

 echo "mv time.table.0 $timedir/amrg.mpi.time.proc0.$pprocs.thread"
 mv time.table.0 $timedir/amrg.mpi.time.proc0.$pprocs.thread

 set procs = `expr 1 +  $procs`

end
