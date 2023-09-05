#!/bin/bash

echo $(date)

set -x

exe=../NOARCH.x

# export CUDA_VISIBLE_DEVICES=2,3
export QUDA_RESOURCE_PATH="/gpfs/alpine/phy131/proj-shared/djmurphy/software/src/cps.quda_oddPC_eofa_merge/cps_pp/tests/eofa_force/run/work"
export OMP_NUM_THREADS=24

now=$(date '+%Y%m%d%H%M%S')
log_file=log-${now}.out
mpirun -n 4 ${exe} ./vmls/do_arg.vml ./vmls/doext_arg.vml ./vmls/evo_arg.vml ./vmls/quda_arg.vml 1.0e-06 ../vmls/eofa_arg.vml -qmp-geom 1 1 1 4 | tee logs/${log_file}

echo $(date)
