#!/bin/bash -eu
#SBATCH -A m4411
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 30:00
#SBATCH --exclusive

#INFO: you need to adjust the paths in the call below
module load PrgEnv-gnu
module load cudatoolkit
module load cray-pmi
#module load python

DEV=/pscratch/sd/r/rohany/legion_s3d/
MECHANISM=hept
KERNEL_PATH="$DEV/kernels/${MECHANISM}_kernels"
export LD_LIBRARY_PATH=$DEV/build/$MECHANISM:$KERNEL_PATH:$DEV/legion/language/build/lib/:$LD_LIBRARY_PATH

GRID=${GRID:-64x64x64}

# export FI_CXI_DEFAULT_CQ_SIZE=13107200
# export FI_CXI_REQ_BUF_MIN_POSTED=10
# export FI_CXI_REQ_BUF_SIZE=25165824
# export FI_MR_CACHE_MONITOR=memhooks
# export FI_CXI_RX_MATCH_MODE=software
export GASNET_OFI_DEVICE_0=cxi0
export GASNET_OFI_DEVICE_1=cxi1
export GASNET_OFI_DEVICE_2=cxi2
export GASNET_OFI_DEVICE_3=cxi3
export GASNET_OFI_DEVICE_TYPE=Node
export GASNET_OFI_RECEIVE_BUFF_SIZE=2M
#export GASNET_OFI_RECEIVE_BUFF_SIZE=single
#export MPICH_MAX_THREAD_SAFETY=multiple
#export MPICH_OFI_NIC_POLICY=NUMA

num_nodes=$SLURM_JOB_NUM_NODES
num_ranks=$(( 4 * SLURM_JOB_NUM_NODES ))

# ranks_x=( [1]=1 [2]=2 [6]=4 [12]=4 [24]=4 [48]=8 [96]=8 [192]=8  [384]=16  [768]=16  [1536]=16)
# ranks_y=( [1]=2 [2]=2 [6]=3 [12]=3 [24]=6 [48]=6 [96]=6 [192]=12 [384]=12  [768]=12  [1536]=24)
# ranks_z=( [1]=2 [2]=2 [6]=2 [12]=4 [24]=4 [48]=4 [96]=8 [192]=8  [384]=8   [768]=16  [1536]=16)

ranks_x=( [1]=2 [2]=2 [4]=4 [8]=4 [16]=4 [32]=8 [64]=8 [128]=8 [256]=16 [512]=16 [1024]=16)
ranks_y=( [1]=2 [2]=2 [4]=2 [8]=4 [16]=4 [32]=4 [64]=8 [128]=8 [256]=8  [512]=16 [1024]=16)
ranks_z=( [1]=1 [2]=2 [4]=2 [8]=2 [16]=4 [32]=4 [64]=4 [128]=8 [256]=8  [512]=8  [1024]=16)


nx=${ranks_x[$num_nodes]}
ny=${ranks_y[$num_nodes]}
nz=${ranks_z[$num_nodes]}

OUTDIR=$DESTINATION/./pwave_x_"$num_nodes"_"$MECHANISM"
mkdir -p $OUTDIR

if [ "$TRACE_CONFIG" = "manual" ]; then
  RHST_TRACE=1
  LEGION_AUTO_TRACE_ARGS=""
  LEGION_TRACE_ARGS="-lg:window 8192 -lg:sched 4096 -lg:no_transitive_reduction"
elif [ "$TRACE_CONFIG" = "auto" ] ; then
  RHST_TRACE=0
  LEGION_AUTO_TRACE_ARGS="-lg:enable_automatic_tracing -lg:auto_trace:batchsize 5000 -level auto_trace=2 -lg:auto_trace:identifier_algorithm multi-scale -lg:auto_trace:multi_scale_factor 500 -lg:auto_trace:repeats_algorithm quick_matching_of_substrings -lg:auto_trace:min_trace_length 25"
  LEGION_TRACE_ARGS="-lg:window 8192 -lg:sched 4096 -lg:no_transitive_reduction"
else
  RHST_TRACE=0
  LEGION_AUTO_TRACE_ARGS=""
  LEGION_TRACE_ARGS=""
fi

./ammonia_s3d.py  \
    -x $DEV/fortran/$MECHANISM/s3d.x -l $DEV/build/$MECHANISM/librhst.so \
    -e RHSF_MAPMODE=allgpu -e GASNET_FREEZE_ON_ERROR=0 -e RHSF_FUSE=200 -e RHST_TRACE=$RHST_TRACE -g $GRID -v 1x1x1 -t 150 \
    -b $DEV/s3d/ -p "$nx"x"$ny"x"$nz" -s slurm -e GASNET_BACKTRACE=1 \
    -e RHST_ARGS="-ll:gpu 1 -ll:cpu 1 -ll:util 3 -ll:bgwork 3 -ll:fsize 32768 -ll:csize 22528 -ll:rsize 1024 -ll:gsize 0 -ll:zsize 640 -ll:stacksize 16 -logfile run_%.log -cuda:lmemresize 1 -cuda:hostreg 0 -lg:prof $num_ranks -lg:prof_logfile prof_%.gz $LEGION_AUTO_TRACE_ARGS $LEGION_TRACE_ARGS" \
    --bomb-file=0d -k $OUTDIR -o $OUTDIR/ammonia.txt --title pressure_wave_test
