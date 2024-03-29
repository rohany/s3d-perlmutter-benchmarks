
# DESTINATION=./window-changes-scaling/man-trace/ TRACE_CONFIG="manual" sbatch -N 1 ./ammonia_job.sh
# DESTINATION=./window-changes-scaling/man-trace/ TRACE_CONFIG="manual" sbatch -N 2 ./ammonia_job.sh
# DESTINATION=./window-changes-scaling/man-trace/ TRACE_CONFIG="manual" sbatch -N 4 ./ammonia_job.sh
# DESTINATION=./window-changes-scaling/man-trace/ TRACE_CONFIG="manual" sbatch -N 8 ./ammonia_job.sh
DESTINATION=./window-changes-scaling/man-trace/ TRACE_CONFIG="manual" sbatch -N 16 ./ammonia_job.sh
# DESTINATION=./window-changes-scaling/auto-trace/ TRACE_CONFIG="auto" sbatch -N 1 ./ammonia_job.sh
# DESTINATION=./window-changes-scaling/auto-trace/ TRACE_CONFIG="auto" sbatch -N 2 ./ammonia_job.sh
# DESTINATION=./window-changes-scaling/auto-trace/ TRACE_CONFIG="auto" sbatch -N 4 ./ammonia_job.sh
# DESTINATION=./window-changes-scaling/auto-trace/ TRACE_CONFIG="auto" sbatch -N 8 ./ammonia_job.sh
DESTINATION=./window-changes-scaling/auto-trace/ TRACE_CONFIG="auto" sbatch -N 16 ./ammonia_job.sh
# DESTINATION=./window-changes-scaling/no-trace/ TRACE_CONFIG="no" sbatch -N 1 ./ammonia_job.sh
# DESTINATION=./window-changes-scaling/no-trace/ TRACE_CONFIG="no" sbatch -N 2 ./ammonia_job.sh
# DESTINATION=./window-changes-scaling/no-trace/ TRACE_CONFIG="no" sbatch -N 4 ./ammonia_job.sh
# DESTINATION=./window-changes-scaling/no-trace/ TRACE_CONFIG="no" sbatch -N 8 ./ammonia_job.sh
DESTINATION=./window-changes-scaling/no-trace/ TRACE_CONFIG="no" sbatch -N 16 ./ammonia_job.sh


