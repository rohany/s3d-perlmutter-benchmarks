trace_configs=("manual" "auto" "no")
nodes=(1 2 4 8 16)
grids=("64x64x64" "128x128x64" "128x128x64")

for g in ${grids[@]}; do
    for c in ${trace_configs[@]}; do
        for n in ${nodes[@]}; do
            GRID="$g" DESTINATION="./sweeptest/$g/$c-trace/" TRACE_CONFIG="$c" sbatch -N $n ./ammonia_job.sh
	done
    done
done
