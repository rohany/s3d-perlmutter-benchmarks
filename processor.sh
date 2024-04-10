trace_configs=("manual" "auto" "no")
grids=("64x64x64" "128x64x64" "128x128x64")

for g in ${grids[@]}; do
    for c in ${trace_configs[@]}; do
        echo "GRID=$g TRACE=$c"
	python iteration.py "./sweeptest/$g/$c-trace/"
    done
done
