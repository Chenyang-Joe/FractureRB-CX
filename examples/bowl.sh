#!/bin/bash

# Record start time
start=$(date +%s.%N)

# Move to script directory (optional safety)
cd "$(dirname "$0")"

# Run the simulation
gdb --args ../build/FractureRB bowl.bullet bowl.csv \
  -o _out/bowl_ \
  -n 500 \
  -i 1e2 \
  -f 5e5 \
  -s 2 \
  --default-solver

# Record end time
end=$(date +%s.%N)

# Calculate duration
duration=$(echo "$end - $start" | bc)

# Format duration
printf "\nStart Time: %s\nEnd Time: %s\nDuration : %.2f seconds\n" "$start" "$end" "$duration"
