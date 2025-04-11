#!/bin/bash

# Record start time in seconds.nanoseconds
start=$(date +%s.%N)

# Move to script's directory (optional)
cd "$(dirname "$0")"

# Run FractureRB simulation
# ../build/FractureRB glass.bullet glass.csv \
gdb --args ../build/FractureRB glass.bullet glass.csv \
  -o _out/glass_ \
  -i 1e5 \
  -f 1e7 \
  -s 2



# Record end time
end=$(date +%s.%N)

# Calculate duration (in seconds)
duration=$(echo "$end - $start" | bc)

# Format duration into HH:MM:SS
duration_sec=$(printf "%.0f" "$duration")
hours=$((duration_sec / 3600))
minutes=$(( (duration_sec % 3600) / 60 ))
seconds=$((duration_sec % 60))

printf "\nStart Time : %s\nEnd Time   : %s\n" "$start" "$end"
printf "Duration   : %02d:%02d:%02d (%.2f seconds)\n" $hours $minutes $seconds "$duration"
