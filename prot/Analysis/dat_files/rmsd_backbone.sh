#!/bin/bash
# RMSD BACKBONE (in ps)

# --- User parameters ---
dt=2.0               # Time spacing between frames in ps (change if different)
outfile="rmsd_bb_time.dat"

rm -f "$outfile"

# Initialize offset (in ps)
offset=0

# Loop through each RMSD file sequentially
for f in npt-1_rmsd.dat npt-2_rmsd.dat npt-3_rmsd.dat npt-4_rmsd.dat npt-5_rmsd.dat npt-6_rmsd.dat npt-7_rmsd.dat; do
  if [[ ! -f "$f" ]]; then
    echo "Warning: File $f not found, skipping."
    continue
  fi

  echo "Processing $f (offset = $offset ps)..."

  # Append adjusted times (add offset in ps)
  awk -v off="$offset" -v dt="$dt" 'NR>1 && NF>=2 {print $1*dt + off, $2}' "$f" >> "$outfile"

  # Compute last frame time in this file
  last=$(awk -v dt="$dt" 'NR>1 && NF>=2 {t=$1*dt} END{print t+0}' "$f")

  # Update offset for next file (in ps)
  offset=$(awk -v o="$offset" -v l="$last" 'BEGIN{print o+l}')
done

echo "Combined RMSD data written to: $outfile"
