set terminal pngcairo size 1000,600 enhanced font 'Arial,14'
set output 'rmsd_backbone.png'

set title "Backbone RMSD vs Time"
set xlabel "Time (ps)"
set ylabel "RMSD (Å)"
set grid

plot "rmsd_bb_time.dat" using 1:2 with lines lw 2 lc rgb "blue" title "Backbone RMSD"

set key top left
