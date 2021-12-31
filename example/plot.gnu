#!/usr/local/bin/gnuplot

set terminal postscript enhanced color 'Helvetica'
set view 0,0
set cbrange [3:10]
unset xtics
unset ytics
set output 'sampled.eps'
splot 'nn_pv_4_layer.dat' u 3:4:($5*100) w pm3d notitle, \
#'nn_sampled_4_layer.dat' u 3:4:($5*100) w p pointtype 7 ps 2 lc rgb "dark-green" notitle
