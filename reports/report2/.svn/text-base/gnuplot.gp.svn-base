set term post enh color solid
set output "./Timings.ps"
set title 'Speedup Vs Number of Processors'
set xlabel 'Number of Processors'
set ylabel 'Speedup factor'
set grid x y
plot "./Timings.dat" u 1:3 w l lt 1 lw 2 lc rgb "green" title "Speedup","./Timings.dat" u 1:1 w l lt 1 lw 2 lc rgb "red" title "Linear"
