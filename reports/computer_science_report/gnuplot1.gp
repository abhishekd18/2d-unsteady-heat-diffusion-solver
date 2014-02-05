set term post enh color solid
set output "./Performance.ps"
set title 'Speedup Vs Optimization levels'
set xlabel 'Optimization level'
set ylabel 'Speedup factor'
set xrange [0:3]
set xtics 1
set grid x y
plot "./performance.out" u 1:2 w linespoints lt 1 lw 2 lc rgb "green" title "Speedup"
