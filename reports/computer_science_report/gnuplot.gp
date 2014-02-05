set term post enh color solid
set output "./parallel.ps"
set title 'Speedup Vs Number of Processors'
set xlabel 'Number of Processors'
set ylabel 'Speedup factor'
set grid x y
set log x
set log y
plot "./parallel.out" u 1:3 w linespoints lt 1 lw 2 lc rgb "green" title "Speedup: Parallel Approach 2","./parallel1.out" u 1:3 w linespoints lt 1 lw 2 lc rgb "blue" title "Speedup: Parallel Approach 3","./parallel.out" u 1:1 w linespoints lt 1 lw 2 lc rgb "red" title "Linear"
