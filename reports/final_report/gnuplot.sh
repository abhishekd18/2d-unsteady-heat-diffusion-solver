echo "set term post enh color solid
set output \"./Timings_1.ps\"
set title 'Speedup Vs Number of Processors'
set xlabel 'Number of Processors'
set ylabel 'Speedup factor'
set grid x y
set log x
set log y
plot \"./Timings_1.dat\" u 1:3 w linespoints lt 1 lw 2 lc rgb \"green\" title \"Speedup\",\
\"./Timings_1.dat\" u 1:1 w linespoints lt 1 lw 2 lc rgb \"red\" title \"Linear\"">gnuplot.gp
gnuplot gnuplot.gp

echo "set term post enh color solid
set output \"./Timings_2.ps\"
set title 'Speedup Vs Number of Processors'
set xlabel 'Number of Processors'
set ylabel 'Speedup factor'
set grid x y
set log x
set log y
plot \"./Timings_2.dat\" u 1:3 w linespoints lt 1 lw 2 lc rgb \"green\" title \"Speedup\",\
\"./Timings_2.dat\" u 1:1 w linespoints lt 1 lw 2 lc rgb \"red\" title \"Linear\"">gnuplot.gp
gnuplot gnuplot.gp

